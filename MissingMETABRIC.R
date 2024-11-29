
## Handling the missing values using mice package

library(survival)
library(mice)
library(BART)
library(glmnet)
library(mboost)
library(tidyverse)
library(AFTrees)
library(BARTTrial)
source("DrawIPCW.R")


METABRIC <- read.csv('METABRIC_RNA_Mutation.csv', header = TRUE)[,c(1:41)]

######################
## data pre-processing
######################

## drop unwanted columns
METABRIC <- METABRIC %>%
  select(-c('patient_id', 'cancer_type', 'cancer_type_detailed', 'her2_status_measured_by_snp6', 'er_status',
            'death_from_cancer', 'X3.gene_classifier_subtype', 'oncotree_code', 'cellularity'))

## combine 4ER+ and 4ER- as 4 in integrative_cluster column
METABRIC$integrative_cluster[METABRIC$integrative_cluster %in% c("4ER+","4ER-")] <- "4"

METABRIC$tumor_subtype <- METABRIC$tumor_other_histologic_subtype
METABRIC$tumor_subtype[METABRIC$tumor_other_histologic_subtype != "Ductal/NST" & METABRIC$tumor_other_histologic_subtype != "Lobular"] <- "Mixed/Other"
METABRIC$tumor_subtype <- factor(METABRIC$tumor_subtype)

METABRIC <- METABRIC %>% select(-tumor_other_histologic_subtype)

METABRIC$tumor_stage <- as.factor(METABRIC$tumor_stage)


METABRIC$neoplasm_histologic_grade <- as.factor(METABRIC$neoplasm_histologic_grade)

## Replace '' with NA
table(METABRIC$type_of_breast_surgery) ## 16 rows empty
#METABRIC$type_of_breast_surgery[METABRIC$type_of_breast_surgery==''] <- NA

## drop rows empty in er_status_measured_by_ihc
table(METABRIC$er_status_measured_by_ihc) ## 25 rows empty
#METABRIC$er_status_measured_by_ihc[METABRIC$er_status_measured_by_ihc==''] <- NA

## make cohort as a factor column
METABRIC$cohort <- as.factor(METABRIC$cohort)

## flip 1 and 0 valves in overall survival
table(METABRIC$overall_survival)
METABRIC$overall_survival <- ifelse(METABRIC$overall_survival == 1, 0, ifelse(METABRIC$overall_survival == 0, 1,
                                                                              METABRIC$overall_survival))

names(METABRIC)[which(names(METABRIC)=="pam50_._claudin.low_subtype")] <- "molecular_subtype"
names(METABRIC)[which(names(METABRIC)=="tumor_subtype")] <- "histol_subtype"
names(METABRIC)[which(names(METABRIC)=="pr_status")] <- "PR_"
names(METABRIC)[which(names(METABRIC)=="brca1")] <- "BRCA1"
names(METABRIC)[which(names(METABRIC)=="brca2")] <- "BRCA2"
names(METABRIC)[which(names(METABRIC)=="palb2")] <- "PALB2"
names(METABRIC)[which(names(METABRIC)=="pten")] <- "PTEN"
names(METABRIC)[which(names(METABRIC)=="tp53")] <- "TP53"
names(METABRIC)[which(names(METABRIC)=="atm")] <- "ATM"
names(METABRIC)[which(names(METABRIC)=="cdh1")] <- "CDH1"
names(METABRIC)[which(names(METABRIC)=="chek2")] <- "CHEK2"
names(METABRIC)[which(names(METABRIC)=="nbn")] <- "NBN"
names(METABRIC)[which(names(METABRIC)=="nf1")] <- "NF1"

## only one value of outcome is zero - I dropped that
which(METABRIC$overall_survival_months == 0)
METABRIC <- METABRIC[METABRIC$overall_survival_months != 0, ]



imputed.metabric <- mice(METABRIC, seed=1378)
completed.metabric <- mice::complete(imputed.metabric, action="long")
table(completed.metabric$.imp)

METABRICList <- list()
for(k in 1:5) {
  METABRICList[[k]] <- subset(completed.metabric, .imp==k)
  METABRICList[[k]] <- METABRICList[[k]] %>% select(-.imp, -.id)
}

#sapply(METABRICM, function(x) sum(is.na(x)))


bart_mod <- bart_dep_mod <- list()
for(i in 1:length(METABRICList)){
  ############################
  ## running different methods
  ############################
  METABRICM <- METABRICList[[i]]
  set.seed(123)
  DATA <- METABRICM[ , !(names(METABRICM) %in% c('overall_survival'))]
  X.train <- model.matrix(overall_survival_months~.-1, data = DATA)
  
  
  # Define the number of iterations and the proportion of data to be used for training
  tau <- 300
  sgrid <- seq(0, tau, by=100)
  ndraws <- 2000
  burnIn <- 500
  Y <- METABRICM$overall_survival_months
  delta <- METABRICM$overall_survival
  
  ###########################################
  ## Draw Censoring weights for both
  ## the independent and dependent cases
  ###########################################
  delta_alpha <- 1
  kappa0 <- 1
  U_tau <- pmin(Y[delta==1], tau)
  Gmat <- matrix(1, nrow=ndraws + burnIn + 1, ncol=length(U_tau))
  for(k in 1:(ndraws + burnIn + 1)) {
    Gmat[k,] <- DrawIPCW(U=Y, delta=delta, Utau=U_tau, sgrid=sgrid,
                         kappa0=kappa0, delta_alpha=delta_alpha)
  }
  Gmat_orig <- 1/sqrt(Gmat)
  
  ## Added small value (1e-100) to the outcome - one of the outcomes is zero 
  cens_bart <- AFTrees(x.train=X.train, y.train=Y, status=1-delta,
                       ndpost=ndraws + burnIn, verbose=FALSE)
  Mucens_draws <- cens_bart$m.train
  GmatDep <- matrix(1, nrow=ndraws + burnIn + 1, ncol=length(U_tau))
  for(k in 1:(ndraws + burnIn)) {
    for(h in 1:length(U_tau)) {
      log.time.points <- log(U_tau[h])
      AA <- (log.time.points - cens_bart$locations[k,] - Mucens_draws[k,h])/cens_bart$sigma[k]
      Cprob <- sum(pnorm(AA, lower.tail=FALSE)*cens_bart$mix.prop[k,])
      GmatDep[k,h] <- 1/Cprob
    }
  }
  GmatDeporig <- 1/sqrt(GmatDep)
  
  #########################################
  ## Get initial estimate of eta_hat
  AFT_try <- survreg(Surv(exp(Y), delta) ~ X.train)
  eta_hat <- AFT_try$scale*AFT_try$scale
  ############################################
  
  ########################################################
  ## Now, do 5-fold cross-validation to get the best
  ## value of eta_hat
  nfolds <- 5
  eta_hat_vals <- c(0.1*eta_hat, 0.25*eta_hat, 0.5*eta_hat, 0.75*eta_hat, eta_hat)
  CVscore <- matrix(NA, nrow=nfolds, ncol=length(eta_hat_vals))
  X.train.obs <- X.train[delta==1,]
  Y.train.obs <- Y[delta==1]
  ntrain <- nrow(X.train.obs)
  fold_memb <- sample(1:nfolds, size=ntrain, replace=TRUE)
  CVscore <- CVscoreDep <- matrix(NA, nrow=nfolds, ncol=length(eta_hat_vals))
  cens_dist <- survfit(Surv(Y, 1-delta) ~ 1)
  GKM <- stepfun(c(0, cens_dist$time), c(1, cens_dist$surv, min(cens_dist$surv)))
  for(u in 1:length(eta_hat_vals)) {
    Gmat <- sqrt(2*eta_hat_vals[u])*Gmat_orig
    GmatDep <- sqrt(2*eta_hat_vals[u])*GmatDeporig
    for(k in 1:nfolds) {
      Y.train_tmp <- Y.train.obs[fold_memb!=k]
      delta.train_tmp <- rep(1, sum(fold_memb!=k))
      X.train_tmp <- X.train.obs[fold_memb!=k,]
      
      Y.test_tmp <- Y.train.obs[fold_memb==k]
      delta.test_tmp <- rep(1, sum(fold_memb==k))
      X.test_tmp <- X.train.obs[fold_memb==k,]
      
      Gmat_train_tmp <- Gmat[,fold_memb!=k]
      Gmat_test_tmp <- Gmat[,fold_memb==k]
      
      GmatDep_train_tmp <- GmatDep[,fold_memb!=k]
      GmatDep_test_tmp <- GmatDep[,fold_memb==k]
      
      ww_dep <- colMeans(1/(GmatDeporig[,fold_memb==k]^2))
      ww <- 1/GKM(Y.test_tmp)
      
      bartmod_tmp <- RMST_BART(Y.train_tmp, delta.train_tmp, X.train_tmp, Gweights=Gmat_train_tmp,
                               x.test=X.test_tmp, tau=tau, k = 2, ndpost=ndraws, nskip=burnIn)
      
      bartmod_dep_tmp <- RMST_BART(Y.train_tmp, delta.train_tmp, X.train_tmp, Gweights=GmatDep_train_tmp,
                                   x.test=X.test_tmp, tau=tau, k = 2,
                                   ndpost=ndraws, nskip=burnIn)
      
      yhat <- bartmod_tmp$yhat.test.mean
      yhat_dep <- bartmod_dep_tmp$yhat.test.mean
      
      CVscore[k, u] <- mean(ww*((Y.test_tmp - yhat)*(Y.test_tmp - yhat)))
      CVscoreDep[k,u] <- mean(ww_dep*((Y.test_tmp - yhat_dep)*(Y.test_tmp - yhat_dep)))
    }
  }
  CVfinal <- colMeans(CVscore)
  CVfinalDep <- colMeans(CVscoreDep)
  eta_hat_star <- eta_hat_vals[which.min(CVfinal)]
  ## There is quite a bit of variability in eta_hat_star_dep depending on the
  ## randomness of the folds in cross-validation.
  ## Setting eta_hat_star_dep = eta_hat_vals[1] is a more stable and interpretable
  ## choice since eta_hat_star_dep equals this on a number of rounds of cross-validation
  ## and provides more regularization than some of the extreme estimates of eta_hat_star_dep
  ## that occurs for some of the train-test splits in cross-validation
  eta_hat_star_dep <- eta_hat_vals[1]
  #######################
  
  ############################
  ## Now using best CV parameters with RMST-BART
  #############################
  
  
  Gmat <- sqrt(2*eta_hat_star)*Gmat_orig
  bart_mod[[i]] <- RMST_BART(Y, delta, x.train=X.train, Gweights=Gmat, tau=tau, k = 2,
                        ndpost=ndraws, nskip=burnIn, ntree=200)
  
  GmatDep <- sqrt(2*eta_hat_star_dep)*GmatDeporig
  bart_dep_mod[[i]] <- RMST_BART(Y, delta, x.train=X.train, Gweights=GmatDep, tau=tau, k = 2.0,
                            ndpost=ndraws, nskip=burnIn, ntree = 200)
}


##############################
## Posterior mean comparison
##############################
## pointwise average - missing values imputation
vectors <- list(bart_mod[[1]]$yhat.train.mean, 
                bart_mod[[2]]$yhat.train.mean, 
                bart_mod[[3]]$yhat.train.mean, 
                bart_mod[[4]]$yhat.train.mean, 
                bart_mod[[5]]$yhat.train.mean)
PostMean_ind_missing <- Reduce("+", vectors) / length(vectors)

all_dropped_indices <- as.vector(unlist(dropped_indices)) ## from NEWMETABRIC file

## since I already remove the row 171 (overall_survival = 0) here before running the missing evaluation
all_dropped_indices_modified <- all_dropped_indices[all_dropped_indices != 170]
# Decrements values greater than 170 by 1
all_dropped_indices_modified <- ifelse(all_dropped_indices_modified > 170, 
                                       all_dropped_indices_modified - 1, 
                                       all_dropped_indices_modified)

# Drop the corresponding elements from PostMean_ind_missing
PostMean_ind_missing_cleaned <- PostMean_ind_missing[-all_dropped_indices_modified]

vectors_dep <- list(bart_dep_mod[[1]]$yhat.train.mean, 
                bart_dep_mod[[2]]$yhat.train.mean, 
                bart_dep_mod[[3]]$yhat.train.mean, 
                bart_dep_mod[[4]]$yhat.train.mean, 
                bart_dep_mod[[5]]$yhat.train.mean)
PostMean_dep_missing <- Reduce("+", vectors) / length(vectors_dep)
# Drop the corresponding elements from PostMean_dep_missing
PostMean_dep_missing_cleaned <-PostMean_dep_missing[-all_dropped_indices_modified]

## computation by dropping missing rows
PostMean_ind <- bart_mod_n$yhat.train.mean ## from NEWMETABRIC file
PostMean_dep <-  bart_dep_mod_n$yhat.train.mean ## from NEWMETABRIC file


library(ggplot2)

# Combine the data into a single data frame for easier plotting
data <- data.frame(
  x = c(PostMean_ind_missing_cleaned, PostMean_dep_missing_cleaned),
  y = c(PostMean_ind, PostMean_dep),
  type = rep(c("Noninformative Censoring", "Informative Censoring"), 
             c(length(PostMean_ind_missing_cleaned), length(PostMean_dep_missing_cleaned)))
)

data$type <- factor(data$type, levels = c("Noninformative Censoring", "Informative Censoring"))
# Create the plot
ggplot(data) +
  geom_point(aes(x = x, y = y), alpha = 0.7) +                
  geom_abline(slope = 1, intercept = 0,  
              linetype = "dashed", color = "red", size = 1) +
  facet_wrap(~type, scales = "free") +    
  labs(
    x = "Missing Values Imputed",
    y = "Missing Values Dropped"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),          
    axis.title = element_text(size = 25),    
    axis.text = element_text(size = 20),     
    strip.text = element_text(size = 21),    
    legend.text = element_text(size = 20),   
    plot.title = element_text(size = 22)     
  )                       

