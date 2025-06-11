
library(survival)
library(BART)
library(glmnet)
library(mboost)
library(tidyverse)
library(AFTrees)
library(BARTTrial)
source("DrawIPCW.R")

## reading data
METABRIC <- read.csv('METABRIC_RNA_Mutation.csv', header = TRUE)[,c(1:41)]

######################
## data pre-processing
######################

## drop unwanted columns
METABRIC <- METABRIC %>%
  select(-c('patient_id', 'cancer_type', 'cancer_type_detailed', 'her2_status_measured_by_snp6', 'er_status',
            'death_from_cancer', 'X3.gene_classifier_subtype', 'oncotree_code', 'cellularity'))
## check number of missing in each column of the dataset
colSums(is.na(METABRIC))

## combine 4ER+ and 4ER- as 4 in integrative_cluster column
METABRIC$integrative_cluster[METABRIC$integrative_cluster %in% c("4ER+","4ER-")] <- "4"

METABRIC$tumor_subtype <- METABRIC$tumor_other_histologic_subtype
METABRIC$tumor_subtype[METABRIC$tumor_other_histologic_subtype != "Ductal/NST" & METABRIC$tumor_other_histologic_subtype != "Lobular"] <- "Mixed/Other"
METABRIC$tumor_subtype <- factor(METABRIC$tumor_subtype)

METABRIC <- METABRIC %>% select(-tumor_other_histologic_subtype)

## track the indices of dropped rows
dropped_indices <- list()
## drop tumor_stage = 0
table(METABRIC$tumor_stage)
## save the indices
dropped_indices$tumor_stage_0 <- which(METABRIC$tumor_stage == 0 & !is.na(METABRIC$tumor_stage))
METABRIC <- subset(METABRIC, !(tumor_stage == 0 & !is.na(tumor_stage)))

## replace tumor_stage = NA with cannot be evaluated
METABRIC$tumor_stage[is.na(METABRIC$tumor_stage)] <- 0
## make tumor_stage as a factor column
METABRIC$tumor_stage <- as.factor(METABRIC$tumor_stage)

## drop rows with missing tumor size
table(METABRIC$tumor_size)
## save the indices
dropped_indices$tumor_size_missing <- which(is.na(METABRIC$tumor_size))
METABRIC <- subset(METABRIC, !is.na(tumor_size))

## replace neoplasm = NA with 0
table(METABRIC$neoplasm_histologic_grade)
METABRIC$neoplasm_histologic_grade[is.na(METABRIC$neoplasm_histologic_grade)] <- 0
## make neoplasm as a factor column
METABRIC$neoplasm_histologic_grade <- as.factor(METABRIC$neoplasm_histologic_grade)

## replace mutation_count = NA with 0
table(METABRIC$mutation_count)
METABRIC$mutation_count[is.na(METABRIC$mutation_count)] <- 0

## drop rows empty in type_of_surgery
table(METABRIC$type_of_breast_surgery) ## 16 rows empty
## save the indices
dropped_indices$type_of_surgery_empty <- which(METABRIC$type_of_breast_surgery == '')
METABRIC <- subset(METABRIC, !(type_of_breast_surgery == ''))

## drop rows empty in er_status_measured_by_ihc
table(METABRIC$er_status_measured_by_ihc) ## 26 rows empty
## save the indices
dropped_indices$er_status_measured_by_ihc_empty <- which(METABRIC$er_status_measured_by_ihc == '')
METABRIC <- subset(METABRIC, !(er_status_measured_by_ihc == ''))

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


############################
## running different methods
############################

## split dataset into training and test set
set.seed(123)
DATA <- METABRIC[ , !(names(METABRIC) %in% c('overall_survival'))]
X.train <- model.matrix(overall_survival_months~.-1, data = DATA)


# Define the number of iterations and the proportion of data to be used for training
#tau <- 300
#sgrid <- seq(0, tau, by=100)
#tau <- 60
#sgrid <- seq(0, tau, by=10)
tau <- 120
sgrid <- seq(0, tau, by=20)
ndraws <- 2000
burnIn <- 500


# Create training and test sets using the selected index numbers
#train_idx <- sample(1:nrow(DATA), round(train_prop * nrow(DATA)), replace = FALSE)
#test_idx <- setdiff(1:nrow(DATA), train_idx)
#X.train <- DATA[train_idx, ]
#X.test <- DATA[test_idx, ]

Y <- METABRIC$overall_survival_months
delta <- METABRIC$overall_survival

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
bart_mod <- RMST_BART(Y, delta, x.train=X.train, Gweights=Gmat, tau=tau, k = 2,
                      ndpost=ndraws, nskip=burnIn, ntree=200)
bart_fitted <- bart_mod$yhat.train.mean
BART_CI <- t(apply(bart_mod$yhat.train, 1, function(x) quantile(x, probs=c(0.025, 0.975))))

GmatDep <- sqrt(2*eta_hat_star_dep)*GmatDeporig
bart_dep_mod <- RMST_BART(Y, delta, x.train=X.train, Gweights=GmatDep, tau=tau, k = 2.0,
                          ndpost=ndraws, nskip=burnIn, ntree = 200)

bart_dep_fitted <- bart_dep_mod$yhat.train.mean
BART_dep_CI <- t(apply(bart_dep_mod$yhat.train, 1, function(x) quantile(x, probs=c(0.025, 0.975))))

## plotting the first 10 repeated variables in one iteration - BART
## replace bart_mod with bart_dep_mod to get the dep plot
VarImp <- sort(colSums(bart_mod$varcount), decreasing=TRUE)[1:10]

library(ggplot2)
# Create a data frame with the numbers and names
VarImpDataF <- data.frame(numbers = c(VarImp)/ndraws,
                          names = c(names(VarImp)))
VarImpDataF <- VarImpDataF[order(VarImpDataF$numbers, decreasing = FALSE),]

library(ggrepel)
# Create the line plot

ggplot(VarImpDataF, aes(x = seq_along(numbers), y = numbers)) +
  geom_line(size = 1, color = "black") +
  geom_point(size = 3, color = "blue") +
  #geom_text(aes(label = names), hjust = -.5, vjust = -.9, size = 6, color = "blue") +
  #geom_text_repel(aes(label = names), size = 6, color = "blue") +
  geom_label_repel(aes(label = names), size = 10, color = "blue") +
  labs(x = "Variable", y = "Mean number of times used") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title=element_text(size=20),
    legend.text=element_text(size=20))

# Create bar charts

VarImp <- sort(colSums(bart_mod$varcount), decreasing=TRUE)[1:10]
VarImpdep <- sort(colSums(bart_dep_mod$varcount), decreasing=TRUE)[1:10]

library(ggplot2)
# Create a data frame with the numbers and names
VarImpDataF <- data.frame(numbers = c(VarImp)/ndraws,
                          names = c(names(VarImp)))
VarImpDataF <- VarImpDataF[order(VarImpDataF$numbers, decreasing = FALSE),]

VarImpDataFdep <- data.frame(numbers = c(VarImpdep)/ndraws,
                          names = c(names(VarImpdep)))
VarImpDataFdep <- VarImpDataFdep[order(VarImpDataFdep$numbers, decreasing = FALSE),]

# Indep - non-informative
## tau = 300
#VarImpDataF$names <- c("Tumor Size", "PR Positive", "Molecular Subtype Luminal A", "Integrative Cluster 7", "Nottingham Prognostic Index", 
#                       "Tumor Stage", "Cohort 2", "Chemotherapy", " Age at Diagnosis", "Cohort 3")

# Dep - informative

#VarImpDataFdep$names <- c("Cohort 2", "ATM", "Molecular Subtype Luminal A", "CHEK 2", "Tumor Size", "TP5 3", " Cohort 3", "NBN", "BRCA 1", "Age at Diagnosis")

# Indep - non-informative
## tau = 60
#VarImpDataF$names <- c("CHEK2", "Neoplasm Histologic Grade 3", "ATM", "Molecular SubtypeNormal", "Er Status Measured by Ihc Positve", 
#                       "Inferred Menopausal State Pre", "Molecular Subtype Her 2", "PR Positive", " NF1", "Molecular Subtype Lum A")

# Dep - informative

#VarImpDataFdep$names <- c("PALB 2", "PR Positive", "Hormone Therapy", "CHEK 2", "Histol Subtype Mixed/Other", "ATM", " Er Status Measured by Ihc Positve", 
#                          "Neoplasm Histologic Grade 3", "Nottingham Prognostic Index", "Chemotherapy")

## tau = 120
VarImpDataF$names <- c("Neoplasm Histologic Grade 3", "Tumor Size", "PTEN", "Cohort 3", "Tumor Stage 1", 
                       "Er Status Measured by Ihc Positve", "Molecular Subtype Lum A", "Chemotherapy", "PR Positive", "Nottingham Prognostic Index")

# Dep - informative

VarImpDataFdep$names <- c("Cohort 2", "Integrative Cluster 5", "CHEK 2", "Tumor Size", "Tumor Stage 1", "Molecular SubtypeLumA", "Chemotherapy", 
                          "PR Positive", "Er Status Measured by Ihc Positve", "Nottingham Prognostic Index")


# stacked bar chart

VarImpDataFdep$Censoring <- "Informative"
VarImpDataF$Censoring <- "Noninformative"
VarImpDataCombined <- rbind(VarImpDataFdep, VarImpDataF)

data <- VarImpDataCombined %>%
  group_by(Censoring) %>%
  arrange(desc(numbers), .by_group = TRUE) %>%
  mutate(index = row_number()) %>%
  ungroup()
data$Censoring <- factor(data$Censoring, levels = c("Noninformative", "Informative"))
ggplot(data, aes(y = rev(factor(index)), x = numbers, fill = Censoring)) +
  geom_bar(stat = "identity", position = "dodge", color = "gray60") +
  geom_text(aes(label = names), position = position_dodge(width = 0.85), color = "gray30", size = 6, fontface = "bold", vjust = 0.5, hjust = 1.1) +
  labs(x = "Mean Number of Times Used", y = "Variable") +
  scale_fill_manual(values = c("Informative" = "azure2", "Noninformative" = "gray80")) +
  coord_cartesian(xlim = c(0, 9.2)) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) 



# individual bar chart
#ggplot(VarImpDataF, aes(y = reorder(names, numbers), x = numbers)) +
#  geom_bar(stat = "identity", fill = 'gray80', width = 0.7) +
#  geom_text(aes(label = names), hjust = 1.1, size = 6, color = "white", fontface = "bold") +  
#  labs(y = "Variable", x = "Mean number of times used") +
#  theme_minimal() +
#  theme(
#    axis.title.x = element_text(size = 20),
#    axis.title.y = element_text(size = 20),
#    axis.text.x = element_text(size = 20),
#    axis.text.y = element_blank(),
#    legend.title = element_text(size = 20),
#    legend.text = element_text(size = 20)
#  )



## Confidence interval for each patient - plot
## replace bart_mod with bart_dep_mod to get the dep plot
means <- bart_mod$yhat.train.mean
df_summary <- data.frame(row = 1:ncol(bart_mod $yhat.train), mean = means)

df_summary$max <- apply(bart_mod $yhat.train, 2, function(x) quantile(x, probs=0.975))
df_summary$min <- apply(bart_mod $yhat.train, 2, function(x) quantile(x, probs=0.025))
## plot sample of lines
#sample_idx <- sample(1:nrow(df_summary), round(0.01 * nrow(df_summary)), replace = FALSE)
#sample_df_summary <- df_summary[sample_idx,]
#sample_long <- reshape2::melt(sample_df_summary, id.vars = "row")

# Sort the rows of df_summary by the mean column
df_summary_sorted <- df_summary[order(df_summary$mean),]


ggplot(df_summary_sorted, aes(y = 1:nrow(df_summary_sorted))) + 
  geom_segment(aes(x = min, xend = max, y = 1:nrow(df_summary_sorted), yend = 1:nrow(df_summary_sorted)), color = "gray") + 
  geom_point(aes(x = mean, y = 1:nrow(df_summary_sorted)), color = "blue", size = 2) +
  #scale_y_continuous(name = "Case number", breaks = c(0, 500, 1000, 1500)) + 
  # dep case
  #scale_x_continuous(name = "Months", limits = c(min(df_summary_sorted$min), max(df_summary_sorted$max))) + 
  ## to make sure they have same x-axis - indep case
  scale_x_continuous(name = "Months", limits = c(min(df_summary_sorted$min), 120)) + 
  labs(y = "Case Number") +
  theme_classic() + 
  theme(axis.text = element_text(size = 18),  
        axis.title = element_text(size = 20),
        #axis.text.y = element_blank(),
       #axis.line.y = element_blank(),
       #axis.ticks.y = element_blank(),
        panel.grid = element_blank(),         
        plot.margin = margin(5, 5, 5, 5))    

##########################################
## Partial Dependence plots
###########################################

Col_ParDep <- c('age_at_diagnosis', 'BRCA1', 'tumor_size', 'nottingham_prognostic_index')
Y <- METABRIC$overall_survival_months
delta <- METABRIC$overall_survival

## Setup Gmatrix and PD grid points before starting the loop:
Gmat <- sqrt(2*eta_hat_star)*Gmat_orig
ngrid <- 40
GridEnd <-  rbind(c(30, 95), c(-3, 3), c(0, 50), c(1, 6.5))
ff <- matrix(NA, nrow = ngrid, ncol = 3)
partial_results <- list()
for(i in 1:length(Col_ParDep)){
  #pp <- seq(min(METABRIC[,Col_ParDep[i]]),max(METABRIC[,Col_ParDep[i]]), length.out = ngrid)
  pp <- seq(GridEnd[i,1], GridEnd[i,2], length.out=ngrid)
  for(k in 1:ngrid){

    Xtmp <- X.train
    Xtmp[,Col_ParDep[i]] <- rep(pp[k], nrow(DATA))
    ## change name of bart_mod here
    bart_mod_tmps <- RMST_BART(Y, delta, x.train=X.train, Gweights=Gmat,
                               x.test=Xtmp, tau=tau, k = 2,
                               ndpost=ndraws, nskip=burnIn)
    ff[k,1] <- mean(bart_mod_tmps$yhat.test.mean)
    ff[k,2] <- pp[k]
    ff[k,3] <- Col_ParDep[i]
    print(c(i, k))
  }
  partial_results[[length(partial_results)+1]] <- ff
}

partial_results_all <- as.data.frame(do.call("rbind", partial_results))
names(partial_results_all) <- c('MeanPrediction', 'Value', 'index')
partial_results_all$xlabels <- rep("", nrow(partial_results_all))
partial_results_all$xlabels[partial_results_all$index =="age_at_diagnosis"] <- "Age at Diagnosis"
partial_results_all$xlabels[partial_results_all$index =="BRCA1"] <- "BRCA1"
partial_results_all$xlabels[partial_results_all$index =="tumor_size"] <- "Tumor Size"
partial_results_all$xlabels[partial_results_all$index =="nottingham_prognostic_index"] <- "Nottingham Prognostic Index"



library(ggplot2)
library(patchwork)

# Create a list to store the individual plots
plots <- list()
# Iterate over each category and create a plot
for (category in unique(partial_results_all[,3])) {
  subset_data <- subset(partial_results_all, partial_results_all[,3] == category)
  plot <- ggplot(subset_data, aes(x = as.numeric(Value), y = as.numeric(MeanPrediction), group = index)) +
    geom_line() +
    xlab(subset_data$xlabels[1])+
    scale_y_continuous(
      breaks = pretty_breaks(n = 4),  
      labels = scales::number_format(accuracy = 0.1)  # One decimal on y-axis
    ) +
    theme_light() +
    theme(axis.title = element_text(size = 22),  # Adjust the size of the axis titles
          axis.text = element_text(size = 20))

  # Set y-axis label for the leftmost plots
  if (category %in% unique(partial_results_all[,3])[c(1,3)]) {
    plot <- plot + ylab("Predicted RMST")
  } else {
    plot <- plot + theme(axis.title.y = element_blank())#,
                         #axis.text.y = element_blank())
    # Remove x-axis labels for the two plots on the right
  }
  plots[[category]] <- plot
}

## Need to change x-axis labels for these plots.
## CHEK2: mRNA expression level
# Combine the individual plots using patchwork
combined_plot <- wrap_plots(plots, ncol = 2)

# Display the combined plot
# Make sure the y-axis limits are visible in each panel.
combined_plot

#########################################################
### 3D plot for tumor_size and nottingham_prognostic_index
#########################################################
Col_ParDep <- c('tumor_size', 'nottingham_prognostic_index')
Gmat <- sqrt(2*eta_hat_star)*Gmat_orig
ngrid <- 40
GridEnd <-  rbind(c(0, 50), c(1, 6.5))

partial_results <- list()

npi_pp <- seq(GridEnd[1,1], GridEnd[1,2], length.out=ngrid)
PP <- expand.grid(seq(GridEnd[1,1], GridEnd[1,2], length.out=ngrid), seq(GridEnd[2,1], GridEnd[2,2], length.out=ngrid) )
lgrid <- nrow(PP)
ff <- matrix(NA, nrow = lgrid, ncol = 3)
for(k in 1:lgrid){
  
  Xtmp <- X.train
  ## npi is column index of NPI
  ## tumor is column index of tumor size
  Xtmp[,Col_ParDep[2]] <- rep(PP[k,2], nrow(DATA))
  Xtmp[,Col_ParDep[1]] <- rep(PP[k,1], nrow(DATA))
  ## change name of bart_mod here
  bart_mod_tmps <- RMST_BART(Y, delta, x.train=X.train, Gweights=Gmat,
                             x.test=Xtmp, tau=tau, k = 2,
                             ndpost=ndraws, nskip=burnIn)
  ff[k,1] <- mean(bart_mod_tmps$yhat.test.mean)
  ff[k,2] <- PP[k,1]
  ff[k,3] <- PP[k,2]
  print(c(i, k))
}




# Separate columns into variables
x <- ff[, 2]  
y <- ff[, 3]  
z <- ff[, 1]  

## While reading from CSV results

#x <- as.numeric(as.vector(unlist(ff[, 2]))) 
#y <- as.numeric(as.vector(unlist(ff[, 3]))) 
#z <- as.numeric(as.vector(unlist(ff[, 1])))  

# Generate a grid for the surface plot
x_vals <- sort(unique(x))
y_vals <- sort(unique(y))
RMST <- matrix(z, nrow = length(x_vals), ncol = length(y_vals), byrow = TRUE)

library(graphics)
persp(x_vals, y_vals, RMST, theta = 135, phi = 30, expand = 0.7, axes = TRUE, ticktype = "detailed",
      zlab = "Predicted RMST", 
      xlab = "Tumor Size", 
      ylab = "Nottingham Prognostic Index", 
      col = "gray95", 
      shade = 0.4, 
      cex.main = 1.5,  
      cex.lab = 1.8, 
      cex.axis = 1.5)    

## Another angel
persp(x_vals, y_vals, RMST, theta = 291, phi = -3, expand = 0.7, axes = TRUE, ticktype = "detailed",
      zlab = "Predicted RMST", 
      xlab = "Tumor Size", 
      ylab = "Nottingham Prognostic Index", 
      col = "gray95", 
      shade = 0.4, 
      cex.main = 1.5,  
      cex.lab = 1.8, 
      cex.axis = 1.5)   

## Using plotly package
#library(plotly)
#p <- plot_ly(x = ~x_vals, y = ~y_vals, z = ~RMST, colors = "Greys") %>%
#  add_surface() %>%
#  layout(scene = list(
#    xaxis = list(title = "Tumor Size"),
#    yaxis = list(title = "Nottingham Prognostic Index"),
#    zaxis = list(title = "Predicted RMST")
#  ))

#htmlwidgets::saveWidget(p, "3Dplot.html")

## 2D-heatmap
data <- read_csv('3D-Predicted-RMST.csv')
ggplot(data) +
  metR::geom_contour_fill(aes(z = `Tumor Size`, x = `Nottingham Prognostic Index`, y = `Predicted RMST`), 
                          na.fill = TRUE, bins = 10) +
  coord_cartesian(expand = FALSE) +
  labs(x = "Tumor Size", y = "Nottingham Prognostic Index")+
  scale_fill_gradient2(low = "white", high = "gray5", midpoint = 120, name = "Predicted RMST") + 
  #scale_fill_distiller(palette = "Spectral") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(linetype = 3, colour = "grey60"),
        axis.text = element_text(colour = 1, size = 20),
        axis.title = element_text(colour = 1, size = 25),
        legend.text = element_text(colour = 1, size = 20),
        legend.title = element_text(colour = 1, size = 25),
        legend.key.size = unit(1, "cm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom")

##############################
## Posterior mean comparison
##############################
PostMean_ind <- bart_mod$yhat.train.mean
PostMean_dep <-  bart_dep_mod$yhat.train.mean
PostMean_data <- data.frame('Independent' = PostMean_ind, 'Dependent' = PostMean_dep)
#write.csv(PostMean_data, 'PostMean_data_METABRIC.train.csv')

# Plot the data
## Might be useful to also include y=x line here.
## Also, maybe make x-axis limits and y-axis limits the same.
ggplot(PostMean_data, aes(x = Independent, y = Dependent)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0,  
              linetype = "dashed", color = "red", size = 1) +
  labs(x = "Noninformative censoring posterior means", y = "Informative censoring posterior means") +
  scale_x_continuous (expand = c(0.061, 0.061)) +
  coord_cartesian(xlim = c(0, 290), ylim = c(0, 290)) +
  theme_light() +
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24))   
