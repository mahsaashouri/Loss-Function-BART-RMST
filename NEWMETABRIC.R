
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


## drop tumor_stage = 0
table(METABRIC$tumor_stage)
METABRIC <- subset(METABRIC, !(tumor_stage == 0 & !is.na(tumor_stage)))

## replace tumor_stage = NA with cannot be evaluated
METABRIC$tumor_stage[is.na(METABRIC$tumor_stage)] <- 0
## make tumor_stage as a factor column
METABRIC$tumor_stage <- as.factor(METABRIC$tumor_stage)

## drop rows with missing tumor size
table(METABRIC$tumor_size)
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
METABRIC <- subset(METABRIC, !(type_of_breast_surgery == ''))

## drop rows empty in er_status_measured_by_ihc
table(METABRIC$er_status_measured_by_ihc) ## 25 rows empty
METABRIC <- subset(METABRIC, !(er_status_measured_by_ihc == ''))

## make cohort as a factor column
METABRIC$cohort <- as.factor(METABRIC$cohort)

## flip 1 and 0 valves in overall survival
table(METABRIC$overall_survival)
METABRIC$overall_survival <- ifelse(METABRIC$overall_survival == 1, 0, ifelse(METABRIC$overall_survival == 0, 1,
                                                                              METABRIC$overall_survival))
############################
## running different methods
############################

## split dataset into training and test set
set.seed(123)
## 400 rows as test set (~20%)
DATA <- METABRIC[ , !(names(METABRIC) %in% c('overall_survival'))]
DATA <- model.matrix(overall_survival_months~.-1, data = DATA)


# Define the number of iterations and the proportion of data to be used for training
n_iterations <- 1
train_prop <- 0.7
tau <- 300
sgrid <- seq(0, tau, by=100)
gam_alph <- 20
sigma <- 1.0
ndraws <- 2000
burnIn <- 100

bcart_fitted <- bart_fitted_ind <- bart_fitted_dep <- AFT_BART_fitted <- list()
# Loop through the iterations
for (j in 1:n_iterations) {
  
  # Randomly select index numbers for the training and test sets
  train_idx <- sample(1:nrow(DATA), round(train_prop * nrow(DATA)), replace = FALSE)
  test_idx <- setdiff(1:nrow(DATA), train_idx)
  
  # Create training and test sets using the selected index numbers
  train.set <- DATA[train_idx, ]
  test.set <- DATA[test_idx, ]
  
  Y <- METABRIC[train_idx,]$overall_survival_months
  delta <- METABRIC[train_idx,]$overall_survival
  
  Y.test <- METABRIC[test_idx,]$overall_survival_months
  rate.test <- 1 + Y.test
  mu.test <- Y.test*pgamma(tau, shape = rate.test+1, rate = rate.test) +
    tau*pgamma(tau, shape = rate.test, rate = rate.test, lower.tail = FALSE)
  
  delta.test <- METABRIC[test_idx,]$overall_survival
  
  ##  RMST BART
  delta_alpha <- 1
  kappa0 <- 1
  U_tau <- pmin(Y[delta==1], tau)
  Gmat <- matrix(1, nrow=ndraws + burnIn + 1, ncol=length(U_tau))
  for(k in 1:(ndraws + burnIn + 1)) {
    Gmat[k,] <- DrawIPCW(U=Y, delta=delta, Utau=U_tau, sgrid=sgrid,
                         kappa0=kappa0, delta_alpha=delta_alpha)
  }
  Gmat <- 1/sqrt(Gmat)
  
  bart_mod <- RMST_BART(Y, delta, train.set, Gweights=Gmat,
                        x.test=test.set, tau=tau, k = 2,
                        ndpost=ndraws, nskip=burnIn, ntree = 200)
  bart_mod$yhat.train <- pmax(bart_mod$yhat.train, 0)
  bart_mod$yhat.test <- pmax(bart_mod$yhat.test, 0)
  bart_mod$yhat.train.mean <- colMeans(bart_mod$yhat.train)
  bart_mod$yhat.test.mean <- colMeans(bart_mod$yhat.test)
  
  bart_fitted_ind[[j]] <- bart_mod
  BART_CI <- t(apply(bart_mod$yhat.test, 1, function(x) quantile(x, probs=c(0.025, 0.975))))
  
  
  cens_bart <- AFTrees(x.train=train.set, y.train=Y, status=1-delta,
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
  GmatDep <- 1/sqrt(GmatDep)
  bart_dep_mod <- RMST_BART(Y, delta, train.set, Gweights=GmatDep,
                            x.test=test.set, tau=tau, k = 2.0,
                            ndpost=ndraws, nskip=burnIn, ntree = 200)
  bart_dep_mod$yhat.train <- pmax(bart_dep_mod$yhat.train, 0)
  bart_dep_mod$yhat.test <- pmax(bart_dep_mod$yhat.test, 0)
  bart_dep_mod$yhat.train.mean <- colMeans(bart_dep_mod$yhat.train)
  bart_dep_mod$yhat.test.mean <- colMeans(bart_dep_mod$yhat.test)
  
  bart_fitted_dep[[j]] <- bart_dep_mod
  BART_dep_CI <- t(apply(bart_dep_mod$yhat.test, 1, function(x) quantile(x, probs=c(0.025, 0.975))))
  
  ## BCART
  bcart_mod <- RMST_BART(Y, delta, train.set, Gweights=Gmat,
                         x.test=test.set, tau=tau, k = 2, ntree=1L,
                         ndpost=ndraws, nskip=burnIn)
  
  bcart_fitted[[j]] <- bcart_mod
  
  
  #### 3. AFT_BART model
  AFT_BART <- abart(train.set, Y, delta, x.test=test.set)
  ndraw_abart <- nrow(AFT_BART$yhat.test)
  AFT_fit_reps <- matrix(NA, nrow=ndraw_abart, ncol=nrow(test.set))
  for(k in 1:ndraw_abart) {
    aft_bart_mu <- AFT_BART$yhat.test[k,]
    aft_bart_sig <- AFT_BART$sigma[k]
    aft_bart_sigsq <- aft_bart_sig*aft_bart_sig
    gt_prob <- pnorm((log(tau) - aft_bart_mu)/aft_bart_sig, lower.tail=FALSE)
    lt_prob <- pnorm((log(tau) - aft_bart_sigsq - aft_bart_mu)/aft_bart_sig)
    
    AFT_fit_reps[k,] <- exp(aft_bart_sigsq/2 + aft_bart_mu)*lt_prob + tau*gt_prob
  }
  #AFT_BART_fitted[[i]] <- colMeans(AFT_fit_reps)
  AFT_BART_fitted[[j]] <- AFT_fit_reps
}


## plotting the first 10 repeated variables in one iteration - BART
VarImp <- tail(sort(colSums(bart_fitted_ind[[1]]$varcount)),10)

library(ggplot2)
# Create a data frame with the numbers and names
VarImpDataF <- data.frame(numbers = c(VarImp)/ndraws,
                          names = c(names(VarImp)))
library(ggrepel)
# Create the plot
ggplot(VarImpDataF, aes(x = seq_along(numbers), y = numbers)) +
  geom_line(size = 1, color = "black") +
  geom_point(size = 3, color = "blue") +
  #geom_text(aes(label = names), hjust = -.5, vjust = -.9, size = 6, color = "blue") +
  #geom_text_repel(aes(label = names), size = 6, color = "blue") +
  geom_label_repel(aes(label = names), size = 6, color = "blue") +
  labs(x = "Variable", y = "Count") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.title=element_text(size=15),
    legend.text=element_text(size=15))

## Average RMST

bart.r <- bart_fitted_ind[[1]]$yhat.test.mean
sqrt(mean((bart.r - mu.test)*(bart.r - mu.test)))
bcart.r <- bcart_fitted[[1]]$yhat.test.mean
sqrt(mean((bcart.r - mu.test)*(bcart.r - mu.test)))
AFT.r <- colMeans(AFT_fit_reps)
sqrt(mean((AFT.r - mu.test)*(AFT.r - mu.test)))


## Confidence interval for each case - plot
means <- bart_fitted_ind[[1]]$yhat.test.mean
ses <- matrixStats::colSds(bart_fitted_ind[[1]]$yhat.test, na.rm=TRUE)
df_summary <- data.frame(row = 1:ncol(bart_fitted_ind[[1]]$yhat.test), mean = means, se = ses)

df_summary$max <- df_summary$mean + 1.96 * df_summary$se
df_summary$min <- df_summary$mean - 1.96 * df_summary$se
df_summary <- subset(df_summary, select = -c(se))
## plot sample of lines
sample_idx <- sample(1:nrow(df_summary), round(0.01 * nrow(df_summary)), replace = FALSE)
sample_df_summary <- df_summary[sample_idx,]
#sample_long <- reshape2::melt(sample_df_summary, id.vars = "row")

# Sort the rows of df_summary by the mean column
df_summary_sorted <- df_summary[order(df_summary$mean),]

# create the plot
plot(1, type='n', xlim=c(min(df_summary_sorted$min), max(df_summary_sorted$max)),
     ylim=c(0.5, nrow(df_summary_sorted)+0.5), ylab='Case number', xlab='Test set credible intervals')

# loop through the rows and draw a line for each
for (i in 1:nrow(df_summary)) {
  segments(df_summary_sorted$min[i], i, df_summary_sorted$max[i], i)
  points(df_summary_sorted$mean[i], i, pch=19, col='blue')
}

## Partial Dependence plots

Col_ParDep <- c('atm', 'palb2', 'pten','nottingham_prognostic_index')
Y <- METABRIC$overall_survival_months
delta <- METABRIC$overall_survival

##  RMST BART
U_tau <- pmin(Y[delta==1], tau)
Gmat <- matrix(1, nrow=ndraws + burnIn + 1, ncol=length(U_tau))
for(k in 1:(ndraws + burnIn + 1)) {
  Gmat[k,] <- DrawIPCW(U=Y, delta=delta, Utau=U_tau, sgrid=sgrid,
                       kappa0=kappa0, delta_alpha=delta_alpha)
}
Gmat <- 1/sqrt(Gmat)

DATA <- METABRIC[ , !(names(METABRIC) %in% c('overall_survival'))]
DATA <- model.matrix(overall_survival_months~.-1, data = DATA)
ff <- matrix(NA, nrow = 100, ncol = 3)
partial_results <- list()
for(i in 1:length(Col_ParDep)){
  pp <- seq(min(METABRIC[,Col_ParDep[i]]),max(METABRIC[,Col_ParDep[i]]), length.out = 100)
  for(k in 1:100){
    
    Xtmp <- DATA
    Xtmp[,Col_ParDep[i]] <- rep(pp[k], nrow(DATA))
    Y <- METABRIC$overall_survival_months
    delta <- METABRIC$overall_survival
    bart_mod <- RMST_BART(Y, delta, x.train=DATA, Gweights=Gmat, x.test=Xtmp, tau=tau, k = 2,
                          ndpost=ndraws, nskip=burnIn)
    ff[k,1] <- mean(colMeans(bart_mod$yhat.test))
    ff[k,2] <- pp[k]
    ff[k,3] <- Col_ParDep[i]
  }
  partial_results[[length(partial_results)+1]] <- ff
}

partial_results_all <- do.call("rbind", partial_results)
colnames(partial_results_all) <- c('MeanPrediction', 'Value', 'index')

library(ggplot2)
library(patchwork)

# Create a list to store the individual plots
plots <- list()

# Iterate over each category and create a plot
for (category in unique(partial_results_all[,3])) {
  subset_data <- as.data.frame(subset(partial_results_all, partial_results_all[,3] == category))
  plot <- ggplot(subset_data, aes(x = as.numeric(Value), y = as.numeric(MeanPrediction), group = index)) +
    geom_line() +
    xlab(subset_data$index)+
    theme_light() +
    theme(axis.title = element_text(size = 22),  # Adjust the size of the axis titles
          axis.text = element_text(size = 20))
  
  # Set y-axis label for the leftmost plots
  if (category %in% unique(partial_results_all[,3])[c(1,3)]) {
    plot <- plot + ylab("Mean Prediction") #+ theme(axis.title.y = element_blank())
  } else {
    plot <- plot + theme(axis.title.y = element_blank(),
                         axis.text.y = element_blank())
    # Remove x-axis labels for the two plots on the right
    #plot <- plot + scale_x_continuous(labels = NULL)
  }
  
  plots[[category]] <- plot
}

# Combine the individual plots using patchwork
combined_plot <- wrap_plots(plots, ncol = 2)

# Display the combined plot
combined_plot


## Posterior mean comparison

PostMean_ind <- c(bart_fitted_ind[[1]]$yhat.train.mean, bart_fitted_ind[[1]]$yhat.test.mean)
PostMean_dep <-  c(bart_fitted_dep[[1]]$yhat.train.mean, bart_fitted_dep[[1]]$yhat.test.mean)
PostMean_data <- data.frame('Independent' = PostMean_ind, 'Dependent' = PostMean_dep)
#write.csv(PostMean_data, 'PostMean_data_METABRIC.train.csv')

# Plot the data

ggplot(PostMean_data, aes(x = Independent, y = Dependent)) +
  geom_point() +
  labs(x = "Independent", y = "Dependent", title = "") +
  theme_light() +
  theme(axis.title = element_text(size = 22),  # Adjust the size of the axis titles
        axis.text = element_text(size = 20))   # Adjust the size of the axis labels


