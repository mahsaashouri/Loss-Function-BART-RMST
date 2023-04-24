
library(survival)
library(BART)
library(glmnet)
library(mboost)
library(tidyverse)

setwd("~/Documents/LossFunctionBART/Loss-Function--BART")
source("A_matrix.R")
source("ChangeMove.R")
source("D_probabilities.R")
source("GPDraw.R")
source("GrowMove.R")
source("LogLik.R")
source("PruneMove.R")
#source("PruneMove_V2.R")
source("RMST_BCART.R")
source("RMST_BART.R")
source("RMST_MHRatio.R")
source("Recursive_A_matrix.R")
source("Transition_Prob.R")
source("fitted_value.R")
source("prior_conditional_on_D_V2.R")
source("tree-configuration.R")
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
n_iterations <- 3
train_prop <- 0.7
sgrid <- seq(0, 4000, by=1)
tau <- 500
gam_alph <- 20
sigma <- 1.0
ndraws <- 500

bcart_fitted <- bart_fitted <- AFT_BART_fitted <- list()
# Loop through the iterations
for (i in 1:n_iterations) {
  
  # Randomly select index numbers for the training and test sets
  train_idx <- sample(1:nrow(DATA), round(train_prop * nrow(DATA)), replace = FALSE)
  test_idx <- setdiff(1:nrow(DATA), train_idx)
  
  # Create training and test sets using the selected index numbers
  train.set <- DATA[train_idx, ]
  test.set <- DATA[test_idx, ]
  
  Y <- METABRIC[train_idx,]$overall_survival_months
  delta <- METABRIC[train_idx,]$overall_survival
  
  Y.test <- METABRIC[test_idx,]$overall_survival_months
  #mu.test <- log(Y.test)
  mu.test <- (Y.test/gam_alph)*pgamma(tau, shape = gam_alph+1, rate = Y.test) + 
    tau*pgamma(tau, shape = gam_alph, rate = Y.test, lower.tail = FALSE)
  
  
  ## BCART
  bcart_mod <- RMST_BCART(Y, delta, train.set, test.set,ndraws=ndraws, tau=tau)
  #bcart_fitted[[i]] <- rowMeans(bcart_mod$fitted.values.test)
  bcart_fitted[[i]] <- bcart_mod
  
  ## BART
  bart_mod <- RMST_BART(Y, delta, train.set, test.set, ndraws=ndraws, tau=tau)
  #bart_fitted[[i]] <- rowMeans(bart_mod$fitted.values.test)
  bart_fitted[[i]] <- bart_mod
  
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
  AFT_BART_fitted[[i]] <- AFT_fit_reps
}

## plotting the first 10 repeated variables in one iteration - BART
VarImp <- tail(sort(colSums(bart_fitted[[1]]$split.vars)),10)

library(ggplot2)
# Create a data frame with the numbers and names
VarImpDataF <- data.frame(numbers = c(VarImp),
                 names = c(names(VarImp)))
# Create the plot
ggplot(VarImpDataF, aes(x = seq_along(numbers), y = numbers)) +
  geom_line(size = 1, color = "black") +
  geom_point(size = 3, color = "blue") +
  geom_text(aes(label = names), hjust = 1.1, vjust = -.5, size = 6, color = "blue") +
  labs(x = "Variable", y = "Count") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.title=element_text(size=15), 
    legend.text=element_text(size=15))

