
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
METABRIC <- read.csv('METABRIC_RNA_Mutation.csv', header = TRUE)[,c(1:31)]

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

index <- sample(1:nrow(DATA), 400)
test.set <- DATA[index,]
train.set <- DATA[-index,]
Y <- METABRIC[-index,]$overall_survival_months
delta <- METABRIC[-index,]$overall_survival

Y.test <- METABRIC[index,]$overall_survival_months
mu.test <- log(Y.test)

## function we need for cox model
CoxExpectedSurv <- function(X, beta_val, H0fn, tau) {
  ## This function computes E( min(T_i, tau) |x_i) for
  ## a cox ph model
  mu.x <- colMeans(X)
  integrand <- function(time, xi, beta_val) {
    nu <- sum(xi*beta_val)
    ans <- exp(-H0fn(time)*exp(nu))
    return(ans)
  }
  nn <- nrow(X)
  fitted_vals <- rep(NA, nn)
  for(k in 1:nn) {
    II <- integrate(integrand, lower=0, upper=tau, xi=X[k,] - mu.x,
                    beta_val=beta_val, subdivisions=500L)
    fitted_vals[k] <- II$value
  }
  return(fitted_vals)
}


## BCART
sgrid <- seq(0, 4000, by=1)
tau <- 500
gam_alph <- 20
sigma <- 1.0

bcart_mod <- RMST_BCART(Y, delta, train.set, test.set,ndraws=500, tau=500, sigma.mu=1.2)
bcart_fitted <- pmin(rowMeans(bcart_mod$fitted.values.test), log(tau))

## Coxph
COXPH.mod <- coxph(Surv(Y, delta) ~ train.set)
coxhaz <- basehaz(COXPH.mod)
H0fn <- approxfun(c(0, coxhaz$time), c(0, coxhaz$hazard),
                  yright=max(coxhaz$hazard))
COXPH <- CoxExpectedSurv(X=test.set, beta_val=COXPH.mod$coefficients,
                         H0fn=H0fn, tau=500)
COXPH_fitted <- pmin(COXPH, log(tau))

## regularized coxph model
RCOXPH <-  glmnet(train.set, Surv(Y, delta), family = "cox", lambda = 1, alpha = 1)
RCOXPH_fitted <- pmin(c(predict(RCOXPH, test.set, type = 'response')), log(tau))

## survival boosting
SBOOST <- glmboost(Surv(Y, delta)~train.set, family = Gehan(), control = boost_control(mstop = 300))
## we have warnings in predict function: 'newdata' had 2000 rows but variables found have 250 rows
SBOOST_fitted <- pmin(c(predict(SBOOST, newdata = data.frame(test.set))), log(tau))

## AFT
AFT <- survreg(Surv(Y, delta) ~ train.set)
## we have warnings in predict function: 'newdata' had 2000 rows but variables found have 250 rows
AFT_fitted <- pmin(predict(AFT, newdata = data.frame(test.set), type = 'response'), log(tau))
#AFT_fitted <- pmin(AFT$linear.predictors, log(tau))

## AFT BART
AFT_BART <- abart(train.set, Y, x.test = test.set, delta)
AFT_BART_fitted <-  pmin(AFT_BART$yhat.test.mean, log(tau))

## AFT null
AFT_null <- survreg(Surv(Y, delta) ~ 1)
AFT_null_fitted <- pmin(predict(AFT_null, newdata = data.frame(test.set), type = 'response'), log(tau))

rmse_bcart <- sqrt(mean((bcart_fitted - mu.test)*(bcart_fitted - mu.test)))
rmse_coxph <- sqrt(mean((COXPH_fitted - mu.test)*(COXPH_fitted - mu.test)))
rmse_rcoxph <- sqrt(mean((RCOXPH_fitted - mu.test)*(RCOXPH_fitted - mu.test)))
rmse_sboost <- sqrt(mean((SBOOST_fitted - mu.test)*(SBOOST_fitted - mu.test)))
rmse_aft <- sqrt(mean((AFT_fitted - mu.test)*(AFT_fitted - mu.test)))
rmse_aft_bart <- sqrt(mean((AFT_BART_fitted - mu.test)*(AFT_BART_fitted - mu.test)))
rmse_aft_null <- sqrt(mean((AFT_null_fitted - mu.test)*(AFT_null_fitted - mu.test)))


