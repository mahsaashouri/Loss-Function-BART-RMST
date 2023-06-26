

library(survival)
library(BART)
library(glmnet)
library(mboost)
library(devtools)

install_github("nchenderson/AFTrees")
library(AFTrees)

##  the RMST_BCART function
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

## Friedman test function

set.seed(123)

# sample function
f.test <- function(x) {10*sin(pi*x[ , 1]*x[ , 2]) + 20*(x[ , 3]-.5)^2+10*x[ , 4]+5*x[ , 5]}

sigma <- 1.0
n <- 250 # 250 or 2000 # number of training observation
n.test <- 2000 # 2000 or 4000 # number of test observation
num_covar <- 100 # 10 or 100 # total number of predictors
ndraws <- 500
sgrid <- seq(0, 10, by=.1)
## choosing this big tau value cause warning
#Warning message:
# In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :
#  collapsing to unique 'x' values
tau <- 5
gam_alph <- 20
nreps <- 3 # number of simulation replications

## function we need for cox model
#CoxExpectedSurv <- function(X, beta_val, H0fn, tau) {
## This function computes E( min(T_i, tau) |x_i) for
## a cox ph model
# mu.x <- colMeans(X)
#integrand <- function(time, xi, beta_val) {
# nu <- sum(xi*beta_val)
#ans <- exp(-H0fn(time)*exp(nu))
#return(ans)
#}
#nn <- nrow(X)
#fitted_vals <- rep(NA, nn)
#for(k in 1:nn) {
# II <- integrate(integrand, lower=0, upper=tau, xi=X[k,] - mu.x,
#                beta_val=beta_val, subdivisions=5000L)
#fitted_vals[k] <- II$value
#}
#return(fitted_vals)
#}

CoxExpectedSurv <- function(X, beta_val, time, H0.vals, tau) {
  ## This function computes E( min(T_i, tau) |x_i) for
  ## a cox ph model
  mu.x <- colMeans(X)
  tpoints <- c(time[time < tau], tau)
  Hpoints <- H0.vals[time < tau]
  nn <- nrow(X)

  fitted_vals <- rep(NA, nn)
  for(k in 1:nn) {
    nu <- sum((X[k,] - mu.x)*beta_val)
    fitted_vals[k] <- sum(exp(-Hpoints*exp(nu))*diff(tpoints))
  }
  return(fitted_vals)
}



cens_prop <- rep(NA, nreps)
rmse_bcart <- rmse_bart <- rmse_coxph <- rmse_rcoxph <- rmse_sboost <- rep(NA, nreps)
rmse_aft <- rmse_aft_bart<- rmse_aft_null <- rmse_ipcw <- rep(NA, nreps)

coverage_bcart <- coverage_bart <- coverage_aft_bart <- rep(NA, nreps)
for(j in 1:nreps) {
  ## training set
  X.train <- matrix(runif(n*num_covar), n, num_covar)
  colnames(X.train) <- paste0('X', 1:num_covar)
  ET.train <- f.test(X.train)
  #mu.train <- digamma(gam_alph) - log(ET.train)
  ## might need to input tau into this calculation?

  T.train <- rgamma(n, shape=gam_alph, rate=ET.train)
  #C.train <- runif(n, min=.5, max=3) ## min = 0.5 or 2
  ## Dependent censoring
  C.train <- runif(n, min=2, max=3) + 5*X.train[ , 1] + 5*X.train[ , 2] + 5*X.train[ , 3] + 5*X.train[ , 4]
  Y.train <- pmin(T.train, C.train)
  delta.train <- ifelse(T.train <= C.train, 1, 0) ## mean delta train 50-60 % or 80-90 %

  ## test set
  X.test <- matrix(runif(n.test*num_covar), n.test, num_covar)
  colnames(X.test) <- paste0('X', 1:num_covar)
  ET.test <- f.test(X.test)
  #mu.test <- digamma(gam_alph) - log(ET.test)
  mu.test <- (ET.test/gam_alph)*pgamma(tau, shape = gam_alph+1, rate = ET.test) +
    tau*pgamma(tau, shape = gam_alph, rate = ET.test, lower.tail = FALSE)
  T.test <- rgamma(n.test, shape=gam_alph, rate=ET.test)
  #C.test <- runif(n.test, min=.5, max=3) ## min = 0.5 or 2
  ## Dependent censoring
  C.test <- runif(n, min=2, max=3) + 5*X.test[ , 1] + 5*X.test[ , 2] + 5*X.test[ , 3] + 5*X.test[ , 4]
  Y.test <- pmin(T.test, C.test)
  delta.test <- ifelse(T.test <= C.test, 1, 0)

  ### 1. AFT linear model
  AFT <- survreg(Surv(Y.train, delta.train) ~ X.train)
  XX <- model.matrix(Y.test ~ X.test)
  aft_linpred <- as.numeric(XX%*%AFT$coefficients)
  aft_sigma <- AFT$scale
  aft_sigsq <- aft_sigma*aft_sigma
  ## exp(aft_linpred) is an approximate fitted value
  ## For RMST, a more precise definition of fitted value is

  gt_prob <- pnorm((log(tau) - aft_linpred)/aft_sigma, lower.tail=FALSE)
  lt_prob <- pnorm((log(tau) - aft_sigsq - aft_linpred)/aft_sigma)
  AFT_fitted <- exp(aft_sigsq/2 + aft_linpred)*lt_prob + tau*gt_prob
  ## do either bootstrap or derive confidence interval directly


  ##### 2. AFT intercept-only model
  AFT_null <- survreg(Surv(Y.train, delta.train) ~ 1)
  aft_linpred <- rep(AFT_null$coefficients[1], nrow(X.test))
  aft_sigma <- AFT_null$scale
  aft_sigsq <- aft_sigma*aft_sigma

  gt_prob <- pnorm((log(tau) - aft_linpred)/aft_sigma, lower.tail=FALSE)
  lt_prob <- pnorm((log(tau) - aft_sigsq - aft_linpred)/aft_sigma)
  AFT_null_fitted <- exp(aft_sigsq/2 + aft_linpred)*lt_prob + tau*gt_prob



  #### 3. AFT_BART model
  AFT_BART <- abart(X.train, Y.train, delta.train, x.test=X.test)
  ndraw_abart <- nrow(AFT_BART$yhat.test)
  AFT_fit_reps <- matrix(NA, nrow=ndraw_abart, ncol=nrow(X.test))
  for(k in 1:ndraw_abart) {
    aft_bart_mu <- AFT_BART$yhat.test[k,]
    aft_bart_sig <- AFT_BART$sigma[k]
    aft_bart_sigsq <- aft_bart_sig*aft_bart_sig
    ## exp(aft_linpred) is an approximate fitted value
    ## For RMST, a more precise definition of fitted value is

    gt_prob <- pnorm((log(tau) - aft_bart_mu)/aft_bart_sig, lower.tail=FALSE)
    lt_prob <- pnorm((log(tau) - aft_bart_sigsq - aft_bart_mu)/aft_bart_sig)

    AFT_fit_reps[k,] <- exp(aft_bart_sigsq/2 + aft_bart_mu)*lt_prob + tau*gt_prob
  }
  AFT_BART_fitted <- colMeans(AFT_fit_reps)
  AFT_BART_CI <- t(apply(AFT_fit_reps, 1, function(x) quantile(x, probs=c(0.025, 0.975))))


  ## 4. Boosting with IPCW weights
  ipcw_weights <- IPCweights(x=Surv(Y.train, delta.train))
  IPW_boost <- glmboost(x=X.train[delta.train==1,], y=pmin(Y.train[delta.train==1], tau),
                        weights = ipcw_weights[delta.train==1])
  IPW_fitted <- as.numeric(predict(IPW_boost, newdata=X.test))

  ## 5. Cox-PH model without penalization of regression coefficients
  COXPH.mod <- coxph(Surv(Y.train, delta.train) ~ X.train)
  #coxhaz <- basehaz(COXPH.mod)
  coxhaz <- survfit(COXPH.mod, x=X.train, y=Surv(Y.train, delta.train))
  H0fn <- approxfun(c(0, coxhaz$time), c(0, coxhaz$cumhaz),
                    yright=max(coxhaz$cumhaz))
  ## for some datasets returns this error: Error in integrate(integrand, lower = 0, upper = tau, xi = X[k, ] - mu.x, :
  #maximum number of subdivisions reached
  COXPH_fitted <- CoxExpectedSurv(X=X.test, beta_val=COXPH.mod$coefficients,
                                  time = coxhaz$time, H0.vals=coxhaz$cumhaz,
                                  tau=tau)

  ### 6. Cox-PH model with lasso penalty
  rcox_tmp <- cv.glmnet(x=X.train, y=Surv(Y.train, delta.train), family = "cox",
                        type.measure = "deviance")
  RCOXPH <-  glmnet(X.train, Surv(Y.train, delta.train), family = "cox",
                    lambda = rcox_tmp$lambda.min)
  Rcoxhaz <- survfit(RCOXPH, x=X.train, y=Surv(Y.train, delta.train))
  RH0fn <- approxfun(c(0, Rcoxhaz$time), c(0, Rcoxhaz$cumhaz),
                     yright=max(Rcoxhaz$cumhaz))
  lasso_betahat <- as.numeric(coef(RCOXPH))
  RCOXPH_fitted <- CoxExpectedSurv(X=X.test, beta_val=lasso_betahat,
                                   time=Rcoxhaz$time,
                                   H0.vals=Rcoxhaz$cumhaz, tau=tau)

  ## 7. RMST BCART
  #bcart_mod <- RMST_BCART(Y.train, delta.train, X.train, X.test,
   #                       ndraws=ndraws, tau=tau)
  ## If doing dependent censoring use:
  bcart_mod <- RMST_BCART(Y.train, delta.train, X.train, X.test,
                          ndraws=ndraws, ipcw="dependent", tau=tau)
  bcart_fitted <- rowMeans(bcart_mod$fitted.values.test)

  bcart_CI <- t(apply(bcart_mod$fitted.values.test, 1,
                    function(x) quantile(x, probs=c(0.025, 0.975))))

  ## 8. RMST BART
  #bart_mod <- RMST_BART(Y.train, delta.train, X.train, X.test,
  #                      ndraws=ndraws, tau=tau)
  ## If doing dependent censoring use:
  bart_mod <- RMST_BART(Y.train, delta.train, X.train, X.test,
                        ndraws=ndraws, ipcw="dependent", tau=tau)
  bart_fitted <- rowMeans(bart_mod$fitted.values.test)
  bart_CI <- t(apply(bart_mod$fitted.values.test, 1,
                     function(x) quantile(x, probs=c(0.025, 0.975))))

  ## Recording RMSE
  rmse_bcart[j] <- sqrt(mean((bcart_fitted - mu.test)*(bcart_fitted - mu.test)))
  rmse_bart[j] <- sqrt(mean((bart_fitted - mu.test)*(bart_fitted - mu.test)))
  rmse_coxph[j] <- sqrt(mean((COXPH_fitted - mu.test)*(COXPH_fitted - mu.test)))
  rmse_rcoxph[j] <- sqrt(mean((RCOXPH_fitted - mu.test)*(RCOXPH_fitted - mu.test)))
  #rmse_sboost[j] <- sqrt(mean((SBOOST_fitted - mu.test)*(SBOOST_fitted - mu.test)))
  rmse_ipcw[j] <- sqrt(mean((IPW_fitted - mu.test)*(IPW_fitted - mu.test)))
  rmse_aft[j] <- sqrt(mean((AFT_fitted - mu.test)*(AFT_fitted - mu.test)))
  rmse_aft_bart[j] <- sqrt(mean((AFT_BART_fitted - mu.test)*(AFT_BART_fitted - mu.test)))
  rmse_aft_null[j] <- sqrt(mean((AFT_null_fitted - mu.test)*(AFT_null_fitted - mu.test)))
  cens_prop[j] <- mean(delta.test) # also record censoring proportion

  ## Recording coverage
  coverage_aft_bart[j] <- mean((mu.test >= AFT_BART_CI[,1]) & (mu.test <= AFT_BART_CI[,2]))
  coverage_bcart[j] <- mean((mu.test >= bcart_CI[,1]) & (mu.test <= bcart_CI[,2]))
  coverage_bart[j] <- mean((mu.test >= bart_CI[,1]) & (mu.test <= bart_CI[,2]))
  ## report the results of coverage as a table.
}


nmethods <- 8
Results <- matrix(NA, nrow=nmethods, ncol=2)
rownames(Results) <- c("AFT Null", "CoxPH", "Cox glmnet", "AFT linear",
                       "ipwc", "AFT BART", "BCART", "BART")
colnames(Results) <- c("Mean RMSE", "Median RMSE")

Results[1,1:2] <- c(mean(rmse_aft_null), median(rmse_aft_null))
Results[2,1:2] <- c(mean(rmse_coxph), mean(rmse_coxph))
Results[3,1:2] <- c(mean(rmse_rcoxph), median(rmse_rcoxph))
Results[4,1:2] <- c(mean(rmse_aft), median(rmse_aft))
Results[5,1:2] <- c(mean(rmse_ipcw), median(rmse_ipcw))
Results[6,1:2] <- c(mean(rmse_aft_bart), median(rmse_aft_bart))
Results[7,1:2] <- c(mean(rmse_bcart), median(rmse_bcart))
Results[8,1:2] <- c(mean(rmse_bart), median(rmse_bart))

round(Results, 4)
