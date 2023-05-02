
library(survival)
library(BART)
library(glmnet)
library(mboost)
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


set.seed(123)


sigma <- 1.0
n <- 250 # 250 or 2000 # number of training observation
n.test <- 2000 # 2000 or 4000 # number of test observation
num_covar <- 10 # 10 or 100 # total number of predictors
ndraws <- 500
sgrid <- seq(0, 10, by=.1)
coef <- runif(num_covar)
Rho <- 0.5
## choosing this big tau value cause warning
#Warning message:
# In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :
#  collapsing to unique 'x' values
tau <- 50
gam_alpha <- 20
nreps <- 1 # number of simulation replications

## function to simulate linear model AR1
sim.reg <- function(nobs, coef, mu, sd, Rho){
  num.var = length(coef)
  beta = as.matrix(coef)
  ## generate data from multivariate normal with AR(1) - using the above function and MASS package
  #H = mvrnorm(n = nobs, mu, Sigma = AR1Cor(num.var, Rho))
  ## generate data from multivariate normal with AR(1) - using 'arima.sim' function
  H = t(replicate(nobs, c(arima.sim(length(coef), model = list(ar = Rho)))))
  Y = exp(H %*% beta)
  return(list(Z = H, Y = c(Y), coeff = coef))
}

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
    fitted_vals[k] <- sum( exp(-Hpoints*exp(nu)*diff(tpoints)) )
  }
  return(fitted_vals)
}



cens_prop <- rep(NA, nreps)
rmse_bcart <- rmse_bart <- rmse_coxph <- rmse_rcoxph <- rmse_sboost <- rep(NA, nreps)
rmse_aft <- rmse_aft_bart<- rmse_aft_null <- rmse_ipcw <- rep(NA, nreps)
for(j in 1:nreps) {
  ## training set
  DataSim <- sim.reg(n, coef = coef, mu = mu, Rho = Rho)
  X.train <- DataSim$Z
  colnames(X.train) <- paste0('X', 1:num_covar)
  T.train <- rgamma(n, shape = gam_alpha, rate = DataSim$Y)
  C.train <- runif(n, min=10, max=50) ## max = 50 or 2000
  Y.train <- pmin(T.train, C.train)
  delta.train <- ifelse(T.train <= C.train, 1, 0) ## mean delta train 50-60 % or 80-90 %
  #mu.train <- digamma(gam_alpha) - log(DataSim$Y)

  ## test set
  DataSim.test <- sim.reg(n.test, coef = coef, mu = mu, Rho = Rho)
  X.test <- DataSim.test$Z
  colnames(X.test) <- paste0('X', 1:num_covar)
  T.test <- rgamma(n.test, shape = gam_alpha, rate = DataSim.test$Y)
  C.test <- runif(n, min=10, max=50) ## max = 50 or 2000
  Y.test <- pmin(T.test, C.test)
  delta.test <- ifelse(T.test <= C.test, 1, 0) ## mean delta train 50-60 % or 80-90 %
  #mu.test <- digamma(gam_alpha) - log(DataSim.test$Y)
  mu.test <- (DataSim.test$Y/gam_alpha)*pgamma(tau, shape = gam_alpha+1, rate = DataSim.test$Y) +
    tau*pgamma(tau, shape = gam_alpha, rate = DataSim.test$Y, lower.tail = FALSE)

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

  ## 4. Boosting with IPCW weights
  ipcw_weights <- IPCweights(x=Surv(Y.train, delta.train))
  IPW_boost <- glmboost(x=X.train[delta.train==1,], y=pmin(Y.train[delta.train==1], tau),
                        weights = ipcw_weights[delta.train==1])
  IPW_fitted <- as.numeric(predict(IPW_boost, newdata=X.test))

  ## 5. Cox-PH model without penalization of regression coefficients
  COXPH.mod <- coxph(Surv(Y.train, delta.train) ~ X.train)
  #coxhaz <- basehaz(COXPH.mod)
  coxhaz <- survfit(COXPH.mod, x=X.train, y=Surv(Y.train, delta.train))

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

  lasso_betahat <- as.numeric(coef(RCOXPH))
  RCOXPH_fitted <- CoxExpectedSurv(X=X.test, beta_val=lasso_betahat,
                                   time=Rcoxhaz$time,
                                   H0.vals=Rcoxhaz$cumhaz, tau=tau)

  ## 7. RMST BCART
  bcart_mod <- RMST_BCART(Y.train, delta.train, X.train, X.test,
                          ndraws=ndraws, tau=tau)
  bcart_fitted <- rowMeans(bcart_mod$fitted.values.test)

  ## 8. RMST BART
  bart_mod <- RMST_BART(Y.train, delta.train, X.train, X.test,
                        ndraws=ndraws, tau=tau)
  bart_fitted <- rowMeans(bart_mod$fitted.values.test)

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
