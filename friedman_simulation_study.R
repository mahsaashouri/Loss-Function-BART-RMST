

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

## Friedman test function

set.seed(123)

# sample function
f.test <- function(x) {10*sin(pi*x[ , 1]*x[ , 2]) + 20*(x[ , 3]-.5)^2+10*x[ , 4]+5*x[ , 5]}

sigma = 1.0
n = 250 # 250 or 2000 # number of training observation
n.test = 2000 # 2000 or 4000 # number of test observation
k = 10 # 10 or 100 # total number of predictors
ndraws <- 500
sgrid <- seq(0, 10, by=.1)
## choosing this big tau value cause warning
#Warning message:
 # In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :
#  collapsing to unique 'x' values
tau <- 500
gam_alph <- 20
nreps <- 5 # number of simulation replications

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


cens_prop <- rep(NA, nreps)
rmse_bcart <- rmse_coxph <- rmse_rcoxph <- rmse_sboost <- rmse_aft <- rmse_aft_bart<- rmse_aft_null <- rep(NA, nreps)
for(j in 1:nreps) {
    ## training set
    X.train <- matrix(runif(n*k), n, k)
    colnames(X.train) <- paste0('X', 1:k)
    ET.train <- f.test(X.train)
    mu.train <- digamma(gam_alph) - log(ET.train)
    ## might need to input tau into this calculation?

    T.train <- rgamma(n, shape=gam_alph, rate=ET.train)
    C.train <- runif(n, min=.5, max=3) ## min = 0.5 or 2
    Y.train <- pmin(T.train, C.train)
    delta.train <- ifelse(T.train <= C.train, 1, 0) ## mean delta train 50-60 % or 80-90 %

    ## test set
    X.test <- matrix(runif(n.test*k), n.test, k)
    colnames(X.test) <- paste0('X', 1:k)
    ET.test <- f.test(X.test)
    mu.test <- digamma(gam_alph) - log(ET.test)
    T.test <- rgamma(n.test, shape=gam_alph, rate=ET.test)
    C.test <- runif(n.test, min=.5, max=3) ## min = 0.5 or 2
    Y.test <- pmin(T.test, C.test)
    delta.test <- ifelse(T.test <= C.test, 1, 0)

    ## RMST-BCART
    bcart_mod <- RMST_BCART(Y.train, delta.train, X.train, X.test, ndraws=500, tau=500, sigma.mu=1.2)
    bcart_fitted <- pmin(rowMeans(bcart_mod$fitted.values.test), log(tau))

    ## Coxph
    COXPH.mod <- coxph(Surv(Y.train, delta.train) ~ X.train)
    coxhaz <- basehaz(COXPH.mod)
    H0fn <- approxfun(c(0, coxhaz$time), c(0, coxhaz$hazard),
                      yright=max(coxhaz$hazard))
    COXPH <- CoxExpectedSurv(X=X.test, beta_val=COXPH.mod$coefficients,
                    H0fn=H0fn, tau=1)
    COXPH_fitted <- pmin(COXPH, log(tau))

    ## regularized coxph model
    RCOXPH <-  glmnet(X.train, Surv(Y.train, delta.train), family = "cox", lambda = 1, alpha = 1)
    RCOXPH_fitted <- pmin(c(predict(RCOXPH, X.test, type = 'response')), log(tau))

    ## survival boosting
    SBOOST <- glmboost(Surv(Y.train, delta.train)~X.train, family = Gehan(), control = boost_control(mstop = 300))
    ## we have warnings in predict function: 'newdata' had 2000 rows but variables found have 250 rows
    SBOOST_fitted <- pmin(c(predict(SBOOST, newdata = data.frame(X.test))), log(tau))

    ## AFT
    AFT <- survreg(Surv(Y.train, delta.train) ~ X.train)
    ## we have warnings in predict function: 'newdata' had 2000 rows but variables found have 250 rows
    AFT_fitted <- pmin(predict(AFT, newdata = data.frame(X.test), type = 'response'), log(tau))
    #AFT_fitted <- pmin(AFT$linear.predictors, log(tau))

    ## AFT BART
    AFT_BART <- abart(X.train, Y.train, x.test = X.test, delta.train)
    AFT_BART_fitted <-  pmin(AFT_BART$yhat.test.mean, log(tau))

    ## AFT null
    AFT_null <- survreg(Surv(Y.train, delta.train) ~ 1)
    AFT_null_fitted <- pmin(predict(AFT_null, newdata = data.frame(X.test), type = 'response'), log(tau))

    rmse_bcart[j] <- sqrt(mean((bcart_fitted - mu.test)*(bcart_fitted - mu.test)))
    rmse_coxph[j] <- sqrt(mean((COXPH_fitted - mu.test)*(COXPH_fitted - mu.test)))
    rmse_rcoxph[j] <- sqrt(mean((RCOXPH_fitted - mu.test)*(RCOXPH_fitted - mu.test)))
    rmse_sboost[j] <- sqrt(mean((SBOOST_fitted - mu.test)*(SBOOST_fitted - mu.test)))
    rmse_aft[j] <- sqrt(mean((AFT_fitted - mu.test)*(AFT_fitted - mu.test)))
    rmse_aft_bart[j] <- sqrt(mean((AFT_BART_fitted - mu.test)*(AFT_BART_fitted - mu.test)))
    rmse_aft_null[j] <- sqrt(mean((AFT_null_fitted - mu.test)*(AFT_null_fitted - mu.test)))
    cens_prop[j] <- mean(delta.test) # also record censoring proportion
}


nmethods <- 7
Results <- matrix(NA, nrow=nmethods, ncol=2)
rownames(Results) <- c("AFT Null", "CoxPH", "Cox glmnet", "AFT linear",
                       "Boosting", "BCART", "AFT BART")
colnames(Results) <- c("Mean RMSE", "Median RMSE")

Results[1,1:2] <- c(mean(rmse_aft_null), median(rmse_aft_null))
Results[2,1:2] <- c(mean(rmse_coxph), mean(rmse_coxph))
Results[3,1:2] <- c(mean(rmse_rcoxph), median(rmse_rcoxph))
Results[4,1:2] <- c(mean(rmse_aft), median(rmse_aft))
Results[5,1:2] <- c(mean(rmse_sboost), median(rmse_sboost))
Results[6,1:2] <- c(mean(rmse_bcart), median(rmse_bcart))
Results[7,1:2] <- c(mean(rmse_aft_bart), median(rmse_aft_bart))

round(Results, 4)

