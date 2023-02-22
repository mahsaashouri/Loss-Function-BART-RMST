

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

## Friedman test function

set.seed(123)

# sample function
f.test <- function(x) {10*sin(pi*x[ , 1]*x[ , 2]) + 20*(x[ , 3]-.5)^2+10*x[ , 4]+5*x[ , 5]}

sigma = 1.0
n = 1000 # number of training observation
k = 10  # total number of predictors
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
rmse_bcart <- rmse_coxph <- rmse_rcoxph <- rmse_sboost <- rmse_aft <- rmst_aft_bart<- rmse_aft_null <- rep(NA, nreps)
for(j in 1:nreps) {
    X.train <- matrix(runif(n*k), n, k)
    colnames(X.train) <- paste0('X', 1:k)
    ET.train <- f.test(X.train)
    mu.train <- digamma(gam_alph) - log(ET.train)
    ## might need to input tau into this calculation?

    T.train <- rgamma(n, shape=gam_alph, rate=ET.train)
    C.train <- runif(n, min=1, max=3)
    Y.train <- pmin(T.train, C.train)
    delta.train <- ifelse(T.train <= C.train, 1, 0)
    bcart_mod <- RMST_BCART(Y.train, delta.train, X.train, ndraws=500, tau=500, sigma.mu=1.2)
    bcart_fitted <- pmin(rowMeans(bcart_mod$fitted.values), log(tau))
    
    ## Coxph
    COXPH.mod <- coxph(Surv(Y.train, delta.train) ~ X.train)
    coxhaz <- basehaz(COXPH.mod)
    H0fn <- approxfun(c(0, coxhaz$time), c(0, coxhaz$hazard),
                      yright=max(coxhaz$hazard))
    COXPH <- CoxExpectedSurv(X=X.train, beta_val=COXPH.mod$coefficients,
                    H0fn=H0fn, tau=1)
    COXPH_fitted <- pmin(COXPH, log(tau))
    
    ## regularized coxph model
    RCOXPH <-  glmnet(X.train, Surv(Y.train, delta.train), family = "cox", lambda = 1, alpha = 1)
    RCOXPH_fitted <- pmin(c(predict(RCOXPH, X.train, type = 'response')), log(tau))
    
    ## survival boosting
    SBOOST <- glmboost(Surv(Y.train, delta.train)~X.train, family = Gehan(), control = boost_control(mstop = 100))
    SBOOST_fitted <- pmin(c(predict(SBOOST)), log(tau))

    ## AFT
    AFT <- survreg(Surv(Y.train, delta.train) ~ X.train)
    AFT_fitted <- pmin(AFT$linear.predictors, log(tau))
    
    ## AFT BART
    AFT_BART <- abart(X.train, Y.train, delta.train)
    AFT_BART_fitted <-  pmin(AFT_BART$yhat.train.mean, log(tau))

    ## AFT null
    AFT_null <- survreg(Surv(Y.train, delta.train) ~ 1)
    AFT_null_fitted <- pmin(AFT_null$linear.predictors, log(tau))

    rmse_bcart[j] <- sqrt(mean((bcart_fitted - mu.train)*(bcart_fitted - mu.train)))
    rmse_coxph[j] <- sqrt(mean((COXPH_fitted - mu.train)*(COXPH_fitted - mu.train)))
    rmse_rcoxph[j] <- sqrt(mean((RCOXPH_fitted - mu.train)*(RCOXPH_fitted - mu.train)))
    rmse_sboost[j] <- sqrt(mean((SBOOST_fitted - mu.train)*(SBOOST_fitted - mu.train)))
    rmse_aft[j] <- sqrt(mean((AFT_fitted - mu.train)*(AFT_fitted - mu.train)))
    rmst_aft_bart[j] <- sqrt(mean((AFT_BART_fitted - mu.train)*(AFT_BART_fitted - mu.train)))
    rmse_aft_null[j] <- sqrt(mean((AFT_null_fitted - mu.train)*(AFT_null_fitted - mu.train)))
    cens_prop[j] <- mean(delta.train) # also record censoring proportion
}




