
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
tau <- 50000
gam_alph <- 20
nreps <- 5 # number of simulation replications

cens_prop <- rep(NA, nreps)
rmse_bcart <- rmse_aft <- rmse_aft_null <- rep(NA, nreps)
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
    bcart_mod <- RMST_BCART(Y.train, delta.train, X.train, ndraws=500, tau=tau, sigma.mu=1.2)
    bcart_fitted <- pmin(rowMeans(bcart_mod$fitted.values), log(tau))

    AFT <- survreg(Surv(Y.train, delta.train) ~ X.train)
    AFT_fitted <- pmin(AFT$linear.predictors, log(tau))

    AFT_null <- survreg(Surv(Y.train, delta.train) ~ 1)
    AFT_null_fitted <- pmin(AFT_null$linear.predictors, log(tau))

    rmse_bcart[j] <- sqrt(mean((bcart_fitted - mu.train)*(bcart_fitted - mu.train)))
    rmse_aft[j] <- sqrt(mean((AFT_fitted - mu.train)*(AFT_fitted - mu.train)))
    rmse_aft_null[j] <- sqrt(mean((AFT_null_fitted - mu.train)*(AFT_null_fitted - mu.train)))
    cens_prop[j] <- mean(delta.train) # also record censoring proportion
}




