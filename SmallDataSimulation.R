
## A way to generate a small dataset to test
##  the RMST_BCART function
source("A_matrix.R")
source("ChangeMove.R")
source("D_probabilities.R")
source("GPDraw.R")
source("GrowMove.R")
source("LogLik.R")
source("PruneMove.R")
source("PruneMove_V2.R")
source("RMST_BCART.R")
source("RMST_Loss.R")
source("RMST_MHRatio.R")
source("Recursive_A_matrix.R")
source("Transition_Prob.R")
source("fitted_value.R")
source("prior_conditional_on_D.R")
source("prior_conditional_on_D_V2.R")
source("tree-configuration.R")


n <- 100
xx <- runif(n)
mean_fn <- rep(NA, n)
mean_fn[xx < 1/3] <- 4
mean_fn[xx >= 1/3 & xx < 2/3] <- 6
mean_fn[xx > 2/3] <- 8

logT <- mean_fn + rnorm(n, mean=0, sd=2)
logCC <- rexp(n, rate=0.2)

Y <- pmin(logT, logCC)
Y <- pmin(Y, 10)
delta <- ifelse(logT <= logCC, 1, 0)
X <- matrix(xx, nrow=n, ncol=1)
colnames(X) <- 'x'
sgrid <- seq(0, 10, by=.1)


test_run <- RMST_BCART(Y, delta, X, ntree=1, ndraws=5, sigma.mu=1.2, muvec=muvec,
                       sgrid=sgrid, alpha=0.95, beta=2, num.risk=0, num.events=0, kappa0=1)




