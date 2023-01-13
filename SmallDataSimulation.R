
## A way to generate a small dataset to test
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


n <- 500
xx <- runif(n)
mean_fn <- rep(NA, n)
mean_fn[xx < 1/3] <- 4
mean_fn[xx >= 1/3 & xx < 2/3] <- 6
mean_fn[xx > 2/3] <- 8

logT <- mean_fn + rnorm(n, mean=0, sd=2)
logCC <- rexp(n, rate=0.15)

Y <- pmin(logT, logCC)
Y <- pmin(Y, 10)
delta <- ifelse(logT <= logCC, 1, 0)
X <- matrix(xx, nrow=n, ncol=1)
colnames(X) <- 'x'
sgrid <- seq(0, 10, by=.1)
muvec <- 0

test_run <- RMST_BCART(Y, delta, X, ntree=1, ndraws=500, sigma.mu=1.2)


pmean <- rowMeans(test_run$fitted.values[,,1])
plot(xx[delta==1], pmean, ylim=c(0, 8))
lines(xx[order(xx)], mean_fn[order(xx)], type="s")


## BART

test <- list(dvec = Dmat[1,], splt.vars = c(), splt.vals = c())
old.tree <- list(test)[rep(1,10)]

test_run_BART <- RMST_BART(Y, delta, X, old.tree = old.tree, ndraws=500, sigma.mu=1.2, muvec=muvec)


