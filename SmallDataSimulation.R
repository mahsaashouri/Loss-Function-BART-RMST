
## A way to generate a small dataset to test
##  the RMST_BCART function
n <- 100
xx <- runif(n)
TT <- rep(NA, n)
TT[xx < 1/3] <- 500 + rnorm(sum(xx < 1/3), mean=0, sd=100)
TT[xx >= 1/3 & xx < 2/3] <- 750 + rnorm(sum(xx >= 1/3 & xx < 2/3), mean=0, sd=100)
TT[xx >= 2/3] <- 1000 + rnorm(sum(xx >= 2/3), mean=0, sd=100)
CC <- rexp(n, rate=0.002)

Y <- pmin(TT, CC)
Y <- pmin(Y, 1100)
delta <- ifelse(TT <= CC, 1, 0)
X <- matrix(xx, nrow=n, ncol=1)

