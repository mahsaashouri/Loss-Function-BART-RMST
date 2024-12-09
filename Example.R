source('RMST_BCART.R')
source('RMST_BART.R')
source('A_matrix.R')
source('GrowMove.R')
source('ChangeMove.R')
source('PruneMove.R')
source('D_probabilities.R')
source('fitted_value.R')
source('GPDraw.R')
source('LogLik.R')
source('prior_conditional_on_D_V2.R')
source('RMST_MHRatio.R')
source('Transition_Prob.R')
source('tree-configuration.R')

library(survival)

## German breast cancer data
data(cancer, package="survival")
head(gbsg)

## define X matrix
X <- gbsg[ , !(names(gbsg) %in% c('status'))]
X <- model.matrix(rfstime~.-1, data = X)
Y <- gbsg$rfstime
delta <- gbsg$status
sgrid <- seq(0, 4000, by=1)
alpha <- .95
beta <- 2
ntree <- 20
kappa0 <- 1
ndraws <- 500
sigma.mu <- 1.2


splt.vals.raw <- c('pid' = 148, 'meno' = 1, 'nodes' = 12, 'pid' = 700, 'size'=15)
splt.vars.raw <- c('pid', 'meno', 'nodes','pid','size')
muvec <- c(0.1, 0.03, 0.2, 0, 1, 4)
dvec <- c(1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2)
tree <-  list(dvec = dvec, splt.vars = splt.vars.raw, splt.vals = splt.vals.raw)
## BCART

RMST_BCART(Y, delta, X, ntree = 1, ndraws = 10, sigma.mu)

## BART

test <- list(dvec = Dmat[1,], splt.vars = c(), splt.vals = c())
old.tree <- list(test)[rep(1,10)]

RMST_BART(Y, delta, X, old.tree, ndraws = 10, sigma.mu)



