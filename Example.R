
source('RMST_BCART.R')
source('A_matrix.R')
source('GrowMove.R')
source('ChangeMove.R')
source('PruneMove.R')
source('D_probabilities.R')
source('fitted_value.R')
source('GPDraw.R')
source('LogLik.R')
source('prior_conditional_on_D.R')
source('RMST_MHRatio.R')
source('Transition_Prob.R')


library(survival)

## German breast cancer data
data(cancer, package="survival")
head(gbsg)

## define X matrix
X <- gbsg[ , !(names(gbsg) %in% c('status'))]
xmat <- model.matrix(rfstime~., data = X)
Y <- gbsg$rfstime
delta <- 1 - gbsg$status
sgrid <- seq(0, 500, by=50)
alpha <- .95
beta <- 2
ntree <- 20
kappa0 <- 1

SS <- ComputeSurvStatistics(sgrid, times=gbsg$rfstime,
                      status=as.numeric(gbsg$status==1))

num.risk <- SS$n.risk
num.events<- SS$n.event

RMST_BCART(Y, delta, X, tree, ndraws, sigma.mu, sgrid, alpha, beta, ntree,
           num.risk, num.events, kappa0)
  

