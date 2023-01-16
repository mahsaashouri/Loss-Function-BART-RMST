
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
n = 10000 # number of training observation   
k = 500  # total number of predictors  

# simulate training set
X.train <- matrix(runif(n*k), n, k)
ET.train <- f.test(X.train)
T.train <- ET.train+sigma*rexp(n)
C.train <- rexp(n)
Y.train <- pmin(T.train, C.train)
delta.train <- ifelse(Y.train <= C.train, 1, 0)

# simulate test set
m <- 1000 # number of test observation
X.test <- matrix(runif(m*k), m, k)
ET.test <- f.test(X.test)
T.test <- ET.test+sigma*rexp(m)
C.test <- rexp(m)
Y.test <- pmin(T.test, C.test)
delta.test <- ifelse(Y.train <= C.train, 1, 0)

sgrid <- seq(0, 10, by=.1)

##run BCART 

train <- RMST_BCART(Y.train, delta, X.train, ntree=1, ndraws=1000, sigma.mu=1.2)





