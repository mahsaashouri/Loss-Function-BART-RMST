
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

# simulate training set
X.train <- matrix(runif(n*k), n, k)
colnames(X.train) <- paste0('X', 1:k)
ET.train <- f.test(X.train)
T.train <- ET.train+sigma*rexp(n)
C.train <- rexp(n, rate = .05)
Y.train <- pmin(T.train, C.train)
delta.train <- ifelse(T.train <= C.train, 1, 0)

# simulate test set
#m <- 100 # number of test observation
#X.test <- matrix(runif(m*k), m, k)
#colnames(X.test) <- paste0('X', 1:k)
#ET.test <- f.test(X.test)
#T.test <- ET.test+sigma*rexp(m)
#C.test <- rexp(m)
#Y.test <- pmin(T.test, C.test)
#delta.test <- ifelse(Y.train <= C.train, 1, 0)

sgrid <- seq(0, 10, by=.1)

##run BCART 

train <- RMST_BCART(Y.train, delta.train, X.train, ntree=1, ndraws=1000, sigma.mu=1.2)
#test <- predict(train, X.test)


plot(rowMeans(train$fitted.values), ET.train)


Y.min <- min(Y.train)
Y.max <- max(Y.train)
plot(rowMeans(train$fitted.values[,,1]), ET.train, asp=1, pch='.',
          xlim=c(Y.min, Y.max), ylab='BCART')
#plot(rowMeans(train$fitted.values[,,1]), ET.train[delta.train==1], asp=1, pch='.',
#     xlim=c(Y.min, Y.max), ylab='BCART')


## coxph model
library(survival)
coxph <- coxph(Surv(Y.train, delta.train) ~ X.train)
plot(survfit(coxph))
coxhaz <- basehaz(coxph)
plot(coxhaz$time, coxhaz$hazard)

## regularized coxph model with glmnet
library(glmnet)
rcoxph <-  glmnet(X.train, Surv(Y.train, delta.train), family = "cox")
#plot(survfit(rcoxph, x = X.train, y = Surv(Y.train, delta.train)))

## basic AFT model
library(survival)
BAFT <- survreg(Surv(Y.train, delta.train) ~ X.train)

tau <- 20

pre1 <- pmin(BAFT$linear.predictors, log(tau))
pre2 <- pmin(exp(BAFT$linear.predictors), tau)
plot(pmin(ET.train, tau), pre2)

## penalized AFT model

## survival boosting - Hothorn paper

## AFT BART
library(BART)
AFTBART <- abart(X.train, Y.train, delta.train)
