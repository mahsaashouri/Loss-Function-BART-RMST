
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
coxph_mod <- coxph(Surv(Y.train, delta.train) ~ X.train)
plot(survfit(coxph_mod))
coxhaz <- basehaz(coxph_mod)
plot(c(0, coxhaz$time), c(0, coxhaz$hazard))
H0fn <- approxfun(c(0, coxhaz$time), c(0, coxhaz$hazard),
                  yright=max(coxhaz$hazard))

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

coxmod_mus <- CoxExpectedSurv(X=X.train, beta_val=coxph_mod$coefficients,
                              H0fn=H0fn, tau=25)

plot(pmin(ET.train, 25), coxmod_mus)

## regularized coxph model with glmnet
library(glmnet)
rcoxph <-  glmnet(X.train, Surv(Y.train, delta.train), family = "cox")


## basic AFT model
library(survival)
BAFT <- survreg(Surv(Y.train, delta.train) ~ X.train)

tau <- 20
pre1 <- pmin(BAFT$linear.predictors, log(tau))
pre2 <- pmin(exp(BAFT$linear.predictors), tau)
plot(pmin(ET.train, tau), pre2)

## penalized AFT model

## survival boosting - Hothorn paper
library(mboost)
SurveBoost <- glmboost(Surv(Y.train, delta.train)~X.train, family = CoxPH(),
                control=boost_control(mstop = 500))
# plot Survival Curves for a Cox Proportional Hazards Model
plot(survFit(SurveBoost))


## AFT BART
library(BART)
AFTBART <- abart(X.train, Y.train, delta.train)
