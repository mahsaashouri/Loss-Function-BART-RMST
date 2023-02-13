
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
ndraws <- 500
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
train <- RMST_BCART(Y.train, delta.train, X.train, ndraws=500, sigma.mu=1.2)

## run BART
test <- list(dvec = Dmat[1,], splt.vars = c(), splt.vals = c())
old.tree <- list(test)[rep(1,5)]

train.BART <- RMST_BART(Y.train, delta.train, X.train, old.tree, ndraws=50, sigma.mu=1.2)

## arrange BART fitted values
fitted.values.m <- matrix(NA, nrow = n, ncol = ndraws)
fitted.values.s <- matrix(NA, nrow = length(old.tree), ncol = n)
for(i in 1:length(old.tree)){
  for(j in 1:ndraws){
    fitted.values.m[,j] <- train.BART$fitted.values[[j]][,i]
  }
  fitted.values.s[i,] <- rowMeans(fitted.values.m)
}

## plot BCART and BART
# BCART
plot(ET.train, rowMeans(train$fitted.values))

plot(rowMeans(train$fitted.values[,,1]), ET.train, asp=1, pch='.',
     xlim=c(min(Y.train), max(Y.train)), ylab='BCART')
# BART
plot(colMeans(fitted.values.s), ET.train)

plot(colMeans(fitted.values.s), ET.train, asp=1, pch='.',
     xlim=c(min(Y.train), max(Y.train)), ylab='BART')




################
## other methods
################

###############
## coxph model
###############
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

###############
## regularized coxph model with glmnet
###############
library(glmnet)
rcoxph <-  glmnet(X.train, Surv(Y.train, delta.train), family = "cox", lambda = 0)

###############
## basic AFT model
###############
library(survival)
BAFT <- survreg(Surv(Y.train, delta.train) ~ X.train)

tau <- 20
pre1 <- pmin(BAFT$linear.predictors, log(tau))
pre2 <- pmin(exp(BAFT$linear.predictors), tau)
plot(pmin(ET.train, tau), pre2)


###############
## survival boosting
###############
## Using mboost package
library(mboost)
SurveMBoost <- glmboost(Surv(Y.train, delta.train)~X.train, family = Gehan(), control = boost_control(mstop = 100))
# Make predictions
predSurveMBoost <- predict(SurveMBoost)#, type = "response")
# Plot the model's performance
plot(survfit(Surv(Y.train, delta.train) ~ predSurveMBoost), lty = 1:2, mark.time = FALSE)

###############
## AFT BART
###############
library(BART)
AFTBART <- abart(X.train, Y.train, delta.train)

