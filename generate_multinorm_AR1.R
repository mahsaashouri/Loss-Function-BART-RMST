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

## correlation matrix for data with AR(1) correlation structure
AR1Cor <- function(n, Rho) {
  expt <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  Rho^expt
}

## function to simulate linear model
sim.reg <- function(nobs, coef, mu, sd, Rho){ 
                    
  num.var = length(coef)  
  beta = as.matrix(coef)
  ## generate data from multivariate normal with AR(1) - using the above function and MASS package
  #H = mvrnorm(n = nobs, mu, Sigma = AR1Cor(num.var, Rho))
  ## generate data from multivariate normal with AR(1) - using 'arima.sim' function
  H = t(replicate(nobs, arima.sim(n = length(coef), model = list(ar = Rho))))
  Y = exp(H %*% beta + rnorm(nobs, 0, sd))
  return(list(Z = H, Y = c(Y), coeff = coef))
}

## example

n = 1000 # number of training observation   
k = 10  # total number of predictors  
coef = c(0.1,0,0.05,0.1,1,0.03,0.9,0,0.25,1)
mu = rep(0,k)
sd = 1
Rho = 0.9

# simulate training set
DataSim <- sim.reg(n = 1000, coef = coef, mu = mu, sd = sd, Rho = Rho) 
X.train <- DataSim$Z
colnames(X.train) <- paste0('X', 1:k)
alpha <- 1
C.train <- rgamma(n, shape =alpha, scale = alpha*DataSim$Y)
Y.train <- pmin(DataSim$Y, C.train)
delta.train <- ifelse(DataSim$Y <= C.train, 1, 0)
sgrid <- seq(0, 10, by=.1)

##run BCART 

train <- RMST_BCART(Y.train, delta.train, X.train, ntree=1, ndraws=1000, sigma.mu=1.2)


Y.min <- min(Y.train)
Y.max <- max(Y.train)
plot(rowMeans(train$fitted.values[,,1]), DataSim$Y, asp=1, pch='.',
     xlim=c(Y.min, Y.max), ylab='BCART')

