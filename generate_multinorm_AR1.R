

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
  H = mvrnorm(n = nobs, mu, Sigma = AR1Cor(num.var, Rho))
  Y = H %*% beta + rnorm(nobs, 0, sd)
  return(list(Z = H, Y = c(Y), coeff = coef))
}

## example

sim.reg(1000, c(0.1,0,0.05), rep(0,3), 0.1,0.8)
