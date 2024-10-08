
## Inputs: responses, predictors, splitting values and variables, mu vector, alpha, beta
## tree structure (row number of Dmat matrix 1:26), number of trees -- we  also need the
## function which computes the tree probabilities

## output: P(C (splitting values), V (splitting variables), mu|D (tree structure))

ProbD <- function(U, X, splt.vals, splt.vars, #muvec
                  alpha, beta, d, tau){
  n <- length(U)
  l <- length(splt.vars)
  PrD <- (DProb(alpha, beta, d))#^(l)
  PrV <- ncol(X)^(-l)
  #PrMu <- prod((1/tau)*dnorm(muvec/tau))
  PrM <- 1; Pr <- c();
  for(i in 1:l){
    xm <- X[,splt.vars[i]]
    Pr[i] <- mean( which(xm == splt.vals[i]) )
    PrM <- Pr[i]*PrM
  }
  #return(PrM*PrMu*PrV*PrD)
  return(PrM*PrV*PrD)
}

## Example

#U <- c( 17, 14, 2, 0.4, 11)
#xmat <- matrix(round(rnorm(20,10),0), ncol = 4)
#colnames(xmat) <- c('x1', 'x2', 'x3', 'x4')
#splt.vals <- c('x1' = xmat[1,1], 'x3' = xmat[2,3], 'x2' = xmat[5,2], 'x4' = xmat[4,4])
#muvec <- c(0.1, 0.03, 0.2, 0, 1, 4)
#splt.vars <- c('x1', 'x3', 'x2','x4')
#dvec <- c(1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2)

#ProbD(U = U, X = xmat, splt.vals = splt.vals, splt.vars = splt.vars, muvec = muvec,
#      alpha = .95, beta = 2, d = 16, tau = 1)
