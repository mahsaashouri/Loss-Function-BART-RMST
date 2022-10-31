
## Inputs: responses, predictors, splitting values and variables, mu vector, alpha, beta
## tree structure (row number of Dmat matrix 1:26), number of trees -- we  also need the 
## function which computes the tree probabilities

## output: P(C (splitting values), V (spliting variables), mu|D (tree structure))

ProbD <- function(logY, xmat, splt.vals, splt.vars, muvec, alpha, beta, d, ntree){
  n <- length(logY)
  l <- length(splt.vars)
  tau <- (max(logY) - min(logY))/(2*sqrt(20))
  PrD <- (DProb(alpha, beta, d))^(l)
  PrV <- (1/ncol(xmat))^(l)
  PrMu <- prod((1/tau)*dnorm(muvec/tau))
  PrM <- 1; Pr <- c();
  for(i in 1:l){
    x <- xmat[,splt.vars[i]]
    Pr[i] <- sum(which(x == x[i]))/n
    PrM <- Pr[i]*PrM
  }
  return(PrM*PrMu*PrV*PrD)
}

## Example

logY <- c( 17, 14, 2, 0.4, 11)
xmat <- matrix(round(rnorm(20,10),0), ncol = 4)
colnames(xmat) <- c('x1', 'x2', 'x3', 'x4')
splt.vals <- c('x1' = .05, 'x3' = .5, 'x2' = 1, 'x4' = .02)
muvec <- c(0.1, 0.03, 0.2, 0, 1, 4)
splt.vars <- c('x1', 'x3', 'x2','x4')
dvec <- c(1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2)

ProbD(logY = logY, xmat = xmat, splt.vals = splt.vals, splt.vars = splt.vars, muvec = muvec, 
      alpha = .95, beta = 2, d = 16, ntree = 20)
