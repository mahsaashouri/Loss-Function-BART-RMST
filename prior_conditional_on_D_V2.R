ProbD <- function(U, X, splt.vals, splt.vars, muvec,
                  alpha, beta, d, tau){
  n <- length(U)
  l <- length(splt.vars)
  PrD <- DProb(alpha, beta, d)
  log_PrMu <- sum( dnorm(muvec, mean=0, sd=tau, log=TRUE) )
  if(l==0) {
     ans <- exp(log_PrMu + log(PrD))
  } else{
     log_PrV <- -l*log( ncol(X) )
     log_PrM <- 0
     for(i in 1:l){
       xm <- X[,splt.vars[i]]
       splt_prob <- mean( which(xm == splt.vals[i]) )
       log_PrM <- log(splt_prob) + log_PrM
     }
     ans <- exp(log_PrM + log_PrMu + log_PrV + log(PrD))
  }
  return(ans)
}

#U <- c( 17, 14, 2, 0.4, 11)
#xmat <- matrix(round(rnorm(20,10),0), ncol = 4)
#colnames(xmat) <- c('x1', 'x2', 'x3', 'x4')
#splt.vals <- c('x1' = xmat[1,1], 'x3' = xmat[2,3], 'x2' = xmat[5,2], 'x4' = xmat[4,4])
#muvec <- c(0.1, 0.03, 0.2, 0, 1, 4)
#splt.vars <- c('x1', 'x3', 'x2','x4')
#dvec <- c(1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2)

#tst <- ProbD(U = U, X = xmat, splt.vals = splt.vals, splt.vars = splt.vars, muvec = muvec,
#      alpha = .95, beta = 2, d = 16, tau = 1)



