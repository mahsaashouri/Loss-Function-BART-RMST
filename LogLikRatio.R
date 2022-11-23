
## Output: log density U|T,G - used to calculate log likelihood ratio in MH-Ratio (alpha)

source('A_matrix_function.R')

LogLik <- function(tree, xmat, U, Gvec, sigma.mu){
  
  DG <- diag(Gvec)
  AT <- AMatrix(xmat, tree$splt.vals, tree$splt.vars, tree$dvec)
  WTG <- t(AT) %*% solve(DG) %*% AT
  WTGDiag <- WTG[row(WTG) == col(WTG)]
  VG <- solve(DG) %*% U
  
  i = 1:length(U)
  FE <- -(1/2)*sum((U[i]^2)/Gvec[i])
  
  Z <- c()
  for(k in 1:nrow(AT)){
    for(j in 1:ncol(AT)){
      tZ <- AT[k,j]*VG[k,]
      sumZ <- tZ +sumZ
    }
    Z[j] <- sumZ
  }
  
  j = 1:ncol(AT)
  SE <- (1/2)*sum(((sigma.mu^2)*(Z[j]^2))/(1+(sigma.mu^2)*WTGDiag[j]))
 return(FE + SE)
}

LogLikRatio <- LogLik(new_tree, xmat, U, Gvec, sigma.mu)/LogLik(old_tree, xmat, U, Gvec, sigma.mu)  
