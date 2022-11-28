
RMST_BCART <- function(U, delta, xmat, tree, ndraws, sigma.mu, sgrid, alpha, beta, ntree, num.risk, num.events, kappa0) {
  ## Skeleton of function for computing
  ## Bayesian CART for the RMST loss function
  
  ## initialize tree
  n <- length(U)
  FittedValues <- matrix(NA, nrow=n, ncol=ndraws)
  old_tree <- tree
  ## old_tree holds the splitting vars, splitting values, and dvec
  for(k in 1:ndraws) {
    # Step 1: Update Tree
    ## Sample one of three move types
    move_type <- sample(1:3, size=1)
    if(move_type==1) {
      proposed_tree <- GrowMove(old_tree, xmat)
    } else if(move_type==2) {
      proposed_tree <- PruneMove(old_tree, xmat)
    } else if(move_type==3) {
      proposed_tree <- ChangeMove(old_tree, xmat)
    }
    MH_ratio <- RMST_MHRatio(U, proposed_tree, old_tree, sigma.mu, Gvec, xmat, move_type, alpha, beta, ntree)
    u <- runif(1)
    if(u <= MH_ratio) {
      new_tree <- proposed_tree
    } else {
      new_tree <- old_tree
    }
    ## Step 2: Update the mu values -  sample from a multivariate Normal distribution
    ## number of terminal nodes in the new tree
    terminal_nodes <- which(new_tree$dvec==2)
    ## get mean and sigma for updating mu values
    DG <- diag(Gvec)
    AT <- AMatrix(xmat, new_tree$splt.vals, new_tree$splt.vars, new_tree$dvec)
    WTG <- t(AT) %*% solve(DG) %*% AT
    VG <- solve(DG) %*% U
    Z <- c()
    for(k in 1:nrow(AT)){
      for(j in 1:ncol(AT)){
        tZ <- AT[k,j]*VG[k,]
        sumZ <- tZ +sumZ
      }
      Z[j] <- sumZ
    }
    mu.vec <- rmvnorm(length(terminal_nodes), (WTG+(sigma.mu^2*diag(1, length(terminal_nodes))))%*%matrix(Z, ncol=1), 
                                               solve(WTG+(sigma.mu^(-2)*diag(1,length(terminal_nodes))))) 
      
    ## Step 3: Update the vector of Gs
    Gvec <- GPDraw(U, sgrid, num.risk, num.events, kappa0) 
    
    ## Record fitted values at each step
    FittedValues[,k] <- FittedValue(xmat, new_tree$splt.vals, new_tree$splt.vars, mu.vec, new_tree$dvec)
  }
  return(FittedValues)
}


