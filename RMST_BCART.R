RMST_BCART <- function(U, delta, xmat, ndraws, mu.p, sigma.mu, sgrid, alpha, beta, ntree, num.risk, num.events, kappa0) {
  ## Skeleton of function for computing
  ## Bayesian CART for the RMST loss function
  
  ## initialize tree
  n <- length(U)
  FittedValues <- matrix(NA, nrow=n, ncol=ndraws)
  old_tree <- Dmat[1,]
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
    
    ## number of terminal nodes in the new tree
    terminal_nodes <- which(new_tree$dvec==2)
    
    ## Step 2: Update the mu values -  sample from a multivariate Normal distribution
    mu.vec <- rmvnorm(length(terminal_nodes), mu.p, sigma.mu) 
      
    ## Step 3: Update the vector of Gs
    Gvec <- GPDraw(U, sgrid, num.risk, num.events, kappa0) 
    
    ## Record fitted values at each step
    FittedValues[,k] <- FittedValue(xmat, new_tree$splt.vals, new_tree$splt.vars, mu.vec, new_tree$dvec)
  }
  return(FittedValues)
}
