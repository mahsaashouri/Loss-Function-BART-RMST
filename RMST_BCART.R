RMST_BCART <- function(U, delta, xmat, ndraws) {
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
    
    ## Step 2: Update the mu values
    mu.vec <- ## sample from a multivariate Normal distribution
      
      ## Step 3: Update the vector of Gs
      Gvec <- SampleG() 
    
    ## Record fitted values at each step
    FittedValues[,k] <- fitted_value_function(new_tree, mu.vec)
  }
  return(FittedValues)
}
