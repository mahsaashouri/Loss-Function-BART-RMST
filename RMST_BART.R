
## function to return the proposed tree

ProposedTree <- function(move_type, old_tree, xmat){
  if(move_type==1) {
    proposed_tree <- GrowMove(old_tree, xmat)
  } else if(move_type==2) {
    proposed_tree <- PruneMove(old_tree, xmat)
  } else if(move_type==3) {
    proposed_tree <- ChangeMove(old_tree, xmat)
  }
  return(proposed_tree)
}

## old.tree is a list of initial trees 

RMST_BART <- function(Y, delta, X, old.tree, ndraws, sigma.mu, muvec,sgrid, alpha, beta, num.risk, num.events, kappa0) {
 
  ## organize data
  if(ncol(X) > 1) {
    xmat <- X[delta==1,]
  } else{
    xmat <- matrix(X[delta==1,], nrow=sum(delta==1), ncol=1)
    colnames(xmat) <- colnames(X)
  }
  U <- Y[delta==1]
  
  ## get Gvec
  SS <- ComputeSurvStatistics(sgrid=sgrid, times=Y, status=1 - delta)
  
  lam.draw <- GPDraw(U=U, sgrid=sgrid, num.risk=SS$n.risk,
                     num.events=SS$n.event, kappa0=1)
  Gvec <- exp(-lam.draw)
  

  n <- length(U)
  tau <- (max(U) - min(U))/(2*sqrt(length(old.tree)))
   
  ## initialize fitted values
  FittedValues <- matrix(NA, nrow = n,  ncol = length(old.tree))
  for(i in 1:length(old.tree)){
    FittedValues[,i] <- FittedValue(xmat, old.tree[[i]]$splt.vals, old.tree[[i]]$splt.vars, muvec, old.tree[[i]]$dvec)
  }
  
  NNodes <- loglikvals <- matrix(NA, nrow = ndraws, ncol = length(old.tree))
  
  Fitted.Values <- list()
  
  for(j in 1:ndraws){
      for(k in 1:length(old.tree)) {
        
      # update tree
      ## sample one of three move types
      
      move_type <- sample(1:3, size=1)
      proposed_tree <- ProposedTree(move_type,  old.tree[[k]], xmat)
      if ("character" %in% class(proposed_tree)){
        while(("character" %in% class(proposed_tree))){
          move_type <- sample(1:3, size=1)
          proposed_tree <- ProposedTree(move_type,  old.tree[[k]], xmat)
        }
      }
      
      # update Rs (U.res)
      U.res <- U - (rowSums(FittedValues) - FittedValues[,k]) 
      
      ## compute the ratio
      MH_ratio <- RMST_MHRatio(U = U.res, new_tree = proposed_tree, old_tree = old.tree[[k]], muvec = muvec, sigma.mu,
                               Gvec, X = xmat, m = move_type, alpha, beta, length(old.tree), tau = tau)
      u <- runif(1)
      if(u <= MH_ratio) {
        new_tree <- proposed_tree
      } else {
        new_tree <- old.tree[[k]]
      }
      
      ## update the mu values -  sample from a normal distribution
      ## number of terminal nodes in the new tree
      terminal_nodes <- which(new_tree$dvec==2)
      
      ## get mean and sigma for updating mu values
      AT <- AMatrix(xmat, new_tree$splt.vals, new_tree$splt.vars, new_tree$dvec)
      WTGDiag <- c(crossprod(AT, 1/Gvec))
      VG <- U/Gvec
      Z <- c(crossprod(AT, VG))
      
      ## update mu values
      
      mu.mean <- Z/(WTGDiag + 1/(sigma.mu*sigma.mu))
      mu.sd <- sqrt(1/(WTGDiag + 1/(sigma.mu*sigma.mu)))
      muvec <- c(mu.mean) + mu.sd*rnorm(length(terminal_nodes))
      
      ## Record fitted values at each step
      FittedValues[,k] <- FittedValue(xmat, new_tree$splt.vals, new_tree$splt.vars, muvec, new_tree$dvec)
      
      ## record number of internal nodes and loglikelihood values
      NNodes[j,k] <- sum(new_tree$dvec==1)
      loglikvals[j,k] <- LogLik(tree=new_tree, X=xmat, U=U.res, Gvec=Gvec, sigma.mu=sigma.mu)
      
      ## replace old_tree with new_tree
      old_tree[[k]] <- new_tree
      }
    Fitted.Values[[length(Fitted.Values)+1]] <- FittedValues
  }
  ans <- list(fitted.values=Fitted.Values, nnodes=NNodes, logliks=loglikvals)
  return(ans)
}
