
RMST_BCART <- function(Y, delta, X, ntree, ndraws, sigma.mu, muvec,sgrid, alpha, beta, num.risk, num.events, kappa0) {
  ## Skeleton of function for computing
  ## Bayesian CART for the RMST loss function
  
  ## organize data
  if(ncol(X) > 1) {
    xmat <- X[delta==1,]
  } else{
    xmat <- matrix(X[delta==1,], nrow=sum(delta==1), ncol=1)
    colnames(xmat) <- colnames(X)
  }
  U <- Y[delta==1]
  
  ## Get Gvec
  SS <- ComputeSurvStatistics(sgrid=sgrid, times=Y, status=1 - delta)
  
  lam.draw <- GPDraw(U=U, sgrid=sgrid, num.risk=SS$n.risk,
                     num.events=SS$n.event, kappa0=1)
  Gvec <- exp(-lam.draw)
  
  ## initialize tree
  n <- length(U)
  
  FittedValues <- array(NA, c(Nrow=n,  Ncol=ndraws, Ntree=ntree))
  NNodes <- loglikvals <- matrix(NA, nrow=ndraws+1, ncol = ntree)
  
  tau <- (max(U) - min(U))/(2*sqrt(ntree))

  for(j in 1:ntree){
    ## initialize trees
    old_tree <- list(dvec = Dmat[1,], splt.vars = c(), splt.vals = c())
    NNodes[1,j] <- sum(old_tree$dvec==1)
    loglikvals[1,j] <- LogLik(tree=old_tree, X=xmat, U=U, Gvec=Gvec, sigma.mu=sigma.mu)
    ## old_tree holds the splitting vars, splitting values, and dvec
    for(k in 1:ndraws) {
      # Step 1: Update Tree
      ## Sample one of three move types
      move_type <- sample(1:3, size=1)
      proposed_tree <- ProposedTree(move_type, old_tree, xmat)
      if ("character" %in% class(proposed_tree)){
        while(("character" %in% class(proposed_tree))){
          move_type <- sample(1:3, size=1)
          proposed_tree <- ProposedTree(move_type, old_tree, xmat)
        }
      }
      MH_ratio <- RMST_MHRatio(U = U, new_tree = proposed_tree, old_tree = old_tree, muvec = muvec, sigma.mu,
                               Gvec, X = xmat, m = move_type, alpha, beta, ntree, tau = tau)
      u <- runif(1)
      if(u <= MH_ratio) {
        new_tree <- proposed_tree
      } else {
        new_tree <- old_tree
      }
      
      ## Step 2: Update the mu values -  sample from a Normal distribution
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
      FittedValues[,k,j] <- FittedValue(xmat, new_tree$splt.vals, new_tree$splt.vars, muvec, new_tree$dvec)
      NNodes[k+1,j] <- sum(new_tree$dvec==1)
      loglikvals[k+1,j] <- LogLik(tree=new_tree, X=xmat, U=U, Gvec=Gvec, sigma.mu=sigma.mu)
      old_tree <- new_tree
    }
  }
  ans <- list(fitted.values=FittedValues, nnodes=NNodes, logliks=loglikvals)
  return(ans)
}