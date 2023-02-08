
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



RMST_BCART <- function(Y, delta, X, X.test=NULL, ndraws=500, sigma.mu=1,
                       alpha=0.95, beta=2, kappa0=1, sgrid=NULL, tau=NULL, burnIn=100) {
  ## skeleton of function for computing
  ## Bayesian CART for the RMST loss function

  ## organize data
  if(ncol(X) > 1) {
    xmat <- X[delta==1,]
  } else{
    xmat <- matrix(X[delta==1,], nrow=sum(delta==1), ncol=1)
    colnames(xmat) <- colnames(X)
  }

  if(is.null(tau)) {
    tau <- max(Y[delta==1])
  }
  U <- pmin(log(Y[delta==1]), tau)
  ## want to add an option for not taking the log
  #U <- pmin(Y, tau)
  ## get Gvec
  if(is.null(sgrid)) {
    sgrid <- c(0, exp(seq(tau/101, tau, length.out=100)))
  }
  SS <- ComputeSurvStatistics(sgrid=sgrid, times=exp(Y), status=1 - delta)
  lam.draw <- GPDraw(eU=exp(U), sgrid=sgrid, num.risk=SS$n.risk,
                     num.events=SS$n.event, kappa0=1)
  Gvec <- exp(-lam.draw)

  ## initialize tree
  n <- length(Y)

  FittedValues <- matrix(NA, nrow=n,  ncol=ndraws+burnIn)
  if(!is.null(X.test)) {
      FittedValues.test <- matrix(NA, nrow = nrow(X.test), ncol = ntree)
  } else {
      FittedValues.test <- NULL
  }
  NNodes <- loglikvals <- rep(NA, ndraws+burnIn+1)


  ## initialize trees
  old_tree <- list(dvec = Dmat[1,], splt.vars = c(), splt.vals = c())
  NNodes[1] <- sum(old_tree$dvec==1)
  loglikvals[1] <- LogLik(tree=old_tree, X=xmat, U=U,
                            Gvec=Gvec, sigma.mu=sigma.mu)
  muvec <- rep(0, NNodes[1])

  for(k in 1:(ndraws + burnIn)) {
     # step 1: Update tree
     ## sample one of three move types

     move_type <- sample(1:3, size=1)
     proposed_tree <- ProposedTree(move_type, old_tree, xmat)
     if ("character" %in% class(proposed_tree)){
       while(("character" %in% class(proposed_tree))){
          move_type <- sample(1:3, size=1)
          proposed_tree <- ProposedTree(move_type, old_tree, xmat)
       }
     }
     ## compute the ratio
     MH_ratio <- RMST_MHRatio(U = U, new_tree = proposed_tree, old_tree = old_tree, sigma.mu,
                              Gvec = Gvec, X = xmat, m = move_type, alpha, beta, ntree, tau = tau)
     u <- runif(1)
     if(u <= MH_ratio) {
        new_tree <- proposed_tree
     } else {
        new_tree <- old_tree
     }

     ## Step 2: update the mu values -  sample from a normal distribution
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
     FittedValues[,k] <- FittedValue(X, new_tree$splt.vals, new_tree$splt.vars,
                                     muvec, new_tree$dvec)
     NNodes[k+1] <- sum(new_tree$dvec==1)
     loglikvals[k+1] <- LogLik(tree=new_tree, X=xmat, U=U, Gvec=Gvec, sigma.mu=sigma.mu)
     old_tree <- new_tree
     if(!is.null(X.test)) {
        FittedValues.test[,k] <- FittedValue(X.test, new_tree$splt.vals, new_tree$splt.vars,
                                             muvec, new_tree$dvec)
     }
  }

  Fitted.Values <- FittedValues[,(burnIn+1):(ndraws+burnIn)]

  ans <- list(fitted.values=Fitted.Values, nnodes=NNodes,
              logliks=loglikvals, fitted.values.test=FittedValues.test)
  return(ans)
}

