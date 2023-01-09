
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

## To do:
##    1. I don't think muvec needs to be input into the function.
##    2. Need to give a burn-in number and discard burn-in iterations in returned result

RMST_BCART <- function(Y, delta, X, ntree, ndraws, sigma.mu, muvec,
                       alpha=0.95, beta=2, kappa0=1, sgrid=NULL, tau=NULL, burnIn = 100) {
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
  U <- pmin(Y[delta==1], tau)

  ## get Gvec
  if(is.null(sgrid)) {
      sgrid <- c(0, exp(seq(tau/101, tau, length.out=100)))
  }
  SS <- ComputeSurvStatistics(sgrid=sgrid, times=exp(Y), status=1 - delta)
  lam.draw <- GPDraw(eU=exp(U), sgrid=sgrid, num.risk=SS$n.risk,
                     num.events=SS$n.event, kappa0=1)
  Gvec <- exp(-lam.draw)

  ## initialize tree
  n <- length(U)

  FittedValues <- array(NA, c(Nrow=n,  Ncol=(ndraws+burnIn), Ntree=ntree))
  NNodes <- loglikvals <- matrix(NA, nrow=(ndraws+burnIn+1), ncol = ntree)

  for(j in 1:ntree){
    ## initialize trees
    old_tree <- list(dvec = Dmat[1,], splt.vars = c(), splt.vals = c())
    NNodes[1,j] <- sum(old_tree$dvec==1)
    loglikvals[1,j] <- LogLik(tree=old_tree, X=xmat, U=U, Gvec=Gvec, sigma.mu=sigma.mu)
    muvec <- rep(0, NNodes[1,j])

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
      MH_ratio <- RMST_MHRatio(U = U, new_tree = proposed_tree, old_tree = old_tree, muvec = muvec, sigma.mu,
                               Gvec, X = xmat, m = move_type, alpha, beta, ntree, tau = tau)
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
      FittedValues[,k,j] <- FittedValue(xmat, new_tree$splt.vals, new_tree$splt.vars, muvec, new_tree$dvec)
      NNodes[k+1,j] <- sum(new_tree$dvec==1)
      loglikvals[k+1,j] <- LogLik(tree=new_tree, X=xmat, U=U, Gvec=Gvec, sigma.mu=sigma.mu)
      old_tree <- new_tree
    }
  }
  Fitted.Values <- array(NA, c(Nrow=n,  Ncol=ndraws, Ntree=ntree))
  for(i in 1:ntree){
    Fitted.Values[,,i] <- FittedValues[,c((burnIn+1):(ndraws+burnIn)),i]
  }
  ans <- list(fitted.values=Fitted.Values, nnodes=tail(NNodes, ndraws), 
              logliks=tail(loglikvals, ndraws))
  return(ans)
}

