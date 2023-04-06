
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

RMST_BART <- function(U, delta, X, X.test=NULL, ndraws=100, transformation="identity",
                      sigma.mu=NULL, ntrees=50, alpha=0.95, beta=2, kappa0=1,
                      sgrid=NULL, tau=NULL, burnIn=100) {
  ## xmat is the matrix that only contains the rows of
  ## X where delta=1
  if(ncol(X) > 1) {
    xmat <- X[delta==1,]
  } else{
    xmat <- matrix(X[delta==1,], nrow=sum(delta==1), ncol=1)
    colnames(xmat) <- colnames(X)
  }
  ## As a default, set tau to maximum of observed survival times
  if(is.null(tau)) {
    tau <- max(U[delta==1])
  }
  U_tau <- pmin(U[delta==1], tau)

  if(is.null(sgrid)) {
    sgrid <- seq(0, tau, length.out=100)
  }
  ## Draw Gvec weights here.
  delta_alpha <- 1
  Gvec <- DrawIPCW(U=U, delta=delta, Utau=U_tau, sgrid=sgrid,
                    kappa0=kappa0, delta_alpha=delta_alpha)

  ## Get KM estimate of censoring distribution and KM inverse censoring weights
  KM_cens <- survfit(Surv(U, 1 - delta) ~ 1)
  GKMfn <- stepfun(c(0, KM_cens$time), c(1, KM_cens$surv, min(KM_cens$surv)))
  GKM_weights <- GKMfn(U_tau)


  ## Setup storage for fitted values
  n <- length(U)
  FittedValues <- matrix(0, nrow = n,  ncol = ntrees)
  if(!is.null(X.test)) {
    FittedValues.test <- matrix(0, nrow = nrow(X.test),  ncol = ndraws)
    FitValue.test.track <- matrix(0, nrow = nrow(X.test), ncol = ntrees)
  } else {
    FittedValues.test <- NULL
  }
  NNodes <- loglikvals <- matrix(NA, nrow = (ndraws+burnIn), ncol = ntrees)

  ## Define centered version of U_tau
  ## and set sigma.mu if its values was not input to the function
  if(transformation=="identity") {
    ## compute muhatb
    muhatb <- mean(U_tau/GKM_weights)
    Y_tau <- U_tau - muhatb
    if(is.null(sigma.mu)) {
      Ymin <- min(U_tau)
      sigma.mu <- (tau - muhatb - Ymin)/(4*sqrt(ntrees))
    }
  } else if(transformation=="log") {
    ## compute muhatb
    muhatb <- mean(log(U_tau)/GKM_weights)
    Y_tau <- log(U_tau) - muhatb
    if(is.null(sigma.mu)) {
      Ymin <- min(U_tau)
      sigma.mu <- (log(tau) - muhatb - log(Ymin))/(4*sqrt(ntrees))
    }
  }

  ## Initialize old.tree list here.
  ##   old.tree is a list of trees with ntrees components
  Dprobs <- rep(NA, 26)
  for(i in 1:26){
    Dprobs[i] <- DProb(alpha = .95, beta = 2, i)
  }
  Dprobs <- Dprobs/sum(Dprobs)
  old.tree <- list()
  for(h in 1:ntrees) {
    old.tree[[h]] <- list()
    tree_structure <- sample(1:26, size=1, prob=Dprobs)
    ## need dvec, splt.vals, and splt.vars components of old.tree
    old.tree[[h]]$dvec <- FindDvec(tree_structure)
    ## Get number of internal nodes
    old.nnodes <- sum(old.tree[[h]]$dvec==1)
    ## sample splt.vars
    if(old.nnodes == 0){
      old.tree[[h]]$splt.vars <- numeric(0)
      old.tree[[h]]$splt.vals <- numeric(0)
    }
    else{
    splt.vars <- sample(colnames(X), size=old.nnodes, replace=TRUE)
    old.tree[[h]]$splt.vars <- splt.vars
    ## sample splt.vals
    splt.vals <- c()
    for(m in 1:length(splt.vars)){
      candidate_splitval <- unique(X[,splt.vars[m]])
      ww <- table(X[,splt.vars[m]])/nrow(X)
      splt.vals[m] <- candidate_splitval[sample(length(candidate_splitval), size=1,prob=ww)]
    }
    old.tree[[h]]$splt.vals <- splt.vals
    }
  }


  Fitted.Values <- matrix(NA, nrow=n, ncol=ndraws)
  for(j in 1:(ndraws+burnIn)){
    for(k in 1:ntrees) {

      # update tree k
      ## sample one of three move types
      move_type <- sample(1:3, size=1)
      proposed_tree <- ProposedTree(move_type,  old.tree[[k]], xmat)
      if("character" %in% class(proposed_tree)){
        while(("character" %in% class(proposed_tree))){
          move_type <- sample(1:3, size=1)
          proposed_tree <- ProposedTree(move_type,  old.tree[[k]], xmat)
        }
      }
      # update tree-k residuals
      Y.res <- Y_tau - (rowSums(FittedValues[delta==1,]) - FittedValues[delta==1,k])

      ## compute the MH ratio
      MH_ratio <- RMST_MHRatio(U = Y.res, new_tree = proposed_tree, old_tree = old.tree[[k]],
                               sigma.mu=sigma.mu, Gvec=Gvec, X = xmat, m = move_type,
                               alpha=alpha, beta=beta)
      u <- runif(1)
      if(u <= MH_ratio) {
        new_tree <- proposed_tree
      } else {
        new_tree <- old.tree[[k]]
      }
      ## update the mu values -  sample from a normal distribution
      ## number of mu values = number of terminal nodes in the new tree
      terminal_nodes <- which(new_tree$dvec==2)

      ## get mean and sigma for updating mu values
      AT <- AMatrix(xmat, new_tree$splt.vals, new_tree$splt.vars, new_tree$dvec)
      WTGDiag <- c(crossprod(AT, 1/Gvec))
      VG <- Y.res/Gvec
      Z <- c(crossprod(AT, VG))

      ## update mu values
      mu.mean <- Z/(WTGDiag + 1/(sigma.mu*sigma.mu))
      mu.sd <- sqrt(1/(WTGDiag + 1/(sigma.mu*sigma.mu)))
      muvec <- c(mu.mean) + mu.sd*rnorm(length(terminal_nodes))

      ## Update "fitted values" for tree k
      FittedValues[,k] <- FittedValue(X, new_tree$splt.vals, new_tree$splt.vars,
                                      muvec, new_tree$dvec)

      ## record number of internal nodes and loglikelihood values
      NNodes[j,k] <- sum(new_tree$dvec==1)
      #loglikvals[j,k] <- LogLik(tree=new_tree, X=xmat, U=U.res, Gvec=Gvec, sigma.mu=sigma.mu)
      if(!is.null(X.test)) {
        FitValue.test.track[,k] <- FittedValue(X.test, new_tree$splt.vals,
                                               new_tree$splt.vars, muvec, new_tree$dvec)
      }
      ## replace old_tree with new_tree
      old.tree[[k]] <- new_tree
    }
    if(j > burnIn) {
      Fitted.Values[,j-burnIn] <- rowSums(FittedValues)
      if(!is.null(X.test)) {
        FittedValues.test[,j-burnIn] <- rowSums(FitValue.test.track)
      }
    }
  }
  ## Add back in muhatb at the end.
  ans <- list(fitted.values=Fitted.Values[(burnIn+1):(ndraws+burnIn)] + muhatb,
              nnodes=tail(NNodes, ndraws),
              logliks=tail(loglikvals, ndraws),
              fitted.values.test=FittedValues.test)
  return(ans)
}
