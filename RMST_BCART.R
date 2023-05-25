
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



RMST_BCART <- function(U, delta, X, X.test=NULL, ndraws=100, transformation="identity",
                       ipcw="independent", sigma.mu=NULL, alpha=0.95, beta=2, kappa0=1,
                       sgrid=NULL, tau=NULL, burnIn=100) {
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
    tau <- max(U[delta==1])
  }
  U_tau <- pmin(U[delta==1], tau)

  if(is.null(sgrid)) {
    sgrid <- seq(0, tau, length.out=100)
  }
  ## Draw Gvec weights here.
  delta_alpha <- 1

  Gmat <- matrix(0, nrow=ndraws + burnIn, ncol=length(U_tau))
  if(ipcw=="independent") {
      for(k in 1:(ndraws + burnIn)) {
          Gmat[k,] <- DrawIPCW(U=U, delta=delta, Utau=U_tau, sgrid=sgrid,
                               kappa0=kappa0, delta_alpha=delta_alpha)
       }
  } else {
       cens_bart <- AFTrees(x.train=X, y.train=U, status=1-delta,
                            ndpost=ndraws + burnIn, verbose=FALSE)
       nd <- ndraws + burnIn
       Mucens_draws <- cens_bart$m.train[,delta==1]
       for(k in 1:nd) {
         for(j in 1:length(U_tau)) {
           log.time.points <- log(U_tau[j])
           AA <- (log.time.points - cens_bart$locations[k,] - Mucens_draws[k,j])/cens_bart$sigma[k]
           Cprob <- sum(pnorm(AA, lower.tail=FALSE)*cens_bart$mix.prop[k,])
           Gmat[k,j] <- 1/Cprob
         }
       }
  }

  ## Get KM estimate of censoring distribution and KM inverse censoring weights
  KM_cens <- survfit(Surv(U, 1 - delta) ~ 1)
  GKMfn <- stepfun(c(0, KM_cens$time), c(1, KM_cens$surv, min(KM_cens$surv)))
  GKM_weights <- GKMfn(U_tau)
  if(sum(GKM_weights == 0) > 0) {
    GKM_weights[GKM_weights==0] <- min(KM_cens$surv[KM_cens$surv > 0])
  }


  ## initialize tree
  n <- length(U)

  FittedValues <- matrix(NA, nrow=n,  ncol=ndraws+burnIn)
  if(!is.null(X.test)) {
      FittedValues.test <- matrix(NA, nrow = nrow(X.test), ncol = ndraws+burnIn)
  } else {
      FittedValues.test <- NULL
  }
  NNodes <- loglikvals <- rep(NA, ndraws+burnIn+1)
  if(transformation=="identity") {
    ## compute muhatb
    muhatb <- mean(U_tau/GKM_weights)
    Y_tau <- U_tau - muhatb
    if(is.null(sigma.mu)) {
      Ymin <- min(U_tau)
      sigma.mu <- (tau - muhatb - Ymin)/4
    }
  } else if(transformation=="log") {
    ## compute muhatb
    muhatb <- mean(log(U_tau)/GKM_weights)
    Y_tau <- log(U_tau) - muhatb
    if(is.null(sigma.mu)) {
      Ymin <- min(U_tau)
      sigma.mu <- (log(tau) - muhatb - log(Ymin))/4
    }
  }

  ## initialize trees
  old_tree <- list(dvec = FindDvec(1), splt.vars = c(), splt.vals = c())
  NNodes[1] <- sum(old_tree$dvec==1)
  loglikvals[1] <- LogLik(tree=old_tree, X=xmat, U=Y_tau,
                            Gvec=Gmat[1,], sigma.mu=sigma.mu)
  muvec <- rep(0, NNodes[1])

  for(k in 1:(ndraws + burnIn)) {
     # step 1: Update tree
     ## sample one of three move types
     Gvec <- Gmat[k,]

     move_type <- sample(1:3, size=1)
     proposed_tree <- ProposedTree(move_type, old_tree, xmat)
     if ("character" %in% class(proposed_tree)){
       while(("character" %in% class(proposed_tree))){
          move_type <- sample(1:3, size=1)
          proposed_tree <- ProposedTree(move_type, old_tree, xmat)
       }
     }
     ## compute the ratio
     MH_ratio <- RMST_MHRatio(U = Y_tau, new_tree = proposed_tree, old_tree = old_tree,
                              sigma.mu=sigma.mu, Gvec=Gvec, X = xmat, m = move_type,
                              alpha=alpha, beta=beta)
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
     VG <- Y_tau/Gvec
     Z <- c(crossprod(AT, VG))

     ## update mu values
     mu.mean <- Z/(WTGDiag + 1/(sigma.mu*sigma.mu))
     mu.sd <- sqrt(1/(WTGDiag + 1/(sigma.mu*sigma.mu)))
     muvec <- c(mu.mean) + mu.sd*rnorm(length(terminal_nodes))

     ## Record fitted values at each step
     FittedValues[,k] <- FittedValue(X, new_tree$splt.vals, new_tree$splt.vars,
                                     muvec, new_tree$dvec)
     NNodes[k+1] <- sum(new_tree$dvec==1)
      loglikvals[k+1] <- LogLik(tree=new_tree, X=xmat, U=Y_tau, Gvec=Gvec, sigma.mu=sigma.mu)
     old_tree <- new_tree
     if(!is.null(X.test)) {
        FittedValues.test[,k] <- FittedValue(X.test, new_tree$splt.vals, new_tree$splt.vars,
                                             muvec, new_tree$dvec)
     }
  }


  Fitted.Values <- FittedValues[,(burnIn+1):(ndraws+burnIn)] + muhatb
  FittedValues.test <- FittedValues.test[,(burnIn+1):(ndraws+burnIn)] + muhatb

  ans <- list(fitted.values=Fitted.Values, nnodes=NNodes,
              logliks=loglikvals, fitted.values.test=FittedValues.test)
  return(ans)
}

