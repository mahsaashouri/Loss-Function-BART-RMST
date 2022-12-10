
RMST_BCART <- function(Y, delta, X, tree, ndraws, sigma.mu, muvec,sgrid, alpha, beta, ntree, num.risk, num.events, kappa0) {
  ## Skeleton of function for computing
  ## Bayesian CART for the RMST loss function
   
  ## organize data
  xmat <- X[delta==1,]
  U <- Y[delta==1]

  ## Get Gvec
  SS <- ComputeSurvStatistics(sgrid=sgrid, times=Y, status=1 - delta)

  lam.draw <- GPDraw(U=U, sgrid=sgrid, num.risk=SS$n.risk,
                     num.events=SS$n.event, kappa0=1)
  Gvec <- exp(-lam.draw)

  ## initialize tree
  n <- length(U)
  FittedValues <- matrix(NA, nrow=n, ncol=ndraws)
  old_tree <- tree
  
  tau <- (max(U) - min(U))/(2*sqrt(ntree))
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
    MH_ratio <- RMST_MHRatio(U = U, new_tree = proposed_tree, old_tree = old_tree, muvec = muvec, sigma.mu,
                             Gvec, X = xmat, m = move_type, alpha, beta, ntree, tau = tau)
    # if(is.na(MH_ratio) == TRUE){
    #   MH_ratio <- 1
    # }
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
    #WTG <- t(AT)%*%solve(diag(Gvec))%*%AT
    WTGDiag <- c(crossprod(AT, 1/Gvec))
    VG <- U/Gvec
    Z <- c(crossprod(AT, VG))

    ## update mu values
    #mu.mean <- (WTG+(sigma.mu^2*diag(1, length(terminal_nodes))))%*%matrix(Z, ncol=1)
    #mu.var <- solve(WTG+(sigma.mu^(-2)*diag(1,length(terminal_nodes))))
    #muvec <- c(mu.mean) + diag(sqrt(mu.var))*rnorm(length(terminal_nodes))
    mu.mean <- Z/(WTGDiag + 1/(sigma.mu*sigma.mu))
    mu.sd <- sqrt(1/(WTGDiag + 1/(sigma.mu*sigma.mu)))
    muvec <- c(mu.mean) + mu.sd*rnorm(length(terminal_nodes))


    ## Record fitted values at each step
    FittedValues[,k] <- FittedValue(xmat, new_tree$splt.vals, new_tree$splt.vars, muvec, new_tree$dvec)
  }
  return(FittedValues)
}


