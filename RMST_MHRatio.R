
RMST_MHRatio <- function(U, new_tree, old_tree, muvec,sigma.mu, Gvec, xmat, m, alpha, beta, ntree){
  ## transition ratio
  old.new <- Transition_Prob(old_tree, new_tree, xmat, m)
  m2 <- ifelse(m == 1, 2, ifelse(m == 2, 1, 3))
  new.old <- Transition_Prob(new_tree, old_tree, xmat, m2)
  ratio.transition <- new.old/old.new
  
  ## log-likelihood ratio
  loglike.old <- LogLik(old_tree, xmat, U, Gvec, sigma.mu)
  loglike.new <- LogLik(new_tree, xmat, U, Gvec, sigma.mu)
  LogLikRatio <- loglike.new/loglike.old
  
  ## prior ratio
  old.prior <- ProbD(U, xmat, old_tree$splt.vals, old_tree$splt.vars, muvec, alpha, beta, old_tree$d, ntree)
  new.prior <- ProbD(U, xmat, new_tree$splt.vals, new_tree$splt.vars, muvec, alpha, beta, new_tree$d, ntree)
  ratio.prior <- new.prior/old.prior
  
  ## MH ratio
  MHratio <- min(1, ratio.transition*LogLikRatio*ratio.prior)
  return(MHratio)
}