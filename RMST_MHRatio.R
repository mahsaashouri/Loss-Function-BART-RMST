
RMST_MHRatio <- function(U, proposed_tree, old_tree, sigma.mu, Gvec, xmat, m, alpha, beta, ntree){
  ## transition ratio
  old.new <- Transition_Prob(old_tree, proposed_tree, xmat, m)
  m2 <- ifelse(m == 1, 2, ifelse(m == 2, 1, 3))
  new.old <- Transition_Prob(proposed_tree, old_tree, xmat, m2)
  ratio.transition <- new.old/old.new
  
  ## log-likelihood ratio
  loglike.old <- LogLik(old_tree, xmat, U, Gvec, sigma.mu)
  loglike.new <- LogLik(proposed_tree, xmat, U, Gvec, sigma.mu)
  LogLikRatio <- loglike.new/loglike.old
  
  ## prior ratio
  old.prior <- ProbD(U, xmat, old_tree$splt.vals, old_tree$splt.vars, muvec, alpha, beta, old_tree$d, ntree)
  new.prior <- ProbD(U, xmat, proposed_tree$splt.vals, proposed_tree$splt.vars, muvec, alpha, beta, proposed_tree$d, ntree)
  ratio.prior <- new.prior/old.prior
  
  ## MH ratio
  MHratio <- min(1, ratio.transition*LogLikRatio*ratio.prior)
  return(MHratio)
}