
RMST_MHRatio <- function(U, new_tree, old_tree, sigma.mu, Gvec,
                         X, m, alpha, beta){
  ## transition ratio
  old.new <- Transition_Prob(old_tree, new_tree, X, m)
  m2 <- ifelse(m == 1, 2, ifelse(m == 2, 1, 3))
  new.old <- Transition_Prob(new_tree, old_tree, X, m2)
  ratio.transition <- new.old/old.new

  ## log-likelihood ratio
  loglike.old <- LogLik(old_tree, X, U, Gvec, sigma.mu)
  loglike.new <- LogLik(new_tree, X, U, Gvec, sigma.mu)
  LikRatio <- exp(loglike.new-loglike.old)

  ## prior ratio

  ## find the dvec row number (in Dmat matrix - all the possible tree structures matrix) for old and new tree
  Dmat <- matrix(NA, nrow = 26, ncol = 15)
  for (k in 1:26) {
    Dmat[k,] <- FindDvec(k)
  }
  d1 <- which(apply(Dmat, 1, function(x) all.equal(x[1:ncol(Dmat)], old_tree$dvec)) == "TRUE")
  d2 <- which(apply(Dmat, 1, function(x) all.equal(x[1:ncol(Dmat)], new_tree$dvec)) == "TRUE")

  old.prior <- ProbD(U, X, old_tree$splt.vals, old_tree$splt.vars, alpha, beta, d1)
  new.prior <- ProbD(U, X, new_tree$splt.vals, new_tree$splt.vars, alpha, beta, d2)
  ratio.prior <- new.prior/old.prior

  ## MH ratio
  MHratio <- min(1, ratio.transition*LikRatio*ratio.prior)
  return(MHratio)
}
