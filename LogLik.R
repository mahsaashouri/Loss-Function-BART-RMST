
## Output: log density U|T,G - used to calculate log likelihood ratio in MH-Ratio (alpha)

source('A_matrix.R')


LogLik <- function(tree, xmat, U, Gvec, sigma.mu){

  AT <- AMatrix(xmat, tree$splt.vals, tree$splt.vars, tree$dvec)
  WTGDiag <- c(crossprod(AT, 1/Gvec))
  VG <- U/Gvec

  FE <- -(1/2)*sum((U*U)/Gvec)

  Z <- c(crossprod(AT, VG))

  SE <- (1/2)*sum(((sigma.mu^2)*(Z*Z))/(1+(sigma.mu^2)*WTGDiag))
  return(FE + SE)
}

xvec1 <- c('x1'=0, 'x2'=0, 'x3'=0, 'x4'=0, 'x5'=0) ## should be assigned to node 4
xvec2 <- c('x1'=0, 'x2'=0, 'x3'=1, 'x4'=0, 'x5'=0) ## should be assigned to node 5
xvec3 <- c('x1'=1, 'x2'=0, 'x3'=1, 'x4'=0, 'x5'=0) ## should be assigned to node 12
xvec4 <- c('x1'=1, 'x2'=0, 'x3'=1, 'x4'=1, 'x5'=0) ## should be assigned to node 13
xvec5 <- c('x1'=0.6, 'x2'=0, 'x3'=1, 'x4'=1, 'x5'=5) ## should be assigned to node 14
xvec6 <- c('x1' = 1, 'x2' = 2, 'x3' = 3, 'x4' = 10, 'x5' = 5) ## should be assigned to node 15
splt.vals.raw <- c('x1' = .05, 'x3' = .5)
splt.vars.raw <- c('x1', 'x3')
dvec <- c(1, 1, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
##
old_tree <- list(dvec = dvec, splt.vars = splt.vars.raw, splt.vals = splt.vals.raw)
X <- rbind.data.frame(xvec1, xvec5, xvec3, xvec6, xvec4, xvec2)
colnames(X) <- c('x1', 'x2', 'x3', 'x4', 'x5')

Gvec <- 1/rexp(6)
U <- rnorm(6)

LogLik(tree=old_tree, xmat=X, U=U, Gvec=Gvec, sigma.mu=1.2)


