

## Inputs: matrix of covariates, Splitting values and variables, mu vectors , and tree structure (from 26 possible structures)
## Output: Fitted mu vector

LeftChild <- function(x){
  y <- c(2, 4, 6, 8, 10, 12, 14)
  return(y[x])
}
RightChild <- function(x){
  y <- c(3, 5, 7, 9, 11, 13, 15)
  return(y[x])
}


FittedValue <- function(X, splt.vals.raw, splt.vars.raw, muvec.raw, dvec){
  if(dvec[1] == 2){
    non_terminal <- FALSE
    mu.hat <- muvec.raw
  }
  else{
    ## update mu vector
    muvec <- replace(dvec, dvec!=2, NA)
    muvec[ which(!is.na(muvec))] <- muvec.raw

    ## update splitting variable vector
    splt.vars <- replace(dvec, dvec!=1, NA)
    splt.vars[ which(!is.na(splt.vars))] <- splt.vars.raw
    ## update splitting value vector

    splt.vals <- replace(dvec, dvec!=1, NA)
    splt.vals[ which(!is.na(splt.vals))] <- splt.vals.raw
  mu.hat <- c()
  for(i in 1:nrow(X)){
    xvec <- X[i,]
    non_terminal <- TRUE
    current.node <- 1
    current.var <- splt.vars[1]
    while(non_terminal){
      if(xvec[current.var] <= splt.vals[current.node]){
        new_node <- LeftChild(current.node)
        new.var <- splt.vars[new_node]
      } else {
        new_node <- RightChild(current.node)
        new.var <- splt.vars[new_node]
      }
      if(dvec[new_node] == 2){
        non_terminal <- FALSE
        mu.hat_ <- muvec[new_node]
      } else {
        current.var <- new.var
        current.node <- new_node
      }
    }
    mu.hat <- c(mu.hat, mu.hat_)
  }
  }
  return(mu.hat)
}

## Example
#xvec1 <- c('x1'=0, 'x2'=0, 'x3'=0, 'x4'=0, 'x5'=0) ## should be assigned to node 4
#xvec2 <- c('x1'=0, 'x2'=0, 'x3'=1, 'x4'=0, 'x5'=0) ## should be assigned to node 5
#xvec3 <- c('x1'=1, 'x2'=0, 'x3'=1, 'x4'=0, 'x5'=0) ## should be assigned to node 12
#xvec4 <- c('x1'=1, 'x2'=0, 'x3'=1, 'x4'=1, 'x5'=0) ## should be assigned to node 13
#xvec5 <- c('x1'=0.6, 'x2'=0, 'x3'=1, 'x4'=1, 'x5'=5) ## should be assigned to node 14
#xvec6 <- c('x1' = 1, 'x2' = 2, 'x3' = 3, 'x4' = 10, 'x5' = 5) ## should be assigned to node 15
#splt.vals.raw <- c('x1' = .05, 'x3' = .5, 'x5' = 1, 'x4' = .02, 'x1'=0.75)
#muvec.raw <- c(0.1, 0.03, 0.2, 0, 1, 4)
#splt.vars.raw <- c('x1', 'x3', 'x5','x4','x1')
#splt.vars.raw <- c(1, 3, 5, 4, 1)
#dvec <- c(1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2)
#xmat <- rbind(xvec1, xvec2, xvec3, xvec4, xvec5, xvec6)
#colnames(xmat) <- c('x1', 'x2', 'x3', 'x4', 'x5')
#FittedValue(xmat, splt.vals.raw, splt.vars.raw, muvec.raw, dvec)
