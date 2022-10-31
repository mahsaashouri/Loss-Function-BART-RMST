

## Inputs: Predictors, Splitting values and variables, mu vectors , and tree structure (from 26 possible structures)
## Output: Fitted mu value

LeftChild <- function(x){
  y <- c(2, 4, 6, 8, 10, 12, 14)
  return(y[x])
}
RightChild <- function(x){
  y <- c(3, 5, 7, 11, 13, 15)
  return(y[x])
}

FittedValue <- function(xvec, splt.vals, splt.vars.raw, muvec.raw, dvec){
  if(dvec[1] == 2){
    non_terminal <- FALSE
    mu.hat <- muvec.raw
  }
  else{ 
    non_terminal <- TRUE
    current.node <- 1
    ## update mu vector
    muvec <- replace(dvec, dvec!=2, NA)
    muvec[ which(!is.na(muvec))] <- muvec.raw
    ## update splitting variable vector
    splt.vars <- replace(dvec, dvec!=1, NA)
    splt.vars[ which(!is.na(splt.vars))] <- splt.vars.raw
    current.var <- splt.vars[1]
  }
  while(non_terminal){
    if(xvec[current.var] <= splt.vals[current.var]){
      new_node <- LeftChild(current.node)
      new.var <- splt.vars[new_node] 
    }
      else{
        new_node <- RightChild(current.node)
        new.var <- splt.vars[new_node] 
      }
    if(dvec[new_node] == 2){
      non_terminal <- FALSE
      mu.hat <- muvec[new_node]
    }
    else{
      current.var <- new.var
      current.node <- new_node
    }
  }
  return(mu.hat)
}

## Example
xvec <- c('x1' = 1, 'x2' = 2, 'x3' = 3, 'x4' = 10, 'x5' = 5)
splt.vals <- c('x1' = .05, 'x3' = .5, 'x5' = 1, 'x4' = .02)
muvec.raw <- c(0.1, 0.03, 0.2, 0, 1, 4)
splt.vars.raw <- c('x1', 'x3', 'x5','x4')
dvec <- c(1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2)
FittedValue(xvec, splt.vals, splt.vars.raw, muvec.raw, dvec)

