
## Randomly pick an internal node and randomly reassign it a splitting rule

ChangeMove <- function(old_tree, X) {
  # old_tree is a list with dvec, splitting values, and splitting variables
  internal_nodes <- which(old_tree$dvec==1)
  if(any(internal_nodes) == FALSE) 
   print('tree cannot change')
  else{
    if(length(internal_nodes) == 1)
      change_node <- internal_nodes
    else
      change_node <- sample(internal_nodes, size=1)
  ## sample new splitting variable
  new_var <- sample(colnames(X), size=1)
  candidate_splitval <- sort(unique(X[,new_var]))
  ww <- table(X[,new_var])/nrow(X)
  new_split <- sample(candidate_splitval, size=1, prob=ww)
  
  ## update splitting variables
  splt.vars <- replace(old_tree$dvec, old_tree$dvec!=1, NA)
  splt.vars[ which(!is.na(splt.vars))] <- old_tree$splt.vars
  
  ## update splitting values
  splt.vals <- replace(old_tree$dvec, old_tree$dvec!=1, NA)
  splt.vals[ which(!is.na(splt.vals))] <-  old_tree$splt.vals
  
  ## change to new splitting variables
  new.splt.vars <- replace(splt.vars, change_node, new_var)
  new.splt.vars <-  new.splt.vars[!is.na(new.splt.vars)]

  ## change to new splitting values
  new.splt.vals <- replace(splt.vals, change_node, new_split)
  new.splt.vals <-  new.splt.vals[!is.na(new.splt.vals)]
  
  ## change to new dvec
  new.dvec <- old_tree$dvec
  
  return(list(dvec = new.dvec, splt.vars = new.splt.vars, splt.vals = new.splt.vals))
  }
}
  

## Example
xvec1 <- c('x1'=0, 'x2'=0, 'x3'=0, 'x4'=0, 'x5'=0) ## should be assigned to node 4
xvec2 <- c('x1'=0, 'x2'=0, 'x3'=1, 'x4'=0, 'x5'=0) ## should be assigned to node 5
xvec3 <- c('x1'=1, 'x2'=0, 'x3'=1, 'x4'=0, 'x5'=0) ## should be assigned to node 12
xvec4 <- c('x1'=1, 'x2'=0, 'x3'=1, 'x4'=1, 'x5'=0) ## should be assigned to node 13
xvec5 <- c('x1'=0.6, 'x2'=0, 'x3'=1, 'x4'=1, 'x5'=5) ## should be assigned to node 14
xvec6 <- c('x1' = 1, 'x2' = 2, 'x3' = 3, 'x4' = 10, 'x5' = 5) ## should be assigned to node 15
splt.vals.raw <- c('x1' = .05, 'x3' = .5, 'x5' = 1, 'x4' = .02, 'x1'=0.75)
splt.vars.raw <- c('x1', 'x3', 'x5','x4','x1')
dvec <- c(1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2)
##
old_tree <- list(dvec = dvec, splt.vars = splt.vars.raw, splt.vals = splt.vals.raw)
X <- rbind.data.frame(xvec1, xvec5, xvec3, xvec6)
colnames(X) <- c('x1', 'x2', 'x3', 'x4', 'x5')

ChangeMove(old_tree, X)


