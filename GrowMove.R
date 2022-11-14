
GrowMove <- function(old_tree, X) {
  ## old_tree is a list with dvec, splitting values, and splitting variables
  terminal_nodes <- which(old_tree$dvec==2)
  ## avoid growing more than d = 3
  terminal_nodes_d <- terminal_nodes[terminal_nodes<8]
  if(any(terminal_nodes_d) == FALSE) 
    print('tree can not grow')
  else{
  grow_node <- sample(terminal_nodes_d, size=1)
  new_var <- sample(colnames(X), size=1)
  
  candidate_splitval <- sort(unique(X[,new_var]))
  ww <- table(X[,new_var])/nrow(X)
  new_split <- sample(candidate_splitval, size=1, prob=ww)
  
  ## set up new tree splitting variables
  new.splt.vars <- c(old_tree$splt.vars[1:(grow_node-1)], new_var, 
                     old_tree$splt.vars[grow_node:length(old_tree$splt.vars)])
  ## set up new tree splitting values
  new.splt.vals <- c(old_tree$splt.vals[1:(grow_node-1)], new_split, 
                     old_tree$splt.vals[grow_node:length(old_tree$splt.vals)])
  names(new.splt.vals) <- c(names(old_tree$splt.vals[1:(grow_node-1)]), new_var, 
                            names(old_tree$splt.vals[grow_node:length(old_tree$splt.vals)]))
  ## set up new tree dvec
  new.dvec <- c(old_tree$dvec[1:(grow_node-1)], 1, old_tree$dvec[(grow_node+1):((grow_node*2)-1)], 2, 2, 
                old_tree$dvec[((grow_node*2)+2):length(old_tree$dvec)])
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
muvec.raw <- c(0.1, 0.03, 0.2, 0, 1, 4)
splt.vars.raw <- c('x1', 'x3', 'x5','x4','x1')
dvec <- c(1, 1, 1, 2, 1, 2, 1, 0, 0, 2, 2, 0, 0, 2, 2)
##
old_tree <- list(dvec = dvec, splt.vars = splt.vars.raw, splt.vals = splt.vals.raw)
X <- rbind.data.frame(xvec1, xvec5, xvec3, xvec6)
colnames(X) <- c('x1', 'x2', 'x3', 'x4', 'x5')

GrowMove(old_tree, X)
