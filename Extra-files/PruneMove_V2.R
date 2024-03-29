
## Randomly pick a parent of two terminal nodes and turn it into a terminal node by collapsing the nodes below it

GetChildIND <- function(x){
  d <- floor(log2(x))+1
  node_vec <- c()
  if(d<4){
    for(i in 1:(4-d)){
      for(j in 0:((2^i)-1)){
        node_vec <- append(node_vec, ((2^i)*x)+j, after = length(node_vec))
      }
    }
  }
  return(node_vec)
}

PruneMove <- function(old_tree, X) {
  ## old_tree is a list with dvec, splitting values, and splitting variables
  internal_nodes <- which(old_tree$dvec==1)
  ## avoid pruning tree with one node
  if(any(internal_nodes) == FALSE)
    print('tree cannot prone')
  else{
    prune_node <- sample(internal_nodes, size=1)

    ## All the possible child nodes for sampled node number
    prune_node_all <- c(prune_node, GetChildIND(prune_node))
    remove_internals <- prune_node_all[prune_node_all<8]

    ## node numbers need to be removed
    remove_vars <- remove_internals[remove_internals %in% internal_nodes]
    #remove_vars <- prune_node_all[!is.na(prune_node_all)]

    ## update splitting variables
    splt.vars <- replace(old_tree$dvec, old_tree$dvec!=1, NA)
    splt.vars[ which(!is.na(splt.vars))] <- old_tree$splt.vars

    ## update splitting values
    splt.vals <- replace(old_tree$dvec, old_tree$dvec!=1, NA)
    splt.vals[ which(!is.na(splt.vals))] <-  old_tree$splt.vals

    ## set up new splitting variables
    new.splt.vars <- splt.vars[-remove_vars]
    new.splt.vars <-  new.splt.vars[!is.na(new.splt.vars)]

    ## set up new splitting values
    new.splt.vals <- splt.vals[-remove_vars]
    new.splt.vals <-  new.splt.vals[!is.na(new.splt.vals)]

    ## set up new dvec
    new.dvec <- old_tree$dvec
    new.dvec[prune_node_all] <- 0; new.dvec[prune_node] <- 2

    return(list(dvec = new.dvec, splt.vars = new.splt.vars, splt.vals = new.splt.vals))
  }
}

## Example
#xvec1 <- c('x1'=0, 'x2'=0, 'x3'=0, 'x4'=0, 'x5'=0) ## should be assigned to node 4
#xvec2 <- c('x1'=0, 'x2'=0, 'x3'=1, 'x4'=0, 'x5'=0) ## should be assigned to node 5
#xvec3 <- c('x1'=1, 'x2'=0, 'x3'=1, 'x4'=0, 'x5'=0) ## should be assigned to node 12
#xvec4 <- c('x1'=1, 'x2'=0, 'x3'=1, 'x4'=1, 'x5'=0) ## should be assigned to node 13
#xvec5 <- c('x1'=0.6, 'x2'=0, 'x3'=1, 'x4'=1, 'x5'=5) ## should be assigned to node 14
#xvec6 <- c('x1' = 1, 'x2' = 2, 'x3' = 3, 'x4' = 10, 'x5' = 5) ## should be assigned to node 15
#splt.vals.raw <- c('x1' = .05, 'x3' = .5, 'x5' = 1, 'x4' = .02, 'x1'=0.75)
#splt.vars.raw <- c('x1', 'x3', 'x5','x4','x1')
#dvec <- c(1, 1, 1, 2, 1, 2, 1, 0, 0, 2, 2, 0, 0, 2, 2)
##
#old_tree <- list(dvec = dvec, splt.vars = splt.vars.raw, splt.vals = splt.vals.raw)
#X <- rbind.data.frame(xvec1, xvec5, xvec3, xvec6)
#colnames(X) <- c('x1', 'x2', 'x3', 'x4', 'x5')

#PruneMove(old_tree, X)

