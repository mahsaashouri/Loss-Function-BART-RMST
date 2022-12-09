
Transition_Prob <- function(old_tree, new_tree, X, m){

  if(m==1){
    terminal_nodes <- which(old_tree$dvec==2)
    ## avoid growing more than d = 4
    terminal_nodes_d <- terminal_nodes[terminal_nodes<8]

    ## match splitting values and variables in old and new trees
    # matchvalue <- match(c(unname(old_tree$splt.vals)), new_tree$splt.vals)
    # diffvalue <- setdiff(1:length(new_tree$splt.vals), matchvalue)
    internal_nodes <- which(old_tree$dvec==1)
    internal_nodes_new <- which(new_tree$dvec==1)

    matchvalue <- match(internal_nodes, internal_nodes_new)
    diffvalue <- setdiff(1:length(internal_nodes_new), matchvalue)
    ## find difference between variables and values
    new_val <- new_tree$splt.vals[diffvalue]
    new_var <- new_tree$splt.vars[diffvalue]

    ## compute prob. of all possible values
    candidate_splitval <- sort(unique(X[,new_var]))
    ww <- table(X[,new_var])/nrow(X)

    ## find the desired prob. for the added value
    #PC <- as.numeric(ww[names(ww)==new_val])
    PC <- as.numeric(na.omit(as.numeric((ww[which(X[,new_var]==new_val)]))))

    PMH <- (1/length(terminal_nodes_d)) * (1/ncol(X)) * PC
  }
  if(m==2){
    terminal_nodes <- which(old_tree$dvec==2)
    internal_nodes <- floor(terminal_nodes/2)
    internal_nodes_dup <- duplicated(internal_nodes)

    ## find only internal nodes with two terminal nodes
    internal_nodes_p <- unique(internal_nodes[internal_nodes_dup])

    PMH <- 1/length(internal_nodes_p)
  }
  if(m==3){
    internal_nodes <- which(old_tree$dvec==1)

    ## match splitting values and variables in old and new trees
    matchvalue <- match(c(unname(old_tree$splt.vals)), new_tree$splt.vals)
    diffvalue <- setdiff(1:length(new_tree$splt.vals), matchvalue)

    ## find difference between variables and values
    new_val <- new_tree$splt.vals[diffvalue]
    new_var <- new_tree$splt.vars[diffvalue]

    ## compute prob. of all possible values
    candidate_splitval <- sort(unique(X[,new_var]))
    ww <- table(X[,new_var])/nrow(X)

    ## find the desired prob. for the changed value
    #PC <- as.numeric(ww[names(ww)==new_val])
    PC <- as.numeric(na.omit(as.numeric((ww[which(X[,new_var]==new_val)]))))

    PMH <- (1/length(internal_nodes)) * (1/ncol(X)) * PC
  }
  return(PMH)
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
dvec <- c(1, 1, 1, 2, 1, 2, 1, 0, 0, 2, 2, 0, 0, 2, 2)
#dvec <- c(2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
##
old_tree <- list(dvec = dvec, splt.vars = splt.vars.raw, splt.vals = splt.vals.raw)
X <- rbind.data.frame(xvec1, xvec5, xvec3, xvec6)
colnames(X) <- c('x1', 'x2', 'x3', 'x4', 'x5')

## m = 1 (GrowMove)
new_tree <- list(dvec = c(1, 1, 1, 2, 1, 1, 1, 0, 0, 2, 2, 2, 2, 2, 2), splt.vars = c("x1", "x3", "x5", "x4", "x5", "x1"),
                 splt.vals = c(0.05, 0.50, 1.00, 0.02, 5.00, 0.75))
Transition_Prob(old_tree, new_tree, X = X, m = 1)
## m = 2 (ProneMove)
new_tree <- list(dvec = c(1, 1, 1, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 2, 2), splt.vars = c("x1", "x3", "x5", "x1"),
                 splt.vals = c(0.05, 0.50, 1.00, 0.75))
Transition_Prob(old_tree, new_tree, X = X, m = 2)
## m = 3 (ChangeMove)
new_tree <- list(dvec = c(1, 1, 1, 2, 1, 2, 1, 0, 0, 2, 2, 0, 0, 2, 2), splt.vars = c("x1", "x3", "x5", "x4", "x5"),
                 splt.vals = c(0.05, 0.50, 1.00, 0.02, 0.00))
Transition_Prob(old_tree, new_tree, X = X, m = 3)

