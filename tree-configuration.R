
FindDvec <- function(tree_structure) {
    ## All possible trees - 2:terminals; 1: non-terminals; 0: otherwise
    dlist <- list()
    ## d = 1
    dlist[[1]] <- c(2, rep(0,14))
    ## d = 2
    dlist[[2]] <- c(1, 2, 2, rep(0, 12))
    dlist[[3]] <- c(1, 1, 2, 2, 2, rep(0, 10))
    dlist[[4]] <- c(1, 2, 1, 0, 0, 2, 2, rep(0, 8))
    dlist[[5]] <- c(1, 1, 1, 2, 2, 2, 2, rep(0, 8))
    ## d= 3
    dlist[[6]] <- c(1, 1, 2, 1, 2, 0, 0, 2, 2, 0, 0, rep(0, 4))
    dlist[[7]] <- c(1, 1, 2, 2, 1, 0, 0, 0, 0, 2, 2, rep(0, 4))
    dlist[[8]] <- c(1, 1, 2, 1, 1, 0, 0, 2, 2, 2, 2, rep(0, 4))
    ## d = 4
    dlist[[9]] <- c(1, 2, 1, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 2, 2)
    dlist[[10]] <- c(1, 2, 1, 0, 0, 1, 2, 0, 0, 0, 0, 2, 2, 0, 0)
    dlist[[11]] <- c(1, 2, 1, 0, 0, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2)
    ## d = 5
    dlist[[12]] <- c(1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0)
    dlist[[13]] <- c(1, 1, 1, 2, 1, 2, 2, 0, 0, 2, 2, 0, 0, 0, 0)
    dlist[[14]] <- c(1, 1, 1, 2, 2, 1, 2, 0, 0, 0, 0, 2, 2, 0, 0)
    dlist[[15]] <- c(1, 1, 1, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 2, 2)
    dlist[[16]] <- c(1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2)
    dlist[[17]] <- c(1, 1, 1, 2, 1, 2, 1, 0, 0, 2, 2, 0, 0, 2, 2)
    dlist[[18]] <- c(1, 1, 1, 2, 1, 1, 2, 0, 0, 2, 2, 2, 2, 0, 0)
    dlist[[19]] <- c(1, 1, 1, 1, 2, 2, 1, 2, 2, 0, 0, 0, 0, 2, 2)
    dlist[[20]] <- c(1, 1, 1, 1, 2, 1, 2, 2, 2, 0, 0, 2, 2, 0, 0)
    dlist[[21]] <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0)
    dlist[[22]] <- c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 0, 0)
    dlist[[23]] <- c(1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 0, 0, 2, 2)
    dlist[[24]] <- c(1, 1, 1, 1, 2, 1, 1, 2, 2, 0, 0, 2, 2, 2, 2)
    dlist[[25]] <- c(1, 1, 1, 2, 1, 1, 1, 0, 0, 2, 2, 2, 2, 2, 2)
    dlist[[26]] <- c(rep(1, 7), rep(2, 8))
    return(dlist[[tree_structure]])
}
#Dmat <- matrix(NA, nrow = 26, ncol = 15)
#for (k in 1:26) {
#  Dmat[k,] <- FindDvec(k)
#}
