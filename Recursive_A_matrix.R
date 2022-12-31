
FindParent <- function(x){
  return(floor(x/2))
}

UpdateA.rec <- function(Amat, splt.vars, r, c){
  p <- FindParent(r)
  if(p > 0){
    i.name = splt.vars[p]
    i = as.numeric(substr(i.name, 2, nchar(i.name)))
    Amat[i, c] = 1
    Amat = UpdateA.rec(Amat, splt.vars, p, c)
  }
  return(Amat)
}



Amatrix <- function(xvec, splt.vars.raw, dvec){

  Amat <- matrix(0, nrow = length(xvec), ncol = length(which(dvec == 2)))
  terminal_indexes = which(dvec == 2)

  ## update splitting variable vector
  splt.vars <- replace(dvec, dvec!=1, NA)
  splt.vars[ which(!is.na(splt.vars))] <- splt.vars.raw

  for(i in 1:ncol(Amat)){
    Amat =  UpdateA.rec(Amat, splt.vars, terminal_indexes[i], i)
  }
  return(Amat)
}



##Example

#xvec6 <- c('x1' = 1, 'x2' = 2, 'x3' = 3, 'x4' = 10, 'x5' = 5) ## should be assigned to node 15
##splt.vals.raw <- c('x1' = .05, 'x3' = .5, 'x5' = 1, 'x4' = .02, 'x1'=0.75)
#splt.vars.raw <- c('x1', 'x3', 'x5','x4','x1')
#dvec <- c(1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2)

#Amatrix(xvec6, splt.vars.raw, dvec)

