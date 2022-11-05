
## RMST regression loss function
## U1, ..., Un: times, delta1, ..., deltan: event indicators, Yi = min(Ui, tau)
## G is the survival function of the censoring time distribution, i.e., G(t) = P {Ci > t}

RMST_Loss <- function(U, tau, delta, Gvec, fitted_vector, J){
  #c <- c(); cvec <- c(); N <- c(); Nvec <- c();
  # for(j in 1:J){
  #   for(i in 1:n){
  #     c[i] <- (1 - delta[i])*ifelse(U[i]>=(j-1) && U[i]<j, 1, 0) 
  #     N[i] <- ifelse(U[i]>(j-1), 1, 0)
  #   }
  #   cvec[j] <- sum(c)
  #   Nvec[j] <- sum(N)
  # }
  # n <- length(U);g <- c();loss.c <- c();
  # for( j in 1:J){
  #   g[j] <- Gvec[j-1] - Gvec[j]
  #   loss.c[j] <- (cvec[j]*log(g[j])) + ((Nvec[j]-cvec[j])*log(1-g[j]))
  # }
  logY <- c(); loss.RMST <- c();
  for(i in 1:n){
    loss.RMST[i] <- (delta[i]/Gvec(U[i]))*(log(U[i])-fitted_vector[i])^2 #+ sum(loss.c)
  }
  return(sum(loss.RMST))
}

