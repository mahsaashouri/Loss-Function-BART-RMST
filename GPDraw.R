GPDraw <- function(eU, sgrid, num.risk, num.events, kappa0) {
  # eU - the vector of observed-follow up times (exp of input to RMST_BCART)
  # sgrid - grid of J+1 points of the form (s_{j-1}, s_{j}] for
  #          the grouped likelihood function
  # num.risk - vector of length J; number at risk at time s_{j-1}
  # num.events - vector of length J; number of events in interval (s_{j-1}, s_{j}]
  # kappa0 - value of kappa0 (like a "mass" parameter for the GP prior)
  J <- length(num.risk)
  #lambda.exp <- sum(num.events)/sum(obs.times) # "Empirical Bayes" choice of base cum. hazard
  #alpha.seq <- c0*lambda.exp*event.times
  #alpha.diff <- c(alpha.seq[1], diff(alpha.seq))
  alpha.diff <- rep(0.01, J)

  lambda.vec <- rgamma(J, shape=alpha.diff + num.events, rate=kappa0 + num.risk - num.events)
  Lambda.vec <- c(0, cumsum(lambda.vec))
  LamFn <- approxfun(sgrid, Lambda.vec)
  GPweights <- LamFn(eU)

  return(GPweights)
}

ComputeSurvStatistics <- function(sgrid, times, status) {
  J <- length(sgrid) - 1
  n.risk <- n.event <- rep(NA, J)
  for(j in 1:J) {
    n.risk[j] <- sum(times >= sgrid[j])
    n.event[j] <- sum(times > sgrid[j] & times <= sgrid[j+1] & status==1)
  }
  return(list(n.risk=n.risk, n.event=n.event))
}

## Example Run:
#library(survival)
#data(cancer)
#sgrid <- seq(0, 500, by=50)
#SS <- ComputeSurvStatistics(sgrid, times=cancer$time,
#                      status=as.numeric(cancer$status==1))

#lam.draw <- GPDraw(U=cancer$time, sgrid=sgrid, num.risk=SS$n.risk,
#                   num.events=SS$n.event, kappa0=1)

#SS <- ComputeSurvStatistics(sgrid=sgrid, times=exp(Y), status=1 - delta)
#
#lam.draw <- GPDraw(U=U, sgrid=sgrid, num.risk=SS$n.risk,
#                   num.events=SS$n.event, kappa0=1)
#Gvec <- exp(-lam.draw)
