LogSumExp <- function(lseq) {
    ## This function computes
    ##.   sum(exp(lseq))

    cstar <- max(lseq)
    ans <- exp(cstar)*sum(exp(lseq - cstar))
    return(ans)
}





LambdaCDF <- function(t, nevent, nrisk, c0, delta_alpha) {

  even_seq <- seq(0, nevent, by=2)
  prob_even <- (nrisk - nevent)/(nrisk - nevent + c0 + even_seq)

  cstar_even <- max(lchoose(nevent, even_seq) + c0*delta_alpha*log(prob_even))
  even_seq <- lchoose(nevent, even_seq) + c0*delta_alpha*log(prob_even) - cstar_even
  sum_even <- exp(cstar_even)*sum( exp(even_seq)   )

  if(nevent > 0) {
    odd_seq <- seq(1, nevent, by=2)
    prob_odd <- (nrisk - nevent)/(nrisk - nevent + c0 + odd_seq)

    cstar_odd <- max(lchoose(nevent, odd_seq) + c0*delta_alpha*log(prob_odd))
    odd_seq <- lchoose(nevent, odd_seq) + c0*delta_alpha*log(prob_odd) - cstar_odd
    sum_odd <- exp(cstar_odd)*sum( exp(odd_seq)   )
  }


  if(nevent > 0) {
    *(sum_even - sum_odd)

  } else {
    *sum_even
  }
  tmp <- (((-1)^(0:nevent))*choose(nevent, 0:nevent))/((nrisk - nevent + c0 + 0:nevent)^(c0*delta_alpha))
  ww <- tmp/sum(tmp)
  print(tmp)
  ff <- pgamma(t, shape=c0*delta_alpha, rate=nrisk - nevent + c0 + 0:nevent)
  FF <- sum(ww*pgamma(t, shape=c0*delta_alpha, rate=nrisk - nevent + c0 + 0:nevent))
  return(FF)
}

InvLambdaCDF <- function(u, nevent, nrisk, c0, delta_alpha) {
  ff <- function(tt, u, nevent, nrisk, c0, delta_alpha) {
    return( LambdaCDF(tt, nevent, nrisk, c0, delta_alpha) - u)
  }
  #ulimit <- abs(qgamma(u, shape=c0*delta_alpha, rate=nrisk - nevent + c0 + 0:nevent))
  ## need a better way to figure out ulimit
  ulimit <- 100
  print(c(nrisk, nevent))
  print(LambdaCDF(0, nevent, nrisk, c0, delta_alpha))
  print(LambdaCDF(100, nevent, nrisk, c0, delta_alpha))

  ans <- uniroot(ff, lower=0, upper=ulimit,
                 u = u, nevent=nevent, nrisk=nrisk, c0=c0, delta_alpha=delta_alpha)$root
  return(ans)
}
DrawLambdas <- function(nevents, nrisks, c0, delta_alpha) {
  nbins <- length(nevents)
  lambda_draw <- rep(NA, nbins)
  for(j in 1:nbins) {
    u <- runif(1)
    lambda_draw[j] <- InvLambdaCDF(u, nevents[j], nrisks[j], c0, delta_alpha)
  }
  return(lambda_draw)
}

DrawIPCW <- function(U, delta, Utau, sgrid, c0, delta_alpha) {
  J <- length(sgrid)
  E <- R <- rep(NA, J-1)
  for(j in 2:J) {
    E[j-1] <- sum((1 - delta)*(U > sgrid[j-1])*(U <= sgrid[j]))
    R[j-1] <- sum(U > sgrid[j-1])
  }
  lambda.draw <- DrawLambdas(nevents=E, nrisks=R, c0=c0, delta_alpha=delta_alpha)
  CumHazFn <- approxfun(sgrid, cumsum(c(0, lambda.draw)))
  weight_ans <- exp(CumHazFn(Utau))
  return(weight_ans)
}


