source("./GSmulti/GSmulti_discr.R")

rn <- function(hyp, M)
  rpois(1, M * exp(hyp))

pmf <- function(x, M, th) dpois(x, M * exp(th))

lr <- function(s, th, M) {
  exp(th * s - (exp(th) - 1) * M)
}

weighted <- function(s, gam, thgam, M) {
  return(sum(gam * exp(thgam * s - (exp(thgam) - 1) * M)))
}
