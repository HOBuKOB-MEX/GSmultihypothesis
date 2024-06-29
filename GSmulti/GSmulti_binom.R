source("./GSmulti/GSmulti_discr.R")

rn <- function(hyp, M)
  rbinom(1, M, exp(hyp) / (1 + exp(hyp)))

pmf <- function(x, M, th) dbinom(x, M, exp(th) / (1 + exp(th)))

lr <- function(s, th, M) {
  exp(th * s - log((1 + exp(th)) / 2) * M)
}

weighted <- function(s, gam, thgam, M) {
  return(sum(gam * exp(thgam * s - log((1 + exp(thgam)) / 2) * M)))
}
