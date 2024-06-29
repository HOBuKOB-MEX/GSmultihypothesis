source("./GSmulti/GSmulti_discr.R")
size.nbinom = 1 # default value of size parameter for negative binomial

rn <- function(hyp, M)
  rnbinom(1, M * size.nbinom, 1 - exp(hyp) / 2)

pmf <- function(x, M, hyp) dnbinom(x, M * size.nbinom, 1 - exp(hyp) / 2)

lr <- function(s, th, M) {
  exp(th * (s) + log((2 - exp(th))) * M * size.nbinom)
}

weighted <- function(s, gam, thgam, M) {
  return(sum(gam * exp(thgam * (s) + log((2 - exp(thgam))) * M * size.nbinom)))
}
