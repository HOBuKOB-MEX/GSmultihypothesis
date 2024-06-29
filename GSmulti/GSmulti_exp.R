library(rootSolve)
library(Bolstad2)
source("./GSmulti/GSmulti_cont.R")
Interval.width = 500

# th is a parameter of exponential distribution in canonical parametrization
#  th is equal to minus the rate of the exponential distribution
#  parametrized exponential family representation:
#  exp(th*x+log(1-th))*h(x) (th negative, h standard exponential density)
lr <- function(s, th, M) {
  exp(th * s + log(1 - th) * M)
}

sum.cdf <- function(pt, x, th, M)
  pgamma(pt - x, rate = (1 - th), shape = M)

sum.pdf <- function(pt, th, M)
  dgamma(pt, rate = (1 - th), shape = M)

intgr <- function(spts, vls, t, M) {
  tmp = sintegral(spts, vls * dgamma(spts - t, rate = 1, shape = M), n.pts = length(spts))$int
  return(tmp)
}

rn <- function(hyp, M) rgamma(1, rate = (1 - hyp), shape = M)

dividing <- function(n, l, th, M) {
  a = array(dim = length(th) - 1)
  for (i in 1:(length(th) - 1)) {
    a[i] = tail(uniroot.all(function(x) {
      a = comp(x, n, i, l, th, M)
      b = comp(x, n, i + 1, l, th, M)
      return((a - b) / (a + b)) }, interval = c(0, Interval.width)), 1) }
  return(a)
}
