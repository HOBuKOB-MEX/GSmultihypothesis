library(rootSolve)
source("./GSmulti/GSmulti_cont.R")

lr <- function(s, th, M) {
  exp(th * s - th^2 / 2 * M)
}

sum.cdf <- function(pt, x, th, M)
  pnorm(pt, mean = th * M + x, sd = sqrt(M))

sum.pdf <- function(pt, th, M)
  dnorm(pt, mean = th * M, sd = sqrt(M))

intgr <- function(spts, vls, t, M) {
  # Backward induction integration
  # numerical integral of \int u(x)*dnorm(x-t,sd=sqrt(M)) dx
  # spts is grid of x values, vls corresponding values of u
  PhiVec = pnorm(spts - t, sd = sqrt(M))
  phiVec = dnorm(spts - t, sd = sqrt(M))
  diffPhi = diff(PhiVec)
  diffphi = diff(phiVec)
  tmp = sum(head(vls, -1) * diffPhi)
  vec = diffPhi * head(spts - t, -1) + M * diffphi
  tmp = tmp - sum((diff(vls) / diff(spts)) * vec)
  return(tmp)
}

rn <- function(hyp, M)   rnorm(1, mean = hyp * M, sd = sqrt(M))

dividing <- function(n, l, th, M) {
  a = array(dim = length(th) - 1)
  for (i in 1:(length(th) - 1)) {
    a[i] = uniroot(function(x) {
      a = comp(x, n, i, l, th, M)
      b = comp(x, n, i + 1, l, th, M)
      return((a - b)) }, interval = c(-1, 1), extendInt = "upX")$root
  }
  return(a)
}
