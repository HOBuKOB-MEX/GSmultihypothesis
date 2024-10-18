library(rootSolve)

monte_carlo_simulation <- function(test, hyp, nMC = 10000) {
  # test simulation
  # hyp = true parameter value
  # returns relative frequency of acceptations
  # the average number of groups and number of observations
  M = test$info$M
  H = length(test$info$M)
  K = length(test$info$th)
  ss = 0
  totaccepted = rep(0, K)
  totn = 0
  totobs = 0
  for (i in 1:nMC) {
    s = 0 # accumulated sum for test run
    n = 0 # number of steps in current run
    runobs = 0
    accepted = rep(0, K) # number of accepting in the run (0 or 1)
    for (stage in 1:(H)) { # run starts
      # generate
      summand = rn(hyp, M[stage])
      s = s + summand # accumulated sum
      n = n + 1  #one step more
      plan = test$data[[stage]]$plan
      runobs = runobs + M[stage]
      if (stage == H) {
        accepted[plan[s + 1]] = accepted[plan[s + 1]] + 1
        break
      }
      stopcond = plan[s + 1] != 0
      if (stopcond)
        accepted[plan[s + 1]] = accepted[plan[s + 1]] + 1
      if (stopcond | stage == H)
        break
    }
    totaccepted = totaccepted + accepted
    totn = totn + n
    totobs = totobs + runobs
    ss = ss + n^2
  }
  OC = (totaccepted) / as.double(nMC)
  totgroups = totn / as.double(nMC)
  obs = totobs / as.double(nMC)
  return(list(OC = OC, groups = totgroups, ESS = obs))
}

comp <- function(s, n, num, l, th, M) {
  k = length(th)
  res = 0
  for (i in 1:k) { if (i == num)next
    res = res + l[i, num] * lr(s, th[i], cumsum(M)[n])
  }
  return(res)
}

dividing <- function(n, l, th, M) {
  a = sapply(0:(maxval * n), function(s) {
    which.min(sapply(1:length(th), function(num) comp(s, n, num, l, th, M))) })
  b = sapply(0:(maxval * n),
             function(s) { min(sapply(1:length(th), function(num) comp(s, n, num, l, th, M))) })
  return(rbind(a, b))
}

back <- function(s, n, M, previous) {
  tmp = 0
  for (i in 0:maxval) {
    prod = previous[2, s + i + 1] * pmf(i, M[n], 0)
    tmp = tmp + prod
  }
  return(tmp + cost_fn(M[n]) * weighted(s, gam, thgam, ifelse(n > 1, cumsum(M)[n - 1], 0)))
}

##########################################
OptPlan <- function(M, l, th, gam, thgam, cost_fn, maxval) {
  ##########################################


  step <- function(n, l, th, M, previous) {
    present = dividing(n - 1, l, th, M)
    for (i in 0:(length(present[1,]) - 1)) {
      tmp = back(i, n, M, previous)
      if (present[2, i + 1] > tmp) {
        present[2, i + 1] = tmp
        present[1, i + 1] = 0
      }
    }
    return(present)
  }

  n = length(M)
  previous = dividing(n, l, th, M)
  test = list()
  if (n > 1)
    repeat {
      test[[n]] = list(n = n, plan = previous[1,])

      previous = step(n, l, th, M, previous)
      if (n == 2)break
      n = n - 1
    }
  test[[1]] = list(n = 1, plan = previous[1,])
  Lagr = (back(0, 1, M, previous))
  return(list(data = test, info = list(M = M, lam = l, th = th, gam = gam, thgam = thgam, maxval = maxval, Lagrangian = Lagr)))
}


PAccept <- function(test, th0, fin) {
  ##########################################

  back <- function(s, n, M, previous) {
    tmp = 0
    for (i in 0:maxval) {
      prod = previous[s + i] * pmf(i, M[n], th0)
      tmp = tmp + prod }
    return(tmp)
  }

  step <- function(n, M, previous) {
    present = test$data[[n - 1]]$plan
    val = array(0, dim = length(present))
    val[which(present == fin)] = 1

    for (i in 1:(length(present))) {
      if (present[i] == 0)
        val[i] = back(i, n, M, previous)
    }
    return(val)
  }

  M = test$info$M
  n = length(M)
  previous = test$data[[n]]$plan
  val = rep(0, length(previous))
  for (i in 1:length(previous)) {
    if (previous[i] == fin)
      val[i] = 1
  }
  if (n > 1)
    repeat {
      val = step(n, M, val)
      if (n == 2)break
      n = n - 1
    }
  return(back(1, 1, M, val))
}

mincomp <- function(s, n, l, th, M) {
  k = length(th)
  vec = comp(s, n, 1, l, th, M)
  for (i in 2:k)
    vec = pmin(comp(s, n, i, l, th, M), vec)
  return(vec)
}

##########################################################
ESS <- function(test, th0, cost_fn) {

  back <- function(s, n, M, previous) {
    tmp = 0
    for (i in 0:maxval) {
      prod = previous[s + i] * pmf(i, M[n], th0)
      tmp = tmp + prod }
    return(tmp)
  }

  step <- function(n, M, previous) {
    present = test$data[[n - 1]]$plan
    val = array(0, length(present))
    for (i in 1:(length(present))) {
      if (present[i] == 0)
        val[i] = cost_fn(M[n]) + back(i, n, M, previous)
    }
    return(val)
  }

  M = test$info$M
  n = length(M)
  previous = test$data[[n]]$plan
  val = array(0, length(previous))
  if (n > 1)
    repeat {
      val = step(n, M, val)
      if (n == 2)break
      n = n - 1
    }
  return(cost_fn(M[1]) + back(1, 1, M, val))
}

