monte_carlo_simulation <- function(test, hyp,
                                   nMC = 10000)
  # Monte Carlo simulation of a test
{
  M = test$info$M
  H = length(test$info$M)
  K = length(test$info$th)
  ss = 0
  totaccepted = rep(0, K)
  totn = 0
  totobs = 0
  # runobs=0
  for (i in 1:nMC) {
    s = 0 #accumulated sum for test run
    n = 0 # number of steps in current run
    runobs = 0
    accepted = rep(0, K) #number of accepting in the run (0 or 1)
    for (stage in 1:(H)) { #run starts
      #generate
      summand = rn(hyp, M[stage])
      s = s + summand #accumulated sum
      n = n + 1  #one step more
      runobs = runobs + M[stage]
      nInt = test$data[[stage]]$nInt
      cInt = test$data[[stage]]$cInt
      a = test$data[[stage]]$div
      if (stage == H) {
        for (i in 1:(K - 1)) {
          if (s < a[i]) {
            accepted[i] = accepted[i] + 1
            break
          }
        }
        if (s > a[K - 1])
          accepted[K] = accepted[K] + 1
        break
      }
      stopcond = TRUE
      for (i in 1:nInt) {
        if (s <= cInt[i, 2] & s >= cInt[i, 1]) stopcond = FALSE
      }
      if (stopcond) {
        for (i in 1:(K - 1)) {
          if (s < a[i]) {
            accepted[i] = accepted[i] + 1
            break
          }
        }
        if (s > a[K - 1])
          accepted[K] = accepted[K] + 1
        break
      }
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

makegrid <- function(left, right, gridsz) {
  nint = ceiling((right - left) / gridsz)
  sval = seq(left, right, length.out = nint + 1)
}

weighted <- function(s, gam, thgam, M) {
  return(sum(gam * lr(s, thgam, M)))
}

comp <- function(s, n, num, l, th, M) {
  k = length(th)
  res = rep(0, length(s))
  for (i in 1:k) { if (i == num)next
    res = res + l[i, num] * sum.pdf(s, th[i], cumsum(M)[n])
  }
  return(res)
}

compg <- function(s, n, num, l, th, M) {
  k = length(th)
  res = rep(0, length(s))
  for (i in 1:k) { if (i == num)next
    res = res + l[i, num] * lr(s, th[i], cumsum(M)[n])
  }
  return(res)
}

mincomp <- function(s, n, l, th, M) {
  k = length(th)
  vec = compg(s, n, 1, l, th, M)
  for (i in 2:k)
    vec = pmin(compg(s, n, i, l, th, M), vec)
  return(vec)
}

##############################################################################
OptPlan <- function(M, l, th, gam, thgam, gridsz, cost_fn, margin = 20) {

  exper0 <- function(x, n, l, th, M) {
    # a=test$data[[n]]$div
    a = dividing(n, l, th, M)
    k = length(th)
    tmp = 0
    ext.div = c(-Inf, a, Inf)
    for (j in 1:k) {
      for (i in 1:k) {
        if (j == i)next
        prob = sum.cdf(ext.div[i + 1], x, th[j], M[n]) - sum.cdf(ext.div[i], x, th[j], M[n])
        tmp = tmp + l[j, i] *
          prob *
          lr(x, th[j], ifelse(n > 1, cumsum(M)[n - 1], 0))

      }
    }
    return(tmp + cost_fn(M[n], n) * ifelse(n > 1, weighted(x, gam, thgam, cumsum(M)[n - 1]), 1))
  }


  exper1 <- function(x, n, l, th, M) {

    k = length(th)
    nInt = test[[n]]$nInt
    cInt = test[[n]]$cInt
    # gridlist=test[[n]]$gridlist
    a = test[[n]]$div
    ext.div = c(-Inf, a, Inf)
    cut = ext.div
    s = 0
    from = 1
    to = max(which(ext.div < cInt[1, 1])) + 1
    cut[to] = cInt[1, 1]
    current = from
    repeat {
      for (j in 1:k) {
        if (j != current) s = s + l[j, current] *
          (sum.cdf(cut[current + 1], x, th[j], M[n]) - sum.cdf(cut[current], x, th[j], M[n])) *
          lr(x, th[j], ifelse(n > 1, cumsum(M)[n - 1], 0))

      }
      current = current + 1
      if (current == to)break
    }
    from = min(which(ext.div > cInt[nInt, 2])) - 1
    to = k + 1
    cut = ext.div
    cut[from] = cInt[nInt, 2]
    current = from
    repeat {
      for (j in 1:k) {
        if (j != current) s = s +
          l[j, current] *
            (sum.cdf(cut[current + 1], x, th[j], M[n]) - sum.cdf(cut[current], x, th[j], M[n])) *
            lr(x, th[j], ifelse(n > 1, cumsum(M)[n - 1], 0))
      }
      current = current + 1
      if (current == to)break
    }
    if (nInt > 1)
      for (int in 2:(nInt)) {
        cut = ext.div
        from = max(which(ext.div < cInt[int - 1, 2]))
        to = min(which(ext.div > cInt[int, 1]))
        cut[from] = cInt[int - 1, 2]
        cut[to] = cInt[int, 1]
        current = from
        repeat {
          for (j in 1:k) {
            if (j != current) s = s + l[j, current] *
              (sum.cdf(cut[current + 1], x, th[j], M[n]) - sum.cdf(cut[current], x, th[j], M[n])) *
              lr(x, th[j], ifelse(n > 1, cumsum(M)[n - 1], 0))
          }
          current = current + 1
          if (current == to)break
        }
      }
    for (int in 1:nInt)
      s = s + intgr(gridlist1[[int]]$sval, gridlist1[[int]]$rhoval, x, M[n])
    return(s + cost_fn(M[n], n) * weighted(x, gam, thgam, cumsum(M)[n - 1]))
  }

  eff0 <- function(x, n, l, th, M)mincomp(x, n - 1, l, th, M) - exper0(x, n, l, th, M)
  eff0.v <- function(x, nd, l, th, M)Vectorize(eff0, vectorize.args = "x")(x = x, n = nd, l = l, th = th, M = M)

  eff1 <- function(x, n, l, th, M)mincomp(x, n - 1, l, th, M) - exper1(x, n, l, th, M)
  eff1.v <- function(x, nd, l, th, M)Vectorize(eff1, vectorize.args = "x")(x = x, n = nd, l = l, th = th, M = M)

  ##############################
  n = length(M)
  k = length(th)
  test = list()

  a = dividing(n, l, th, M)
  test[[n]] = list(n = n, div = a)
  if (n > 1) {
    a = (dividing(n - 1, l, th, M))
    cInt = array(dim = c(k - 1, 2))
    gridlist = list()
    left = min(a)
    right = max(a)
    limits = uniroot.all(eff0.v, lower = left - margin, upper = right + margin, nd = n, l = l, th = th, M = M)
    nInt = length(limits) / 2
    current = 1
    cInt = array(dim = c(k - 1, 2))
    for (i in 1:nInt) {
      cInt[i, 1] = limits[current]
      cInt[i, 2] = limits[current + 1]
      current = current + 2
      sval = makegrid(cInt[i, 1], cInt[i, 2], gridsz)
      rhoval = sapply(sval, function(x)(exper0(x, n, l = l, th = th, M = M)))
      gridlist[[i]] = list(sval = sval, rhoval = rhoval)
    }
    test[[n - 1]] = list(n = n - 1, div = a, nInt = nInt, cInt = cInt)
    gridlist1 = gridlist
    n = n - 1
    if (n > 1)
    {
      repeat {
        a = dividing(n - 1, l, th, M)
        gridlist = list()
        left = min(a)
        right = max(a)
        limits = (uniroot.all(eff1.v, lower = left - margin, upper = right + margin, nd = n, l = l, th = th, M = M))
        nInt = length(limits) / 2

        current = 1
        cInt = array(dim = c(k - 1, 2))
        for (i in 1:nInt) {
          cInt[i, 1] = limits[current]
          cInt[i, 2] = limits[current + 1]
          current = current + 2
          sval = makegrid(cInt[i, 1], cInt[i, 2], gridsz)
          rhoval = sapply(sval, function(x)(exper1(x, n, l = l, th = th, M = M)))
          gridlist[[i]] = list(sval = sval, rhoval = rhoval)
        }
        test[[n - 1]] = list(n = n - 1, div = a, nInt = nInt, cInt = cInt)
        gridlist1 = gridlist
        if (n == 2)break
        n = n - 1
      }
    }
  }
  if (length(M) > 1)  Lagr = cost_fn(M[1], 1) + exper1(0, 1, l = l, th = th, M = M)
  else Lagr = exper0(0, 1, l, th, M)
  return(list(data = test, info = list(theta = th, lambda = l, M = M, gridsz = gridsz, vartheta = thgam, gamma = gam, const_fn = cost_fn, Lagr = Lagr)))
}

################################################################
ESS <- function(test, th0, gridsz, cost_fn) {

  asn0 <- function(x, n, th0, M) {
    return(cost_fn(M[n], n) * lr(x, th0, ifelse(n > 1, cumsum(M)[n - 1], 0)))
  }

  asn1 <- function(x, n, th0, M) {
    k = length(test$info$th)
    nInt = test$data[[n]]$nInt
    cInt = test$data[[n]]$cInt
    tmp1 = cost_fn(M[n], n) * lr(x, th0, ifelse(n > 1, cumsum(M)[n - 1], 0))
    for (i in 1:nInt)
      tmp1 = tmp1 + intgr(gridlist1[[i]]$sval, gridlist1[[i]]$rhoval, x, M[n])
    return(tmp1)
  }

  n = length(test$info$M)
  k = length(test$info$th)
  M = test$info$M
  if (n > 1) {
    nInt = test$data[[n - 1]]$nInt
    cInt = test$data[[n - 1]]$cInt
    gridlistnew = list()
    for (i in 1:nInt) {
      sval = makegrid(cInt[i, 1], cInt[i, 2], gridsz)
      rhoval = sapply(sval, function(x)asn0(x, n, th0, M))
      gridlistnew[[i]] = list(sval = sval, rhoval = rhoval)
    }
    gridlist1 = gridlistnew
    n = n - 1
    if (n > 1) {
      repeat {
        nInt = test$data[[n - 1]]$nInt
        cInt = test$data[[n - 1]]$cInt
        gridlistnew = list()
        for (i in 1:nInt) {
          sval = makegrid(cInt[i, 1], cInt[i, 2], gridsz)
          rhoval = sapply(sval, function(x)asn1(x, n, th0, M))
          gridlistnew[[i]] = list(sval = sval, rhoval = rhoval)
        }
        gridlist1 = gridlistnew
        if (n == 2)break
        n = n - 1
      }
    }
  }
  return(ifelse(length(M) > 1, asn1(0, 1, th0, M), asn0(0, 1, th0, M)))
}

DBCPlan <- function(M, l, th, gam, thgam, gridsz, cost_fn, margin = 20) {

  exper1 <- function(x, n, l, th, M) {
    cost_fn(M[n], n) * weighted(x, gam, thgam, cumsum(M)[n - 1])
  }


  eff1 <- function(x, n, l, th, M)mincomp(x, n - 1, l, th, M) - exper1(x, n, l, th, M)
  eff1.v <- function(x, nd, l, th, M)Vectorize(eff1, vectorize.args = "x")(x = x, n = nd, l = l, th = th, M = M)

  ##############################
  n = length(M)
  k = length(th)
  test = list()

  a = dividing(n, l, th, M)

  test[[n]] = list(n = n, div = a)

  a = dividing(n - 1, l, th, M)
  cInt = array(dim = c(k - 1, 2))

  left = uniroot(eff1, c(a[1] - margin, a[1]), n = n, l = l, th = th, M = M, extendInt = "upX")$root
  right = uniroot(eff1, c(a[k - 1], a[k - 1] + margin), n = n, l = l, th = th, M = M, extendInt = "downX")$root

  limits = sort(uniroot.all(eff1.v, lower = left - margin, upper = right + margin, nd = n, l = l, th = th, M = M))
  nInt = length(limits) / 2

  current = 1
  cInt = array(dim = c(k - 1, 2))
  for (i in 1:nInt) {
    cInt[i, 1] = limits[current]
    cInt[i, 2] = limits[current + 1]
    current = current + 2
  }
  test[[n - 1]] = list(n = n - 1, div = a, nInt = nInt, cInt = cInt)

  n = n - 1
  if (n > 1)
  {
    repeat {
      a = dividing(n - 1, l, th, M)
      left = uniroot(eff1, c(a[1] - margin, a[1]), n = n, l = l, th = th, M = M, extendInt = "upX")$root
      right = uniroot(eff1, c(a[k - 1], a[k - 1] + margin), n = n, l = l, th = th, M = M, extendInt = "downX")$root
      limits = sort(uniroot.all(eff1.v, lower = left - margin, upper = right + margin, nd = n, l = l, th = th, M = M))
      nInt = length(limits) / 2

      current = 1
      cInt = array(dim = c(k - 1, 2))
      for (i in 1:nInt) {
        cInt[i, 1] = limits[current]
        cInt[i, 2] = limits[current + 1]
        current = current + 2
      }

      test[[n - 1]] = list(n = n - 1, div = a, nInt = nInt, cInt = cInt)
      if (n == 2)break
      n = n - 1
    }
  }

  return(list(data = test, info = list(theta = th, lambda = l, M = M, gridsz = gridsz, vartheta = thgam, gamma = gam, const_fn = cost_fn)))

}

##########################################################################
PAccept <- function(test, th0, fin, gridsz) {

  prob0 <- function(x, n, M)
  {
    a = test$data[[n]]$div
    k = length(test$info$th)
    #   tmp=array(0,dim=c(k,k))
    s = 0
    ext.div = c(-Inf, a, Inf)
    for (current in 1:k) {
      if (current == fin) {
        s = s + (sum.cdf(ext.div[current + 1], x, th0, M[n]) - sum.cdf(ext.div[current], x, th0, M[n])) *
          lr(x, th0, ifelse(n > 1, cumsum(M)[n - 1], 0))
      }
    }
    return(s)
  }

  prob1 <- function(x, n, M) {
    k = length(test$info$th)
    nInt = test$data[[n]]$nInt
    cInt = test$data[[n]]$cInt
    a = test$data[[n]]$div
    ext.div = c(-Inf, a, Inf)
    cut = ext.div
    s = 0
    from = 1
    to = max(which(ext.div < cInt[1, 1])) + 1
    cut[to] = cInt[1, 1]
    current = from
    repeat {
      if (current == fin) {
        s = s + (sum.cdf(cut[current + 1], x, th0, M[n]) - sum.cdf(cut[current], x, th0, M[n])) *
          lr(x, th0, ifelse(n > 1, cumsum(M)[n - 1], 0))
      }
      current = current + 1
      if (current == to)break
    }
    from = min(which(ext.div > cInt[nInt, 2])) - 1
    to = k + 1
    cut = ext.div
    cut[from] = cInt[nInt, 2]
    current = from
    repeat {
      if (current == fin) {
        s = s + (sum.cdf(cut[current + 1], x, th0, M[n]) - sum.cdf(cut[current], x, th0, M[n])) *
          lr(x, th0, ifelse(n > 1, cumsum(M)[n - 1], 0))
      }
      current = current + 1
      if (current == to)break
    }
    if (nInt > 1)
      for (int in 2:nInt) {
        cut = ext.div
        from = max(which(ext.div < cInt[int - 1, 2]))
        to = min(which(ext.div > cInt[int, 1]))
        cut[from] = cInt[int - 1, 2]
        cut[to] = cInt[int, 1]
        current = from
        repeat {
          if (current == fin) {
            s = s + (sum.cdf(cut[current + 1], x, th0, M[n]) - sum.cdf(cut[current], x, th0, M[n])) *
              lr(x, th0, ifelse(n > 1, cumsum(M)[n - 1], 0))
          }
          current = current + 1
          if (current == to)break
        }
      }
    for (int in 1:nInt)
      s = s + intgr(gridlist1[[int]]$sval, gridlist1[[int]]$rhoval, x, M[n])
    return(s)
  }

  n = length(test$info$M)
  k = length(test$info$th)
  M = test$info$M
  if (n > 1) {
    nInt = test$data[[n - 1]]$nInt
    cInt = test$data[[n - 1]]$cInt
    gridlistnew = list()
    for (i in 1:nInt) {
      sval = makegrid(cInt[i, 1], cInt[i, 2], gridsz)
      rhoval = sapply(sval, function(x)prob0(x, n, M))
      gridlistnew[[i]] = list(sval = sval, rhoval = rhoval)
    }
    gridlist1 = gridlistnew
    n = n - 1
    if (n > 1)
      repeat {
        nInt = test$data[[n - 1]]$nInt
        cInt = test$data[[n - 1]]$cInt
        gridlistnew = list()
        for (i in 1:nInt) {
          sval = makegrid(cInt[i, 1], cInt[i, 2], gridsz)
          rhoval = sapply(sval, function(x)prob1(x, n, M))
          gridlistnew[[i]] = list(sval = sval, rhoval = rhoval)
        }
        gridlist1 = gridlistnew
        if (n == 2)break
        n = n - 1
      }
  }
  if (length(M) > 1)
    return(prob1(0, 1, M))
  else   return(prob0(0, 1, M))
}

