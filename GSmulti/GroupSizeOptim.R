#Bayes
source("./GSmulti/GSmulti_norm.R")

# 2 hypotheses group sequential test with sample size standard deviation
# alternate between Lagrange minimization and alpha and beta fitting
# within the loop over number of groups and group overhead cost
delta=0.1
th = c(-delta, delta)
alpha=c(0.05,0.1)

thgam = th
gridsz =0.02
#FSS=324.970443165642# ESS for 1 group
print("Hypotheses about normal means (unitary variance)", quote = F)
print(th)
k = length(th)
gam = rep(0.5,2)
l=c( 690,690)
cost_fn <- function(x,k)pergroup+x
print("#################################################", quote = F)
print("Example for normal distribution", quote = F)
print("#################################################", quote = F)


eq3 <- function(x,th, gam, thgam, cost_fn, gridsz, alpha) {
  #FSS finding
  print(x)
   lam=array(dim=c(2,2))
   lam[2,1]=x[2]
   lam[1,2]=x[1]
  M=x[3]
 test=   OptPlan(M, lam, th, gam, thgam, gridsz, cost_fn)
  al = sapply(1:2,function(i)1-PAccept(test, th[i], i, gridsz))
  tmp = max(abs(al - alpha) / alpha)
print(al)
   print(tmp)
  return(tmp)
}

FSS_fn<-function(alpha,th=th,gam=gam,gridsz=gridsz,thgam=thgam,cost_fn=cost_fn){
x=c(300,300,100)
res1=optim(x,eq3,alpha=alpha,th=th,gam=gam,gridsz=gridsz,thgam=thgam,cost_fn=cost_fn,method="Nelder-Mead",control=list(abstol=0.00001))
x=res1$par
  return(x)
}

 #FSS=324.870443165642
monte_carlo_simulation <- function(test, hyp,
                                   nMC = 10000)
  # Monte Carlo simulation of a test
{
  M = test$info$M
  H = length(test$info$M)
  K = length(test$info$th)
  ss = 0
  totmss=0
  totaccepted = rep(0, K)
  totn = 0
  totobs = 0
  # runobs=0
  for (i in 1:nMC) {
    s = 0 #accumulated sum for test run
    n = 0 # number of steps in current run
    runobs = 0
    mss=0
    accepted = rep(0, K) #number of accepting in the run (0 or 1)
    for (stage in 1:(H)) { #run starts
      #generate
      summand = rn(hyp, M[stage])
      s = s + summand #accumulated sum
      n = n + 1  #one step more
      runobs = runobs + M[stage]
      mss=mss+M[stage]^2
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
    totmss=totmss+runobs^2
    totn = totn + n
    totobs = totobs + runobs
    ss = ss + n^2
  }
  OC = (totaccepted) / as.double(nMC)
  totgroups = totn / as.double(nMC)
  obs = totobs / as.double(nMC)
  mvar=totmss/as.double(nMC)-obs^2
  return(list(OC = OC, groups = totgroups, ESS = obs,mvar=mvar))
}

monte_carlo_simulation <- function(test, hyp,cost_fn=function(x,k)x,
                                   nMC = 10000)
  # Monte Carlo simulation of a test
{
  M = test$info$M
  H = length(test$info$M)
  K = length(test$info$th)
  ss = 0
  totmss=0
  totaccepted = rep(0, K)
  totn = 0
  totobs = 0
  totcost=0
  # runobs=0
  for (i in 1:nMC) {
    s = 0 #accumulated sum for test run
    n = 0 # number of steps in current run
    runobs = 0
    cost=0
    mss=0
    accepted = rep(0, K) #number of accepting in the run (0 or 1)
    for (stage in 1:(H)) { #run starts
      #generate
      summand = rn(hyp, M[stage])
      s = s + summand #accumulated sum
      n = n + 1  #one step more
      runobs = runobs + M[stage]
      cost=cost+cost_fn(M[stage],k)
      mss=mss+M[stage]^2
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
    totmss=totmss+runobs^2
       totn + n
    totobs = totobs + runobs
    totcost=totcost+cost
    ss = ss + n^2
  }
  OC = (totaccepted) / as.double(nMC)
  totgroups = totn / as.double(nMC)
  obs = totobs / as.double(nMC)
  mvar=totmss/as.double(nMC)-obs^2
  return(list(OC = OC, groups = totgroups, ESS = obs,sd=sqrt(mvar),cost=totcost/as.double(nMC)))
}

sec_mom_SS<- function(test, th0, gridsz) {

  asn0 <- function(x, n, th0, M) {
    return(((cumsum(M)[n])^2-(cumsum(M)[n-1])^2) * lr(x, th0, ifelse(n > 1, cumsum(M)[n - 1], 0)))
  }

  asn1 <- function(x, n, th0, M) {
    k = length(test$info$th)
    nInt = test$data[[n]]$nInt
    cInt = test$data[[n]]$cInt
    tmp1 = ((cumsum(M)[n])^2-ifelse(n>1,(cumsum(M)[n-1])^2,0)) * lr(x, th0, ifelse(n > 1, cumsum(M)[n - 1], 0))
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
minlagr <- function(x, lam,th, gam, thgam, cost_fn, gridsz) {
 # if(any(x<10))return(1000)
  test= OptPlan(x, lam, th, gam, thgam, gridsz, cost_fn)
tmp=test$info$Lagr
 # print(tmp)
  return(tmp)

}

eq <- function(x,M, th, gam, thgam, cost_fn, gridsz, alpha) {
  #FSS finding
 # print(x)


  if(any(x<=cost_fn(M[1],1)))return(1000)

   lam=array(dim=c(2,2))
   lam[2,1]=x[2]
   lam[1,2]=x[1]
 #res=optim(M,minlagr,lam=lam,th=th,gam=gam,gridsz=gridsz,thgam=thgam,cost_fn=cost_fn,method="Nelder-Mead",control=list(abstol=0.000001))
 #x=res$par
 test=   OptPlan(M, lam, th, gam, thgam, gridsz, cost_fn)
  al = sapply(1:2,function(i)1-PAccept(test, th[i], i, gridsz))
  tmp = max(abs(al - alpha) / alpha)
print(al)
 #  print(tmp)
  return(tmp)
}

eq1 <- function(x,M,th, gam, thgam, cost_fn, gridsz, alpha) {
  #FSS finding
   print(x)
   lam=array(dim=c(2,2))
   lam[2,1]=x[2]
   lam[1,2]=x[1]
 test=   OptPlan(M, lam, th, gam, thgam, gridsz, cost_fn)
  al = sapply(1:2,function(i)1-PAccept(test, th[i], i, gridsz))

  tmp = max(abs(al - alpha) / alpha)
 print(al)
   print(tmp)
  return(tmp)
}
FSS_all=FSS_fn(alpha=alpha,th=th,gam=gam,gridsz=gridsz,thgam=thgam,cost_fn=function(x,k)x)
#FSS=(qnorm(alpha[1])+qnorm(alpha[2]))^ 2/(2*delta)^ 2
print(paste("FSS",FSS_all[3]))
l=FSS_all[1:2]
     for (ngr in c(4)){
      for(pergroup in c(0,1,2,3,5,10,20))

       {

print(date())
         print(ngr)
         print(pergroup)
M=rep(FSS_all[3]*2/ngr,ngr)
         lam=array(dim=c(2,2))
  lam[2,1]=l[2]
  lam[1,2]=l[1]

repeat{
  #l=c(800,800)
#l=c(435.6556, 831.1367)
 #M=rep(FSS*3/ngr,ngr)
  lam=array(dim=c(2,2))
  lam[2,1]=l[2]
  lam[1,2]=l[1]

  res=optim(l,eq,M =M,alpha=alpha,th=th,gam=gam,gridsz=gridsz,thgam=thgam,cost_fn=cost_fn,method="Nelder-Mead",control=list(abstol=0.001))
  l=res$par
     print(l)
lam=array(dim=c(2,2))
  lam[2,1]=l[2]
  lam[1,2]=l[1]

 res=optim(M,minlagr,lam=lam,th=th,gam=gam,gridsz=gridsz,thgam=thgam,cost_fn=cost_fn,method="Nelder-Mead",control=list(abstol=0.0001))
 x=res$par
  old=res$value

  print(x)
  print(old)
  M=x

 test=   OptPlan(M, lam, th, gam, thgam, gridsz, cost_fn)
  al = sapply(1:2,function(i)1-PAccept(test, th[i], i, gridsz))
       print(al)
  if(max(abs(al-alpha)/alpha)<0.001)break
 # if(abs(old-test$info$Lagr)<0.01)break
}
           test = OptPlan(M, lam, th, gam, thgam, gridsz, cost_fn)

  eff=sum(gam*sapply(length(thgam),function(i)ESS(test, thgam[i],gridsz,cost_fn)))/cost_fn(FSS_all[3],1)
print(paste("eff",eff))
stdev=sum(gam*sapply(length(thgam),function(i) (sqrt(sec_mom_SS(test, thgam[i],gridsz)-ESS(test, thgam[i],gridsz,cost_fn=function(x, k=1) x^2)))))/FSS_all[3]
print(stdev)
#t=ngr* test$info$M[1]/FSS
s1 = c("Bayes",pergroup,al, ngr,  lam, test$info$M, FSS_all[3], eff,stdev)
    write.table(rbind(s1), "./GSmulti.csv", append = TRUE, sep = ",", col.names = FALSE)
print(date())

x=l
  lam=array(dim=c(2,2))
  lam[2,1]=x[2]
  lam[1,2]=x[1]

M1=rep(sum(M)/ngr,ngr)
res1=optim(l,eq1,M=M1,alpha=al,th=th,gam=gam,gridsz=gridsz,thgam=thgam,cost_fn=cost_fn,method="Nelder-Mead",control=list(abstol=0.001))
l=res1$par
lam[1,2]=l[1]
lam[2,1]=l[2]
test1 = OptPlan(M1, lam, th, gam, thgam, gridsz, cost_fn)
al=array(dim=2)
           al[1]=   (1- PAccept(test1, th[1], 1, gridsz))
           al[2]=    (1-PAccept(test1, th[2], 2, gridsz))
# sapply(1:3,function(i)ESS(test, th[i], gridsz, cost_fn))
dev=sum(gam*sapply(length(thgam),function(i) (sqrt(sec_mom_SS(test1, thgam[i],gridsz)-ESS(test1, thgam[i],gridsz,function(x,k=1)x^2))))/FSS_all[3])

  eff1=sum(gam*sapply(length(thgam),function(i)ESS(test1, thgam[i],gridsz,cost_fn)))/cost_fn(FSS_all[3],k=1)
  ESSv=sapply(length(thgam),function(i)ESS(test, thgam[i],gridsz,function(x,k=1)x))
print(paste("eff1",eff1))
stdev=sum(gam*sapply(length(thgam),function(i) (sqrt(sec_mom_SS(test1, thgam[i],gridsz)-ESS(test1, thgam[i],gridsz,function(x,k=1)x^2)))))/FSS_all[3]
print(stdev)
print(eff1/eff)
s1 = c( pergroup,al,ngr,M1,  lam, test$info$M,  eff1,stdev,ESSv)
    write.table(rbind(s1), "./GSmulti.csv", append = TRUE, sep = ",", col.names = FALSE)
 write.table(eff1/eff, "./GSmulti.csv", append = TRUE, sep = ",", col.names = FALSE)
       }
     }
  