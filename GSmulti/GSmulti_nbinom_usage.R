source("./GSmulti/GSmulti_nbinom.R")
cost_fn <- function(x)x
print("#################################################", quote = F)
print("Example for negative binomal distribution", quote = F)
print("#################################################", quote = F)
size.nbinom = 3 # overwrite the default value of size parameter
M = rep(20, 5) # 5 groups of 20 observations each
th = c(-2.5, -1.5, -1) #natural parametrization of the negative binomial distribution
print("Hypotheses about prob parameter", quote = F)
print(1 - exp(th) / 2) # dnbinom prob parameter FYI
print(paste("size parameter", size.nbinom), quote = F)
k = length(th)
gam = rep(1, k) / k
thgam = th
maxval = 50
lam = array(200, dim = c(length(th), length(th))) # matrix of Lagrange multipliers
for (i in 1:k) lam[i, i] = NA # not used
print("Lagrangian multipliers")
print(lam)
for (i in 1:k) lam[i, i] = 0
test = OptPlan(M, lam, th, gam, thgam, cost_fn, maxval)
alpha = array(0, dim = (c(length(th), length(th))))
for (l in 1:k) {
  for (m in 1:k) alpha[l, m] = PAccept(test, th[l], m)
}
print("Operating characteristics", quote = F)
print(alpha)
print("#################################################", quote = F)
print(paste("Calculated Lagrangian", sum(alpha * lam) + sum(gam * sapply(thgam, function(x)ESS(test, x, cost_fn)))), quote = F)
print(paste("internal Lagrangian", test$info$Lagr), quote = F)

for (i in 1:length(th)) {
  print("#################################################", quote = F)

  print(paste("Under hypothesis ", i), quote = F)
  print("#################################################", quote = F)

  res = monte_carlo_simulation(test, hyp = th[i], nMC = 100000)
  print(paste("ESS", ESS(test, th[i], cost_fn)), quote = F)
  print(paste("simulated", res$ESS), quote = F)
  # res=sapply(1:length(th),function(x)PAccept(test,th[i],x))
  # print(res)
  print(paste("Operating characteristics"), quote = F)
  print(alpha[i,])
  print(paste("simulated"), quote = F)
  print(res$OC)
  print(paste("Expected number of groups (simulated)", res$groups), quote = F)
}
