source("./GSmulti/GSmulti_exp.R")
cost_fn <- function(x)x
print("#################################################", quote = F)
print("Example for exponential distribution", quote = F)
print("#################################################", quote = F)
M = rep(50, 5) # 5 groups of 50 observations each
th = c(-0.3, 0, 0.3)
#  vector of canonical parameter
print("Hypotheses about the rate parameter", quote = F)
rate = 1 - th # vector of rate parameters
print(rate)
k = length(th)
gam = rep(1, k) / k
thgam = th
gridsz = 0.05
lam = array(590, dim = c(length(th), length(th))) # matrix of Lagrange multipliers
for (i in 1:k) lam[i, i] = NA # not used
print("Lagrangian multipliers")
print(lam)
for (i in 1:k) lam[i, i] = 0
test = OptPlan(M, lam, th, gam, thgam, gridsz, cost_fn)
alpha = array(0, dim = (c(length(th), length(th))))
for (l in 1:k) {
  for (m in 1:k) alpha[l, m] = PAccept(test, th[l], m, gridsz)
}
print("Operating characteristics", quote = F)
print(alpha)
print("#################################################", quote = F)
print(paste("Calculated Lagrangian", sum(alpha * lam) + sum(gam * sapply(thgam, function(x)ESS(test, x, gridsz, cost_fn)))), quote = F)
print(paste("internal Lagrangian", test$info$Lagr), quote = F)
for (i in 1:length(th)) {
  print("#################################################", quote = F)
  print(paste("Under hypothesis ", i), quote = F)
  print("#################################################", quote = F)
  res = monte_carlo_simulation(test, hyp = th[i], nMC = 100000)
  print(paste("ESS", ESS(test, th[i], gridsz, cost_fn)), quote = F)
  print(paste("simulated", res$ESS), quote = F)
  print(paste("Operating characteristics"), quote = F)
  print(alpha[i,])
  print(paste("simulated"), quote = F)
  print(res$OC)
  print(paste("Expected number of groups (simulated)", res$groups), quote = F)
}
