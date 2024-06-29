source("./GSmulti/GSmulti_pois.R")
cost_fn <- function(x)x
print("#################################################", quote = F)
print("Example for Poisson distribution", quote = F)
print("#################################################", quote = F)
M = rep(40, 10) # 10 groups of 40 observations each
th = c(-0.2, 0, 0.2)
print("Hypotheses about mean", quote = F)
print(exp(th)) # dnbinom prob parameter FYI
k = length(th)
gam = rep(1, k) / k
thgam = th
maxval = 90 # cropping level
lam = array(180, dim = c(length(th), length(th))) # matrix of Lagrange multipliers
lam[2, 1] = 350
lam[2, 3] = 350
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
