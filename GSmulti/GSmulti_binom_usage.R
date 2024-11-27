source("./GSmulti/GSmulti_binom.R")

cost_fn <- function(x, k) x
print("#################################################", quote = F)
print("Example for binomial distribution", quote = F)
print("#################################################", quote = F)
M = rep(20, 5) # 5 groups of 20 observations each
prob = c(0.4, 0.6, 0.8) # binom prob parameter
print("Hypotheses about prob parameter", quote = F)
th = log(prob / (1 - prob)) # th natural parameter FYI
print(prob)
k = length(th)
gam = rep(1, k) / k
thgam = th
maxval = max(M)
mylam = c(70, 258, 105)
lam = array(c(NA, mylam[2], mylam[3], mylam[1], NA, mylam[3], mylam[1], mylam[2], NA), dim = c(length(th), length(th)))
# matrix of Lagrange multipliers
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
  print(paste("Operating characteristics"), quote = F)
  print(alpha[i,])
  print(paste("simulated"), quote = F)
  print(res$OC)
  print(paste("Expected number of groups (simulated)", res$groups), quote = F)
}

