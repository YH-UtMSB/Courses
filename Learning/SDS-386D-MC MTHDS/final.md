## HW 6 & Final

### The Bayes Factor MTHD

```r
rm(list = ls())
library(truncnorm)
set.seed(12345)
# Prepare dataset.
datapath = "E:/R/final.dat"

raw_dat = read.table(datapath)
N = dim(raw_dat)[1]

y = as.numeric(raw_dat$V1)
X = matrix(c(rep(1,N), raw_dat$V2, raw_dat$V3), ncol = 3)
x1 = as.numeric(raw_dat$V2)
x2 = as.numeric(raw_dat$V3)

z = exp(y) - 1
# set negative elements in z as 1e-5
z = (z >= 0) * z + (z < 0) * 1e-5

# Computes joint likelihood k(data, theta| model)
jt.lkhd = function(theta){
  
  l = length(theta) - 2
  b = theta[-1]; lbd = theta[1]
  
  Xp = X[,1:(l+1)]
  
  mu = log(1 + exp(Xp %*% b))
  
  jt_lkhd = prod(dnorm(y, mu, lbd^(-.5))) * prod(dnorm(b, 0, 5)) * dgamma(lbd, .5, .5)
  
  return(jt_lkhd)
  
}



Ap = function(theta.f, theta.t){
  # Computing acceptance probability alpha
  # for jumping from theta.f to theta.t
  # theta = [lbd, b0, b1, (b2)]
  l = length(theta.f) - 2
  
  lbd.f = theta.f[1]; lbd.t = theta.t[1]
  
  p.n = dtruncnorm(lbd.f, 0, Inf, lbd.t) * jt.lkhd(theta.t)
  p.d = dtruncnorm(lbd.t, 0, Inf, lbd.f) * jt.lkhd(theta.f)
  
  r = p.n / p.d

  
  alp = min(1, r)
  
  return(alp)
}


param_init = function(m){
  # compute theta_init with least square estimate
  
  # to do this, must set z to be global variable
  mod = paste0("log(z) ~ 1 + ")
  vars = paste0("x", 1:m)
  mod_expr = paste0(mod, paste0(vars, collapse = " + "))
  
  fit = lm(mod_expr)
  beta = coef(fit)
  
  z.hat = predict(fit)
  y.hat = log(1 + exp(z.hat))
  
  sig = sd((y - y.hat))
  lbd = sig^(-2)
  
  theta_init = c(lbd, beta)
  varname = c("lbd", paste0("b", 0:m))
  names(theta_init) = varname
  
  return(theta_init)
}


theta.prop = function(theta.f){
  
  # propose new parameters
  lbd.n = rtruncnorm(1, a = 0, b = Inf, theta.f[1], 1)
  betas.n = rnorm((m+1), theta[-1], 1)
  
  theta.t = c(lbd.n, betas.n)
  names(theta.t) = varname
  
  return(theta.t)
  
}





# Running gibbs sampling
sta = 10000
Ite = 100000; st = Ite - sta + 1; ed = Ite

result = matrix(rep(0,4), nrow = 2)
rownames(result) = c("M = 1", "M = 2")
colnames(result) = c("margin_lkhd", "Pr(M|data)")

params = list("Theta1" = c(), "Theta2" = c())

for (m in c(1,2)) {

  comp = paste0("Theta", m)
  Theta = c()
  
  theta = param_init(m)
  varname = names(theta)
  
  for (iter in 1:Ite) {
    
    # propose new parameters
    theta.f = theta
    theta.t = theta.prop(theta.f)
    
    
    # Run MH on thetas
    ap = Ap(theta.f, theta.t)
    u = runif(1)
    theta = theta.t * (u <= ap) + theta.f * (u > ap)
    
    Theta = cbind(Theta, theta)
    
  }
  
  params[[comp]] = Theta[,st:ed]
  
  # Plot traceplots
  par(mfrow = c(2,2))
  for (idx in 1:(m+2)) {
    plot(Theta[idx,], type = "l", main = varname[idx])
  }
  par(mfrow = c(1,1))
  
  
  # Obtain theta_tilde
  theta.tl = apply(Theta, 1, mean)
  names(theta.tl) = varname
  
  
  # Estimate Pi_hat
  alp.pi = apply(Theta, 2, function(theta.f) Ap(theta.f, theta.tl))
  q.pi = dtruncnorm(theta.tl["lbd"], 0, Inf, Theta[1,]) * 
    apply(Theta[-1,], 2, function(betas) prod(dnorm(theta.tl[-1], betas, 1)))
  pi.hat = mean(alp.pi * q.pi)
  
  
  # Estimate alpha_hat
  Alp = c()
  for (j in 1:sta) {
    
    theta.pr = theta.prop(theta.tl)
    alp.i = Ap(theta.tl, theta.pr)
    
    Alp = c(Alp, alp.i)
    
  }
  alp.hat = mean(Alp)
  
  
  # Estimate f(theta| data, model)
  f.hat = pi.hat / alp.hat
  mg_lkhd = jt.lkhd(theta.tl) / f.hat
  
  P_model = mg_lkhd / 2
  
  result[m, ] = c(mg_lkhd, P_model)
  
}

result[, "Pr(M|data)"] = result[, "Pr(M|data)"] / sum(result[, "Pr(M|data)"])
result
BF = result["M = 1", "margin_lkhd"] / result["M = 2", "margin_lkhd"]
BF

```
