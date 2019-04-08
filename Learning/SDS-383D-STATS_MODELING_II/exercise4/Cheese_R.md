#### Half finished ...

```r
library(lme4)
library(tidyverse)
library(lattice)
library(readr)

cheese = read_csv("E:/R/cheese.csv")

cheese$store = factor(cheese$store)
cheese = cheese[order(cheese$store),]

# collect the number of data point for each store
N = as.numeric(summary(cheese$store))
n = length(N)


edpt = cumsum(N)
stpt = c(1,(edpt+1)[-n])

m = 50
tr_idx = c(matrix(rep(stpt, times = m), nrow = m, byrow = TRUE) + (0:(m-1)))


# split dataset into y (training set) and y.va (validation set)
Q = cheese$vol[tr_idx]; Q.va = cheese$vol[-tr_idx]
P = cheese$price[tr_idx]; P.va = cheese$price[-tr_idx]
D = cheese$disp[tr_idx]; D.va = cheese$disp[-tr_idx]

y = log(matrix(Q, nrow = n, byrow = TRUE))
x = log(matrix(P, nrow = n, byrow = TRUE))
D = matrix(D, nrow = n, byrow = TRUE)

N1 = sum(D)


# initialize parameters
# sigerr := sigma
alp.0 = 10; alp.1 = 0; eta.2 = eta.3 = 0; mu = -2; sigerr = fg = 1
phi.2 = rep(0, n); phi.3 = rep(0, n)

# set up parameter collectors
Alp.0 = Alp.1 = Eta.2 = Eta.3 = Mu = Sig = Fg = c()
Phi.2 = Phi.3 = c()

# running Gibbs sampling

Ite = 10000

for (t in 1:Ite) {
  
  # frequently used variables
  sigsq = sigerr^2
  
  # compute depdendent factors: beta.2 and beta.3
  beta.2 = eta.2 * sigerr * phi.2 + mu
  beta.3 = eta.3 * sigerr * phi.3
  
  
  # sample alp.0 from p(alp.0| ... ) then record
  a = y - alp.1 * D - beta.2 * x - beta.3 * D * x
  mean.alp0 = mean(a)
  sd.alp0 = (sigsq / (m*n))^(0.5)
  alp.0 = rnorm(1, mean.alp0, sd.alp0)
  Alp.0 = c(Alp.0, alp.0)
  
  
  # sample alp.1 from p(alp.1| ... ) then record
  a = y - alp.0 - beta.2 * x - beta.3 * D * x
  mean.alp1 = mean(a * D)
  sd.alp1 = (sigsq / N1)^(0.5)
  alp.1 = rnorm(1, mean.alp1, sd.alp1)
  Alp.1 = c(Alp.1, alp.1)
  
  
  # sample mu from p(mu| ... ) then record
  mij = (y - alp.0 - alp.1 * D - (beta.2 - mu) * x - beta.3 * D * x) / x
  mean.mu = sum(mij * x^2) / sum(x^2)
  sd.mu = (sigerr^2 / sum(x^2))^(0.5)
  mu = rnorm(1, mean.mu, sd.mu)
  Mu = c(Mu, mu)
  
  
  # sample eta.2 from p(eta.2| ... ) then record
  mu.p = mu / sigerr
  H2 = (y - alp.0 - alp.1 * D - beta.3 * D * x) / sigerr
  prec.eta2 = sum((phi.2 * x)^2) + 1
  mean.eta2.denom = sum((H2 - mu.p * x) * (phi.2 * x))
  mean.eta2 = mean.eta2.denom / prec.eta2
  eta.2 = rnorm(1, mean.eta2, prec.eta2^(-0.5))
  Eta.2 = c(Eta.2, eta.2)
  
  
  # update beta.2
  beta.2 = eta.2 * sigerr * phi.2 + mu
  
  
  # sample eta.3 from p(eta.3| ... ) then record
  H3 = (y - alp.0 - alp.1 - beta.2 * x) / sigerr
  prec.eta3 = sum((phi.3 * x)^2) + 1
  mean.eta3.denom = sum(H3 * phi.3 * x * D)
  mean.eta3 = mean.eta3.denom / prec.eta3
  eta.3 = rnorm(1, mean.eta3, prec.eta3^(-0.5))
  Eta.3 = c(Eta.3, eta.3)
  
  
  # update beta.3
  beta.3 = eta.3 * sigerr * phi.3
  
  
  # sample phi.2 from p(phi2| ... ) then record
  prec.phi2 = apply((eta.2 * x)^2, 2, sum) + fg
  Mphi = eta.2 * x * (y - alp.0 - alp.1 * D - beta.3 * D * x - mu * x) / sigerr
  mean.phi2 = apply(Mphi, 2, sum)
  Phi.2 = rnorm(n, mean.phi2, prec.phi2^(-0.5))
  Phi.2 = cbind(Phi.2, phi.2)
  
  
  # sample phi.3 from p(phi3| ... ) then record
  prec.phi3 = apply(D * (eta.3 * x)^2, 2, sum) + fg
  Mphi = D * eta.3 * x * (y - alp.0 - alp.1 - beta.2 * x) / sigerr
  mean.phi3 = apply(Mphi, 2, sum)
  phi.3 = rnorm(n, mean.phi3, prec.phi3^(-0.5))
  Phi.3 = cbind(Phi.3, phi.3)
  
  
  # sample fg from p(fg | ... ) then record, fg := g^(-2)
  fg.sp = n + .5
  fg.rt = (sum(phi.2^2 + phi.3^2) + 1) / 2
  fg = rgamma(1, fg.sp, fg.rt)
  g = fg^(-.5)
  Fg = c(Fg, fg)
  
  
  # sample 1/sigsq from p(1/sigsq| ... ) then record
  tau.2 = abs(eta.2) * g; tau.3 = abs(eta.3) * g
  fsgmsq.sp = m*n / 2 + n + 2
  fsgmsq.rt = 0.5 * (
    sum((y - alp.0 - alp.1 * D - beta.2 * x - beta.3 * D * x)^2) + 
      sum((beta.2 - mu)^2) / (tau.2^2) +
      sum(beta.3^2) / (tau.3^2)
  )
  fsgmsq = rgamma(1, fsgmsq.sp, fsgmsq.rt)
  sigerr = fsgmsq^(-0.5)
  Sig = c(Sig, sigerr)
  
}


# cross validation ... to be continued


```
