# This is a collaborative work by
# Shuying Wang and Yilin He.


rm(list = ls())
set.seed(1022)

library(readr)
library(mosaic)
library(truncnorm)

polls = read_csv(".../path/to/polls.csv")
polls = polls[complete.cases(polls),]

# set factors
polls$state = factor(polls$state); polls$age = factor(polls$age); polls$edu = factor(polls$edu)

polls = polls[order(polls$state), ]


all_states = levels(polls$state); Mst = as.numeric(summary(polls$state))
all_edu = levels(polls$edu); Med = as.numeric(summary(polls$edu))
all_age = levels(polls$age); Mag = as.numeric(summary(polls$age))

x6 = (polls$weight - mean(polls$weight)) / sd(polls$weight)

y = polls$bush
N = length(y)


# initialize model parameters
mu.st = rep(0, length(Mst))
mu.ed = rep(0, length(Med))
mu.ag = rep(0, length(Mag))
beta = matrix(rep(0, 7), nrow = 1, ncol = 7)

x0 = rep(1, N)
x1 = x2 = x3 = x4 = x5 = rep(0, N)
X = cbind(x0, x1, x2, x3, x4, x5, x6)
lbd = 1; phi = 0

# start running the chain
Ite = 5000; sta = 1000; st = Ite - sta + 1; ed = Ite

Acc = c()
Phi = c()
Beta = matrix(nrow = sta, ncol = 7)
Mu.st = matrix(nrow = sta, ncol = 49)
Mu.ag = Mu.ed = matrix(nrow = sta, ncol = 4)
X1 = X2 = X3 = X4 = X5 = matrix(nrow = sta, ncol = N)

for (iter in 1:Ite){
  z = ifelse(y == 1, rtruncnorm(1, 0, Inf, X %*% t(beta), 1), rtruncnorm(1, -Inf, 0, X %*% t(beta), 1))
  
  # Sample beta from N(inv(A)B, inv(A))
  # A = X'X + lbd * I
  # B = sum (zi * Xi) + lbd * phi
  A.inv = solve(t(X) %*% X + lbd * diag(7))
  B = apply(z * X, 2, sum) + lbd * phi
  beta = rmvnorm(1, A.inv %*% B, A.inv)
  
  
  # Set up Mu matrix
  Mu = cbind(mu.st[polls$state], mu.ed[polls$edu], mu.ag[polls$age], (2*polls$female-1), (2*polls$black-1))
  
  # Sample xi from N(inv(Ax)Bxi, inv(Ax))
  # Ax = beta * beta' + I
  # Bxi = beta * (zi - b0 - b6 * x6i) + mui
  b = beta[2:6]
  Ax.inv = solve(b %*% t(b) + diag(5))
  Bx = t(Mu + (z - beta[1] - beta[7] * x6) %*% t(b))
  X = apply(Ax.inv %*% Bx, 2, function(mean) rmvnorm(1, mean, Ax.inv))
  
  x.st = X[1,]; x.ed = X[2,]; x.ag = X[3,]; x.fm = X[4,]; x.bk = X[5,]
  X = cbind(x0, t(X), x6)
  
  # Sample mu.st, mu.ed, mu.ag
  mean.st = sapply(all_states, function(s) mean(x.st[which(polls$state == s)]))
  var.st = 1 / Mst
  mu.st = rnorm(length(Mst), mean.st*var.st, var.st^(0.5))
  
  mean.ed = sapply(all_edu, function(e) mean(x.ed[which(polls$edu == e)]))
  var.ed = 1 / Med
  mu.ed = rnorm(length(Med), mean.ed*var.ed, var.ed^(0.5))
  
  mean.ag = sapply(all_age, function(a) mean(x.ag[which(polls$age == a)]))
  var.ag = 1 / Mag
  mu.ag = rnorm(length(Mag), mean.ag*var.ag, var.ag^(0.5))
  
  phi = rnorm(1, mean(beta), (7 * lbd)^(-0.5))
  lbd = rgamma(1, 4, 0.5 + 0.5 * sum((beta - phi)^2))
  
  if(iter >= st){
    y.hat = c(X %*% t(beta)) > 0
    acc = 1 - sum(abs(y - y.hat)) / N
    Acc = c(Acc, acc)
    print(acc)
    
    Mu.st[iter - st + 1, ] = mu.st
    Mu.ed[iter - st + 1, ] = mu.ed
    Mu.ag[iter - st + 1, ] = mu.ag
    
    Phi = c(Phi, phi)
    Beta[iter - st + 1, ] = beta 
    
    X1[iter - st + 1, ] = x.st
    X2[iter - st + 1, ] = x.ed
    X3[iter - st + 1, ] = x.ag
    X4[iter - st + 1, ] = x.fm
    X5[iter - st + 1, ] = x.bk
  }
  print(iter)
}

Result = cbind(Acc, Phi, Mu.ed, Mu.ag, Beta)
write.csv(Mu.st, file = "Mu.st.csv")
write.csv(Result, file = "results.csv")
write.csv(X1, file = "state_score.csv")
write.csv(X2, file = "edu_score.csv")
write.csv(X3, file = "age_score.csv")
write.csv(X4, file = "female_score.csv")
write.csv(X5, file = "black_score.csv")

