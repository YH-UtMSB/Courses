rm(list = ls())

set.seed(1022)
library(readr)
library(mosaic)
library(truncnorm)

polls = read_csv(".../path/to/polls.csv")
polls = polls[complete.cases(polls),]



# Use geometry average male/female ratio as the test statistic

# Test the independence of state and gender

df.s = table(polls$state, polls$female)
Nm.s = df.s[,1]; Nf.s = df.s[,2]
idx = intersect(which(Nm.s > 0), which(Nf.s > 0))

test.s = (prod(Nm.s[idx] / Nf.s[idx]))^(1 / length(idx))

# Do the permutation test
H0.s = do(1000)*{
  df.s = table(polls$state, shuffle(polls$female))
  Nm.s = df.s[,1]; Nf.s = df.s[,2]
  idx = intersect(which(Nm.s > 0), which(Nf.s > 0))
  if(length(idx) >= 1){
    (prod(Nm.s[idx] / Nf.s[idx]))^(1 / length(idx))
  }
}
H0.s = H0.s$result
p.s = length(c(which(H0.s < test.s), which(H0.s > (2*mean(H0.s)-test.s)))) / length(H0.s)


# Test the independence of age and gender

df.a = table(polls$age, polls$female)
Nm.a = df.a[,1]; Nf.a = df.a[,2]
idx = intersect(which(Nm.a > 0), which(Nf.a > 0))

test.a = (prod(Nm.a[idx] / Nf.a[idx]))^(1 / length(idx))

# Do the permutation test
H0.a = do(1000)*{
  df.a = table(polls$age, shuffle(polls$female))
  Nm.a = df.a[,1]; Nf.a = df.a[,2]
  idx = intersect(which(Nm.a > 0), which(Nf.a > 0))
  if(length(idx) >= 1){
    (prod(Nm.a[idx] / Nf.a[idx]))^(1 / length(idx))
  }
}
H0.a = H0.a$result
p.a = length(c(which(H0.a < test.a), which(H0.a > (2*mean(H0.a)-test.a)))) / length(H0.a)


# Test the independence of edu and gender

df.e = table(polls$edu, polls$female)
Nm.e = df.e[,1]; Nf.e = df.e[,2]
idx = intersect(which(Nm.e > 0), which(Nf.e > 0))

test.e = (prod(Nm.e[idx] / Nf.e[idx]))^(1 / length(idx))

# Do the permutation test
H0.e = do(1000)*{
  df.e = table(polls$edu, shuffle(polls$female))
  Nm.e = df.e[,1]; Nf.e = df.e[,2]
  idx = intersect(which(Nm.e > 0), which(Nf.e > 0))
  if(length(idx) >= 1){
    (prod(Nm.e[idx] / Nf.e[idx]))^(1 / length(idx))
  }
}
H0.e = H0.e$result
p.e = length(c(which(H0.e < test.e), which(H0.e > (2*mean(H0.e)-test.e)))) / length(H0.e)



# Test the independence of ethnicity and gender

df.k = table(polls$black, polls$female)
Nm.k = df.k[,1]; Nf.k = df.k[,2]
idx = intersect(which(Nm.k > 0), which(Nf.k > 0))

test.k = (prod(Nm.k[idx] / Nf.k[idx]))^(1 / length(idx))

# Do the permutation test
H0.k = do(1000)*{
  df.k = table(polls$black, shuffle(polls$female))
  Nm.k = df.k[,1]; Nf.k = df.k[,2]
  idx = intersect(which(Nm.k > 0), which(Nf.k > 0))
  if(length(idx) >= 1){
    (prod(Nm.k[idx] / Nf.k[idx]))^(1 / length(idx))
  }
}
H0.k = H0.k$result
p.k = length(c(which(H0.k < test.k), which(H0.k > (2*mean(H0.k)-test.k)))) / length(H0.k)



par(mfrow = c(2,2))

hist(
  H0.s, breaks = 20, 
  main = "State",
  xlab = "geoavg of m/f ratio",
  probability = TRUE
)
abline(v = test.s, col = "red")

hist(
  H0.a, breaks = 20, 
  main = "Age",
  xlab = "geoavg of m/f ratio",
  probability = TRUE
)
abline(v = test.a, col = "red")

hist(
  H0.e, breaks = 20, 
  main = "Edu level",
  xlab = "geoavg of m/f ratio",
  probability = TRUE
)
abline(v = test.e, col = "red")

hist(
  H0.k, breaks = 20, 
  main = "Ethnicity",
  xlab = "geoavg of m/f ratio",
  probability = TRUE
)
abline(v = test.k, col = "red")

par(mfrow = c(1,1))


cat(
  "empirical p-values:\n",
  "p-state: ", p.s, "\n",
  "p-age: ", p.a, "\n", 
  "p-edu: ", p.e, "\n",
  "p-black: ", p.k, "\n"
)