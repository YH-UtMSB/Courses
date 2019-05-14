## Load in data and results

X.F = read.csv(".../path/to/female_score.csv")
X.F = X.F[ , -1]
X.B = read.csv(".../path/to/black_score.csv")
X.B = X.B[ , -1]
X.S = read.csv(".../path/to/state_score.csv")
X.S = X.S[ , -1]
X.E = read.csv(".../path/to/edu_score.csv")
X.E = X.E[ , -1]
X.A = read.csv(".../path/to/age_score.csv")
X.A = X.A[ , -1]
Results = read.csv(".../path/to/Results.csv")
Beta = Results[1:1000, 12:18]
Mu.s = read.csv(".../path/to/mu.st.csv")

polls = read.csv(".../path/to/polls.csv")

polls = polls[complete.cases(polls),]
polls$state = factor(polls$state); polls$age = factor(polls$age); polls$edu = factor(polls$edu)
polls = polls[order(polls$state), ]


## Centralize the scores and beta
SD.F = apply(X.F, 1, sd)
Mean.F = apply(X.F, 1, mean)
SD.B = apply(X.B, 1, sd)
Mean.B = apply(X.B, 1, mean)
SD.S = apply(X.S, 1, sd)
Mean.S = apply(X.S, 1, mean)
SD.E = apply(X.E, 1, sd)
Mean.E = apply(X.E, 1, mean)
SD.A = apply(X.A, 1, sd)
Mean.A = apply(X.A, 1, mean)

X.F = t(apply(X.F, 1, function(x) (x - mean(x)) / sd(x)))
X.B = t(apply(X.B, 1, function(x) (x - mean(x)) / sd(x)))
X.S = t(apply(X.S, 1, function(x) (x - mean(x)) / sd(x)))
X.E = t(apply(X.E, 1, function(x) (x - mean(x)) / sd(x)))
X.A = t(apply(X.A, 1, function(x) (x - mean(x)) / sd(x)))

Beta[ , 1] = Beta[ , 1] + Beta[ , 2] * Mean.S + Beta[ , 3] * Mean.E+ Beta[ , 4] * Mean.A+ Beta[ , 5] * Mean.F + Beta[ , 6] * Mean.B
Beta[ , 2] = Beta[ , 2] * SD.S
Beta[ , 3] = Beta[ , 3] * SD.E
Beta[ , 4] = Beta[ , 4] * SD.A
Beta[ , 5] = Beta[ , 5] * SD.F
Beta[ , 6] = Beta[ , 6] * SD.B


## Take the beta and x scores of the last (1000th) iteration 

x.f = X.F[1000, ]
x.b = X.B[1000, ]
x.s = X.S[1000, ]
x.e = X.E[1000, ]
x.a = X.A[1000, ]
x.w = (polls$weight - mean(polls$weight)) / sd(polls$weight)
X = cbind(rep(1, 2015), x.s, x.e, x.a, x.f, x.b, x.w)
beta = unlist(Beta[1000, ])


## plot the gender and ethnicity scores
plot(x.f, x.b, type = 'n', xlab = 'Gendar score', ylab = 'ethnicity score')

female = polls$female
black = polls$black
b_m.ind = which(female == 0 & black == 1)
b_f.ind = which(female == 1 & black == 1)
nb_m.ind = which(female == 0 & black == 0)
nb_f.ind = which(female == 1 & black == 0)

text(x = x.f[b_m.ind], y = x.b[b_m.ind], labels = '\\MA', vfont = c('sans serif', 'bold'), col = 'blue')
text(x = x.f[b_f.ind], y = x.b[b_f.ind], labels = '\\VE', vfont = c('sans serif', 'bold'), col = 'green')
text(x = x.f[nb_m.ind], y = x.b[nb_m.ind], labels = '\\MA', vfont = c('sans serif', 'bold'), col = 'red')
text(x = x.f[nb_f.ind], y = x.b[nb_f.ind], labels = '\\VE', vfont = c('sans serif', 'bold'), col = 'yellow')
legend('topleft', cex = 0.5, legend = c('black female', 'black male', 'nonblack female', 'nonblack male'),
        col = c('green', 'blue', 'yellow', 'red'), pch = 1)


## Compute maginal effects of beta
P = dnorm(X %*% beta, 0, 1)

me.s = beta[2] * P
me.e = beta[3] * P
me.a = beta[4] * P
me.f = beta[5] * P
me.b = beta[6] * P
me.w = beta[7] * P

mean(me.s)
mean(me.e)
mean(me.a)
mean(me.f)
mean(me.b)
mean(me.w)


## Mean of Beta
beta_mean = colMeans(Beta)


## Compare mu_state with the true support rate

# Function for counting number of reverse-pairs 
mergearr = function(arr, mid){
  
  curr = 0
  temparr = c()
  last = length(arr)
  i = 1; j = mid + 1
  count = 0
  
  while (i<=mid & j<=last) {
    curr = curr + 1
    if(arr[i] <= arr[j]){
      temparr[curr] = arr[i]
      i = i + 1
    }
    else{
      temparr[curr] = arr[j]
      j = j + 1
      count = count + mid - i + 1
    }
  }
  
  if(i <= mid){
    while (i <= mid) {
      curr = curr + 1
      temparr[curr] = arr[i]
      i = i + 1
    }
  }
  else{
    while (j <= last) {
      curr = curr + 1
      temparr[curr] = arr[j]
      j = j + 1
    }
  }
  
  result = list(
    "sortarr" = temparr,
    "numrev" = count
  )
  
  return(result)
}


mergesort = function(arr, first, last){
  
  if(first == last){
    numrev = 0
    sortarr = arr[first]
  }
  else{
    mid = floor((first + last) / 2)
    lresult = mergesort(arr, first, mid)
    rresult = mergesort(arr, (mid+1), last)
    
    combarr = c(lresult$sortarr, rresult$sortarr)
    
    mresult = mergearr(combarr, (mid - first + 1))
    
    numrev = lresult$numrev + rresult$numrev + mresult$numrev
    sortarr = mresult$sortarr
  }
  
  result = list(
    "sortarr" = sortarr,
    "numrev" = numrev
  )
  
  return(result)
  
}

mu.s = 0 - as.numeric(paste(unlist(Mu.s[1003, -1])))
rate.s = as.numeric(paste(unlist(Mu.s[1004, -1])))
Mst = as.numeric(summary(polls$state))
ind.s = as.numeric(paste(unlist(Mu.s[1001, -1])))
state = cbind(mu.s, rate.s, Mst, ind.s, levels(polls$state))
state = state[order(mu.s), ]

nInv = mergesort(state[ ,2], 1, 49)

plot(state[ , 1], state[ , 2], xlab = 'mu.s', ylab = 'support rate', pch = 20, cex = 0.8, col = 'blue')
abline(lm(state[ , 2] ~ state[  ,1]), col = 'red')
summary(lm(state[,2] ~ state[,1]))


