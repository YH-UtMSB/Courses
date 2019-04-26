## Cheese-elastic Case

Package requirment:

```r
library(lme4)
library(tidyverse)
library(lattice)
library(readr)

rm(list=ls())
```

To start with, I fit a Bayesian linear model as well as several MLR models through OLS. I expect the Bayesian performs at least among the top 30% of the models fitted in the frequentist fashion. So I split the cheese data set into a bigger training half and a smaller evaluation half. By inspecting the data size for each store, I fix the training sample size for each store at 50. Namely, I gather the first 50 data for each store to form the training set, and what remain become the eval set.

```r

cheese = read_csv("E:/R/cheese.csv")

cheese$store = factor(cheese$store)
cheese = cheese[order(cheese$store),]

# collect the number of data point for each store
N = as.numeric(summary(cheese$store))
n = length(N)


# Gather the first 50 data in each store to form the training set
# by specifying starting points for each store 
# and count 49 samples afterwards for each starting point.
edpt = cumsum(N)
stpt = c(1,(edpt+1)[-n])

m = 50
tr_idx = c(matrix(rep(stpt, times = m), nrow = m, byrow = TRUE) + (0:(m-1)))


# Split dataset into training set and eval set for OLS model.
# We're about to train the model on the training set and evaluate its performance on the eval set.
# Specifically, we're going to compare the performace of the model fitted via gibbs sampling and 
# the ones by solving OLS. 
cheese.tr = cheese[tr_idx,]
cheese.ev = cheese[-tr_idx,]


# Prepare training & eval data for bayesian model.
# Collecting training data is just a walk in the park...
y = log(cheese.tr$vol); x = log(cheese.tr$price); D = cheese.tr$disp
y = matrix(y, nrow = n, byrow = TRUE)
x = matrix(x, nrow = n, byrow = TRUE)
D = matrix(D, nrow = n, byrow = TRUE)

```

Since "rectangulating" y's, x's and D's in the evaluation into matrices makes life a lot easier, I append 0's to the short rows to make them as long as the longest.  

```r

```
