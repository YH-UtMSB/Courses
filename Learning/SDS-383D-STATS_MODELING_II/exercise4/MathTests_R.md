## Math Tests

#### Part I

This part is all about examining the phenomenon that the extreme school-level averages tend to be at schools where fewer students were sampled.


```r
# mathtest problem
suppressMessages('mosaic')
library(readr)
library(mosaic)

datapath = 'E:/R/mathtest.csv'
mathtest = read_csv(datapath)
# order the data by school
mathtest = mathtest[order(mathtest$school),]

# collect schools
sch = unique(mathtest$school)

# Ensure that school is treated as a categorical variable
mathtest$school = factor(mathtest$school)

# This won't work without mosaic loaded
schoolmeans = mean(mathscore~school, data=mathtest)
schoolsizes = as.numeric(summary(mathtest$school))

# take the extreme point
maxmean = max(schoolmeans)
maxmean.size = schoolsizes[which.max(schoolmeans)]

# problem 1 -- Notice the extremes tend to be for schools with low sample sizes
plot(schoolsizes, schoolmeans, pch="+")
points(maxmean.size, maxmean)

```
![](fig/ex04-1-1-fig1.jpeg)
