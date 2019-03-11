Programming Practices in Exercise 3
================
Yilin He

### Curve fitting by linear smoothing -- part (B)

Basic setup:

1.  *x*'s are sampled from 𝒩(0, 4); *f*(*x*)=*x*<sup>3</sup> + *x*<sup>2</sup> + *x*.
2.  Error *ε* is sampled from 𝒩(0, 1).
3.  Bandwidth: try *h* = .25, .5, 1, 2.

Black dots are the training data, red curves are fitted values.

``` r
# Bandwidth: h = .5, 1, 2, 4, 8
Bw = c(.25, .5, 1, 2)

# ground truth
f = function(x){
  return(x^3 + x^2 + x)
}

# prepare dataset
N = 500
x.tr = rnorm(N, mean = 0, sd = 2)
x.tr = x.tr - mean(x.tr)
e = rnorm(N, mean = 0, sd = 1)
y.tr = f(x.tr) + e

x.fit = seq(-6, 6, .1)

par(mfrow=c(2,2))
for (h in Bw) {
  
  # compute weights; weights are stored in matrix W
  # W[i,j] = K((x.fit[i] - x.tr[j])/h) / norm_const
  W = matrix(rep(x.tr, length(x.fit)), nrow = length(x.fit), byrow = TRUE)
  W = exp( - 0.5 * ((W - x.fit) / h) ^ 2 )
  W = W / apply(W, 1, sum)
  
  # compute integrated y's
  y.fit = W %*% y.tr
  
  plot(x.tr, y.tr, pch = ".", main = paste0("h = ", h))
  lines(x.fit, y.fit, col="red")
  
}
```

![](exercise03_R_files/figure-markdown_github/unnamed-chunk-1-1.png)

``` r
par(mfrow=c(1,1))
```
