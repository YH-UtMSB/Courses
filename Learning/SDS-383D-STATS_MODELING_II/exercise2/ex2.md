Programming Practices in Exercise 2
================
Yilin He

R Markdown
----------

``` r
# this line disables messages upon loading package 'mosaic'
suppressMessages(library(mosaic))
library(readr)
# library(mosaic)

gdpgrowth = read_csv("E:/R/gdpgrowth.csv")
```

    ## Parsed with column specification:
    ## cols(
    ##   CODE = col_character(),
    ##   COUNTRY = col_character(),
    ##   GR6096 = col_double(),
    ##   DENS60 = col_double(),
    ##   COAST65 = col_double(),
    ##   POPGR6090 = col_double(),
    ##   EAST = col_integer(),
    ##   DEF60 = col_double(),
    ##   LGDP60 = col_double(),
    ##   EDUC60 = col_double(),
    ##   LIFE60 = col_double()
    ## )

``` r
y = gdpgrowth$GR6096; n = length(y)

x0 = rep(1., n); x1 = gdpgrowth$DEF60
X = matrix(c(x0, x1), ncol = 2, byrow = FALSE)

# set up prior parameters: Lambda, K, m
L = diag(n); K = diag(2); m = c(0.3,0.2)

L.star = t(X) %*% L %*% X + K
m.star = solve(L.star) %*% (t(X) %*% L %*% y + K %*% m)

m.star
```

    ##            [,1]
    ## [1,] 0.01571719
    ## [2,] 0.19318186

``` r
coef(lm(y~x1))
```

    ## (Intercept)          x1 
    ##  0.01176803  0.20650591

``` r
plot(x1, y)
abline(m.star, col = "red")
abline(lm(y~x1))
```

![](ex2_files/figure-markdown_github/unnamed-chunk-1-1.png)
