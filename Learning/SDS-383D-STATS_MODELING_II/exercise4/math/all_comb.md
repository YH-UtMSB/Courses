All combinations
================
Yilin He

#### It is a dynamic programming example to display all possible combinations of digits from 1 to k.

``` r
all_combin = function(k){
  
  # "Com" stands for "combinations",
  # its rows host all the k_choose_t possible t-element tuples.
  Com = c()
  num_elem = 0
  
  while(num_elem <= k){
    
    if(length(Com) == 0){
      cat("null\n")
      Com = matrix(c(1:k), nrow = k, byrow = TRUE)
    }
    else{
      N = dim(Com)[1]
      for (i in 1:N) {
        cat(Com[i,], "\n")
      }
      
      # num_expd is the number of t+1-element tuple prefixed by the current tuple
      num_expd = k - Com[,num_elem]
      num_derive = sum(num_expd)
      
      # replicate each prefix "num_prod" times, making up the head of t+1-element tuples
      head = unlist(lapply(c(1:N), function(t) cbind(rep(Com[t,], times = num_expd[t]))))
      head = matrix(head, nrow = num_derive, byrow = TRUE)
      
      num_expd = num_expd[which(num_expd != 0)]
      # the collection of "last element"s 
      appender = unlist(lapply(num_expd, function(ee) cbind(c(1:ee))))
      
      
      length(appender) == num_derive
      
      if(num_derive > 0){
        Com = cbind(head, head[,num_elem] + appender)
      }
      else{}
    }
    
    num_elem = num_elem + 1
    
  }
  
}

# See all combinations when k = 4
all_combin(4)
```

    ## null
    ## 1 
    ## 2 
    ## 3 
    ## 4 
    ## 1 2 
    ## 1 3 
    ## 1 4 
    ## 2 3 
    ## 2 4 
    ## 3 4 
    ## 1 2 3 
    ## 1 2 4 
    ## 1 3 4 
    ## 2 3 4 
    ## 1 2 3 4

I'll use a similar code to check performance of all the possible random effects models.
