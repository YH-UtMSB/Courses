library(MatchIt)

set.seed(1023)
data = lalonde

# To get Fisher stats T for above/under high-school subgroup, run the following
# two lines respectively:
data = subset(lalonde, nodegree==0)
# data = subset(lalonde, nodegree==1)

W = data$treat; Y.obs = data$re78

get_neyman = function(W, Y.obs){
  
  ATE = mean(Y.obs[which(W==1)]) - mean(Y.obs[which(W==0)])
  Sc = sd(Y.obs[which(W==0)])
  St = sd(Y.obs[which(W==1)])
  Vney = Sc^2 / length(which(W==0)) + St^2 / length(which(W==1))
  SDney = sqrt(Vney)
  
  results = c(ATE, SDney)
  names(results) = c("ATE", "SDney")
  
  return(results)
}

est = get_neyman(W, Y.obs)
