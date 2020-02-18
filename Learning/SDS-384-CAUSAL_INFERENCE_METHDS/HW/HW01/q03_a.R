library(MatchIt)
set.seed(1023)

data = lalonde

# To get Fisher stats T for above/under high-school subgroup, run the following
# two lines respectively:
data = subset(lalonde, nodegree==0)
# data = subset(lalonde, nodegree==1)

W = data$treat; Y.obs = data$re78; X = data$re74
N = length(W); Nt = sum(W)

# the number of simulations
n_sim = 1000

Y.mis = Y.obs
Y1 = W * Y.obs + (1. - W) * Y.mis
Y0 = (1. - W) * Y.obs + W * Y.mis


get_stat = function(Wl, Y1, Y0, C = 0, stat="dif"){
  Y.obs = Wl * Y1 + (1. - Wl) * Y0
  
  if(stat == 'dif'){
    Yt = Y.obs[which(Wl==1)]
    Yc = Y.obs[which(Wl==0)]
    return(abs(mean(Yt) - mean(Yc) - C))
  }
  else if(stat == 'rank')
    V.obs = Y.obs
  else{
    G.obs = Y.obs - X
    V.obs = G.obs
  }
  
  # obtain rank stat
  value = sort(unique(V.obs))
  counts = matrix(table(V.obs))
  counts.l = cumsum(counts) - counts
  R.set = counts.l + 0.5 * (counts - N) 
  R = c()
  for(i in 1:N){
    R = c(R, R.set[which(value == V.obs[i])])
  }
  Rt = R[which(Wl==1)]
  Rc = R[which(Wl==0)]
  
  tr = abs(mean(Rt) - mean(Rc))
  
  return(tr)
}

# the 3 stats associated with currently observed assignments & outcomes
r_stat = get_stat(W, Y1, Y0, stat = 'rank')
rg_stat = get_stat(W, Y1, Y0, stat = 'rgain')
dif_stat = get_stat(W,Y1, Y0, stat = 'dif')

t.r = t.rg = t.dif = c()

for (l in 1:n_sim) {
  
  Wl = sample(W)
  tr = get_stat(Wl, Y1, Y0, stat = 'rank')
  trg = get_stat(Wl, Y1, Y0, stat = 'rgain')
  tdif = get_stat(Wl, Y1, Y0, stat = 'dif')
  
  t.r = c(t.r, tr)
  t.rg = c(t.rg, trg)
  t.dif = c(t.dif, tdif)
}


p.tr = length(which(t.r>r_stat)) / n_sim
p.trg = length(which(t.rg>rg_stat)) / n_sim
p.dif = length(which(t.dif>dif_stat)) / n_sim