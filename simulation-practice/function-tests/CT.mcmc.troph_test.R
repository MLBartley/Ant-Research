##
## Simulation Study = CTMC 
## parameters: {gamma_h, gamma_l, lambda_h, lambda_l, X_t, M}
## 20 June 2016
##
##########################################################

## needed for remote access, not sure if needed every time?
chooseCRANmirror(ind = 27)
install.packages("ctmcmove", dependencies = T)
library("ctmcmove")
library("gtools")
library("mvtnorm")


Time = 1 * 60 * 60

#two state transition probabilities for X_t
M = matrix(c(.99, .01, .01, .99), 2, 2)


#high/low trophallaxis rate states
x = rep(NA,T)
x[1] = 1

for(t in 2:Time){
  ## sample latent state
  x[t] = sample(1:2, 1, prob = M[x[t - 1],])

}

#check
plot(x, type = 'l')

#choose known parameters for gammas and lambdas
known = c(.15, .2, .1, .3)
gamma.high.known = known[1] + known[2]

#calculate P matrices for KNOWN values - for high/low states

R_H = matrix(0, nrow = 14 + 1, ncol = 14 + 1)
rownames(R_H) <- 0:14
colnames(R_H) <- 0:14

for(i in 1:nrow(R_H)){
  
  if( i %% 2 != 0 & i != nrow(R_H) & i != (nrow(R_H) - 1)){
    R_H[i, i + 2] = gamma.high.known  #gamma_high
  }
  
  if( i %% 2 != 0 & i != 1 & i != 2){
    R_H[i, i - 2] = (i - 1)/2 * known[3] #lambda_high 
  }
  
}

P_H.known = Pctmc(Q = R_H, t = 1)
rownames(P_H.known) <- 0:14
colnames(P_H.known) <- 0:14


R_L = matrix(0, nrow = 14 + 1, ncol = 14 + 1)
rownames(R_L) <- 0:14
colnames(R_L) <- 0:14

for(i in 1:nrow(R_L)){
  
  if( i %% 2 != 0 & i != nrow(R_L) & i != (nrow(R_L) - 1)){
    R_L[i, i + 2] = known[2]  #gamma_low
  }
  
  if( i %% 2 != 0 & i != 1 & i != 2){
    R_L[i, i - 2] = (i - 1)/2 * known[4] #lambda_low 
  }
  
}

P_L.known = Pctmc(Q = R_L, t = 1)
rownames(P_L.known) <- 0:14
colnames(P_L.known) <- 0:14


#use P matrices to sample data, N 
N = rep(NA, Time)

N[1] = 0

for(i in 2:Time){
  if(x[i] == 1){
    N[i] = sample(0:14, 1, prob = P_H.known[N[i-1] + 1, ])
  }else{
    N[i] = sample(0:14, 1, prob = P_L.known[N[i-1] + 1, ])
  }
 
}

#visualize, looks alright
plot(N, type = 'l')
points(x, col = 'red')

data = N

#hyperparameters
a = .1
b = .8
c = .2
d =  .8
r = .2
q = .5

#tuning parameter
tau = matrix( c(.05, 0, 0, 0,
                0, .05, 0, 0,
                0, 0, .05, 0, 
                0, 0, 0, .05), nrow = 4, ncol = 4)

n.mcmc = 10000

theta = matrix(c(9999, 1, 1, 9999), 2, 2)

test1 = CT.mcmc.troph()