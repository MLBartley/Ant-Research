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


Time = 4 * 60 * 60

#two state transition probabilities for X_t
M = matrix(c(.99, .01, .99, .01), 2, 2)


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
known = c(0, .0325, .016, .016)
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
#
high4 <- read.csv("./Data/Colony1_trophallaxis_high_density_4hr.csv")
prep.high = prep.troph.pairs(high4)

N.high = prep.high$pairs2
par(mfrow = c(2,1))

plot(N.high, main = "High Density", ylab = "# Interactions", xlab = "", type = "l", col = col[2])
plot(N, type = 'l',main = "Simulated Data", xlab = '', ylab = "# Interaction ")

par(mfrow = c(1,1))


points(x, col = 'red')

data = N

#hyperparameters
a = .001
b = .8
c = .03
d =  .8
r = .016
q = .5

#tuning parameter
tau = matrix( c(.05, 0, 0, 0,
                0, .05, 0, 0,
                0, 0, .05, 0, 
                0, 0, 0, .05), nrow = 4, ncol = 4)

n.mcmc = 5000

theta = matrix(c(800, 200, 200, 800), 2, 2)

#homes

params = matrix(NA, nrow = 4, ncol = n.mcmc)
rownames(params) <- c("gamma_high^tilde", "gamma_low", "lambda_high", "lambda_low")
colnames(params) <- 1:n.mcmc

X.param = matrix(NA, nrow = Time, 
                 ncol = n.mcmc, byrow = T)

# probability matrix for 2 state discrete time Markov Chain

#stores all M matrices for every iteration in mcmc
M.param = matrix(data = NA, nrow = 4, 
                 ncol = n.mcmc, byrow = T)
rownames(M.param) <- c("HH", "HL", "LH", "LL")
colnames(M.param) <- 1:n.mcmc

#is replaced every loop
M = matrix(c(0.99, 0.01, 0.01, 0.99),byrow = T,  2, 2)
rownames(M) <- c("High", "Low")
colnames(M) <- c("High", "Low")


#initialize

params[1:2, 1] = known[1:2]/2 #recall, known[1] is gamma.high^tilde
params[3:4, 1] = known[3:4]/2

X.param[, 1] = x #1 - high state, 2 - low state

M.param[, 1] = as.vector(t(M))


log.fullcond = function(params, P_L, P_H, data, X.param){
  
  sumP_H = 0
  sumP_L = 0
  
  for(t in 2:length(data)){
    if(X.param[t, l-1] == 1){
      sumP_H = sumP_H + log(P_H[data[t - 1] + 1, data[t] + 1])
    }else{
      sumP_L = sumP_L + log(P_L[data[t - 1] + 1, data[t] + 1])
    }
  }
  
  
  loglike = sumP_L + 
    sumP_H +
    dgamma(params[1], a, b, log = T) + 
    dgamma(params[2], c, d, log = T) + 
    dgamma(params[3], r, q, log = T) +
    dgamma(params[4], r, q, log = T)
  
  
  
  return(loglike)
}

accept = 0
sigma = NA


 

#compile estimates - stop code midway version below
X.est = matrix(data = rep(NA, Time), nrow = Time, ncol = 1)


gamma.high.tilde.est = mean(params[1, ])
gamma.low.est = mean(params[2, ])
lambda.high.est = mean(params[3, ])
lambda.low.est = mean(params[4, ])

gamma.high.est = mean(params[1,]) + mean(params[2])

estimate = c(gamma.high.tilde.est, gamma.low.est, lambda.high.est, lambda.low.est, gamma.high.est)

for(t in 1:Time ){
  X.est[t, 1] = mean(X.param[t, ])  
}

###
###
### Results
###
###

#######
#accept ratio and estimated vs known
#######

accept/n.mcmc
known
estimate
round(estimate, digits = 3)

#plot the estimation runs.

col = c("#120d08", "#bc5356", "#538bbc", "#53bc84")


pdf(file = paste("./output/", Sys.time(), ".pdf", sep = ""))

#Parameters
plot(0,0,xlab="MCMC Runs",
     ylab="Rates (per minute)",
     ylim=c(0,(max(params[1, ] + params[2, ]) * 60)), 
     xlim=c(0,n.mcmc), 
     type="n",
     cex.lab = 1)
lines(1:n.mcmc, 60 * (params[1, ] + params[2, ]), col = col[1])
lines(1:n.mcmc, 60 * params[2, ], col = col[2])
abline(h = (known[1] + known[2]) * 60, col = col[1])
abline(h = known[2] * 60, col = col[2])

lines(1:n.mcmc, (60 * params[3, ]), col = col[3])
lines(1:n.mcmc, (60 * params[4, ]), col = col[4])
abline(h = known[3] * 60, col = col[3])
abline(h = known[4] * 60, col = col[4])

#X params

#Single X
X = X.param[sample(1:Time, 1), ]
plot(0, 0, xlab = "MCMC Runs", ylab = "Single X", ylim = c(0,max(X)), 
     xlim = c(0,n.mcmc), type = "n", cex.lab = 1)
lines(1:n.mcmc, X, col = col[4])

#States over time
plot(x, type = "l")
lines(X.est, type = "l", cex.lab = 1, col = col[2])


#M
plot(0,0,xlab="MCMC Runs", ylab = "M", ylim = c(0, max(M.param)), xlim=c(0,n.mcmc), 
     type="n", cex.lab = 1)
for(i in 1:(4)){
  lines(1:n.mcmc, M.param[i, ], col = col[i])
}

dev.off()




###
### Code for visuals midway through running - helpful for debugging
###

col = c("#120d08", "#bc5356", "#538bbc", "#53bc84")


#Parameters
 plot(0,0,xlab="MCMC Runs",
     ylab="Rates (per minute)",
     ylim=c(0,(max(params[1, 1:(l-1)] + params[2, 1:(l-1)]) * 60)), 
     xlim=c(0,n.mcmc), 
     type="n",
     cex.lab = 1)
lines(1:(l-1), 60 * (params[1, 1:(l-1)] + params[2, 1:(l-1)]), col = col[1])
lines(1:(l - 1), 60 * params[2, 1:(l-1)], col = col[2])
abline(h = (known[1] + known[2]) * 60, col = col[1])
abline(h = known[2] * 60, col = col[2])

lines(1:(l-1), (60 * params[3, 1:(l-1)]), col = col[3])
lines(1:(l-1), (60 * params[4, 1:(l-1) ]), col = col[4])
abline(h = known[3] * 60, col = col[3])
abline(h = known[4] * 60, col = col[4])

#States over time
#
X.est = matrix(data = rep(NA, Time), nrow = Time, ncol = 1)

for(t in 1:Time ){
  X.est[t, 1] = mean(X.param[t, 1:(l-1)])  
}


plot(x, type = "l")
lines(X.est, type = "l", cex.lab = 1, col = col[2])


#M
plot(0,0,xlab="MCMC Runs", ylab = "M", ylim = c(0, max(M.param[1:(l-1)])), xlim=c(0,n.mcmc), 
     type="n", cex.lab = 1)
for(i in 1:(4)){
  lines(1:(l-1), M.param[i, 1:(l-1)], col = col[i])
}
abline(h = c(.8, .2))

