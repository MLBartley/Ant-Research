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
known = c(.06, .01, .02, .005)

#calculate P matrices for KNOWN values - for high/low states

R_H = matrix(0, nrow = 14 + 1, ncol = 14 + 1)
rownames(R_H) <- 0:14
colnames(R_H) <- 0:14

for(i in 1:nrow(R_H)){
  
  if( i %% 2 != 0 & i != nrow(R_H) & i != (nrow(R_H) - 1)){
    R_H[i, i + 2] = known[1]  #gamma_high
  }
  
  if( i %% 2 != 0 & i != 1 & i != 2){
    R_H[i, i - 2] = (i - 1) * known[3] #lambda_high 
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
    R_L[i, i - 2] = (i - 1) * known[4] #lambda_low 
  }
  
}

P_L.known = Pctmc(Q = R_L, t = 1)
rownames(P_L.known) <- 0:14
colnames(P_L.known) <- 0:14


#use P matrices to sample data, N 
N = rep(NA, Time)

N[1] = 0

for(i in 2:Time){
  if(x[i-1] == 1){
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
a = .06
b = .8
r = .01
q = .5

#tuning parameter
tau = rep(.01, 4)

n.mcmc = 1000

theta = matrix(c(999, 1, 1, 999), 2, 2)

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
M = matrix(c(0.995, 0.005, 0.005, 0.995),byrow = T,  2, 2)
rownames(M) <- c("High", "Low")
colnames(M) <- c("High", "Low")


#initialize

params[1:2, 1] = c(known[1] - known[2], known[2])/2
params[3:4, 1] = known[3:4]/2

X.param[, 1] = rep(1, Time) #1 - high state, 2 - low state

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
    dgamma(params[2], a, b, log = T) + 
    dgamma(params[3], r, q, log = T) +
    dgamma(params[4], r, q, log = T)
  
  return(loglike)
}

accept = 0

for(l in 2:n.mcmc){
  
  # print out every 10 iterations completed
  if( l %% 10 == 0 ) cat(paste("iteration", l, "complete\n")) 
  
  #
  #update

  #adaptive tuning parameter
  if(l < n.mcmc/2 & l %% 100 == 0){

    sigma = c(0, 0, 0, 0)
    for(v in 1:4){
      sigma[v] = ((2.38 ^ 2) / 4) * var(log(params[v, 1:(l - 1)]))
    }

    if(sigma[1] != 0){
      tau[1] = sigma[1]
    } 
    if(sigma[2] != 0){
      tau[2] = sigma[2]
    } 
    if(sigma[3] != 0){
      tau[3] = sigma[3]
    } 
    if(sigma[4] != 0){
      tau[4] = sigma[4]
    }

  }
   
  proposal = rnorm(4, mean = log(params[, l - 1]), sd = tau)

  theta.star = exp(proposal)
  
  gamma.high = theta.star[1] + theta.star[2] #calculate correct gamma_high
  
  
  #calculate P* matrices - for high/low states
  
  R_H.star = matrix(0, nrow = max(data) + 1, ncol = max(data) + 1)
  rownames(R_H.star) <- 0:max(data)
  colnames(R_H.star) <- 0:max(data)
  
  for(i in 1:nrow(R_H.star)){
    
    if( i %% 2 != 0 & i != nrow(R_H.star) & i != (nrow(R_H.star) - 1)){
      R_H.star[i, i + 2] = gamma.high  #gamma_high
    }
    
    if( i %% 2 != 0 & i != 1 & i != 2){
      R_H.star[i, i - 2] = (i - 1) * theta.star[3] #lambda_high 
    }
    
  }
  
  P_H.star = Pctmc(Q = R_H.star, t = 1)
  rownames(P_H.star) <- 0:max(data)
  colnames(P_H.star) <- 0:max(data)
  
  
  R_L.star = matrix(0, nrow = max(data) + 1, ncol = max(data) + 1)
  rownames(R_L.star) <- 0:max(data)
  colnames(R_L.star) <- 0:max(data)
  
  for(i in 1:nrow(R_L.star)){
    
    if( i %% 2 != 0 & i != nrow(R_L.star) & i != (nrow(R_L.star) - 1)){
      R_L.star[i, i + 2] = theta.star[2]  #gamma_low
    }
    
    if( i %% 2 != 0 & i != 1 & i != 2){
      R_L.star[i, i - 2] = (i - 1) * theta.star[4] #lambda_low 
    }
    
  }
  
  P_L.star = Pctmc(Q = R_L.star, t = 1)
  rownames(P_L.star) <- 0:max(data)
  colnames(P_L.star) <- 0:max(data)
  
  
  
  #calculate P matrices for PREVIOUS values - for high/low states
  
  R_H = matrix(0, nrow = max(data) + 1, ncol = max(data) + 1)
  rownames(R_H) <- 0:max(data)
  colnames(R_H) <- 0:max(data)
  
  for(i in 1:nrow(R_H)){
    
    if( i %% 2 != 0 & i != nrow(R_H) & i != (nrow(R_H) - 1)){
      R_H[i, i + 2] = params[1, l - 1] + params[2, l - 1]  #gamma_high
    }
    
    if( i %% 2 != 0 & i != 1 & i != 2){
      R_H[i, i - 2] = (i - 1) * params[3, l - 1] #lambda_high 
    }
    
  }
  
  P_H = Pctmc(Q = R_H, t = 1)
  rownames(P_H) <- 0:max(data)
  colnames(P_H) <- 0:max(data)
  
  R_L = matrix(0, nrow = max(data) + 1, ncol = max(data) + 1)
  rownames(R_L) <- 0:max(data)
  colnames(R_L) <- 0:max(data)
  
  for(i in 1:nrow(R_L)){
    
    if( i %% 2 != 0 & i != nrow(R_L) & i != (nrow(R_L) - 1)){
      R_L[i, i + 2] = params[2, l-1]  #gamma_low
    }
    
    if( i %% 2 != 0 & i != 1 & i != 2){
      R_L[i, i - 2] = (i - 1) * params[4, l-1] #lambda_low 
    }
    
  }
  
  P_L = Pctmc(Q = R_L, t = 1)
  rownames(P_L) <- 0:max(data)
  colnames(P_L) <- 0:max(data)
  
  
  
  #calculate probability
  MHprob = exp(log.fullcond(theta.star, P_L.star, P_H.star, data, X.param) -
                 log.fullcond(params[, l - 1], P_L, P_H, data, X.param))
  
  
  #accept/reject
  
  if(runif(1) < MHprob){
    accept = accept + 1
    params[, l] = theta.star[]
    #needed for X
    P_H = P_H.star
    P_L = P_L.star
  }else{
    params[, l] = params[, l - 1]
  }
  
  
  ##X Parameters
  
  ### simulation - keep same
  
  #X.param[, l] = x #comment out when simulating X
  ###################
  
  
  m = matrix(data = 0, nrow = 2, ncol = 2)
  # number states going from i to j, refreshes every run


  prob_H = P_H[data[1] + 1, data[2] + 1] * .5 *
    M[X.param[1, l - 1], X.param[2, l - 1]]


  prob_L = P_L[data[1] + 1, data[2] + 1] * .5 *
    M[X.param[1, l - 1], X.param[2, l - 1]]

  X.param[1, l] = sample(x = (1:2), size = 1, prob = c(prob_H, prob_L))

  m[X.param[1, l], X.param[1, l]] = m[X.param[1, l], X.param[1, l]] + 1
  # 
  # 
  for(t in 2:(Time - 1)){

    prob_H = P_H[data[X.param[t-1]] + 1, data[X.param[t]] + 1]  *
      M[X.param[t, l - 1], X.param[t-1, l-1]] *
      M[X.param[t + 1, l - 1], X.param[t, l-1]]

    prob_L = P_L[data[X.param[t-1]] + 1, data[X.param[t]] + 1]  *
      M[X.param[t, l - 1], X.param[t-1, l-1]] *
      M[X.param[t + 1, l - 1], X.param[t, l-1]]

   X.param[t, l] = sample(x = (1:2), 1,  prob = c(prob_H, prob_L))

    m[X.param[t - 1, l], X.param[t, l]] = m[X.param[t - 1, l],
                                            X.param[t,l]] + 1
  }

  prob_H = P_H[data[X.param[Time - 1]] + 1, data[X.param[Time]] + 1]  *
    M[X.param[Time, l - 1], X.param[Time - 1, l - 1]]

  prob_L = P_L[data[X.param[Time - 1]] + 1, data[X.param[Time]] + 1]  *
    M[X.param[Time, l - 1], X.param[Time - 1, l - 1]]

   X.param[Time, l] = sample(x = (1:2), 1,  prob = c(prob_H, prob_L))

  m[X.param[Time - 1, l], X.param[Time, l]] = m[X.param[Time - 1, l],
                                                X.param[Time, l]] + 1
  # 
  #M matrix parameter

  # M[1, ] = (rdirichlet(n = 1 , alpha = theta[1, ] + m[1, ]))
  # M[2, ] = (rdirichlet(n = 1 , alpha = theta[2, ] + m[2, ]))

  M.param[, l] = as.vector(t(M))
  
}




#compile estimates
X.est = matrix(data = rep(NA, Time), nrow = Time, ncol = 1)


gamma.high.est = mean(params[1, ] + params[2, ])
gamma.low.est = mean(params[2, ])
lambda.high.est = mean(params[3, ])
lambda.low.est = mean(params[4, ])

estimate = c(gamma.high.est, gamma.low.est, lambda.high.est, lambda.low.est)

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


#plot the estimation runs.

col = c("#120d08", "#bc5356", "#538bbc", "#53bc84")

pdf(file = paste("./simulation-practice/output/", Sys.time(), ".pdf", sep = ""))

#gamma
plot(0,0,xlab="MCMC Runs",
     ylab="Gamma (per minute)",
     ylim=c(0,(max(params[1, ] + params[2, ]) * 60)), 
     xlim=c(0,n.mcmc), 
     type="n",
     cex.lab = 1)
lines(1:n.mcmc, 60 * (params[1, ] + params[2, ]), col = col[1])
lines(1:n.mcmc, 60 * params[2, ], col = col[2])
abline(h = known[1] * 60)
abline(h = known[2] * 60)

 #lambda

plot(0,0,xlab = "MCMC Runs", ylab = "Lambda (per minute)", 
     ylim = c(0, max(params) * 60), xlim = c(0, n.mcmc), 
     type = "n", cex.lab = 1)

lines(1:n.mcmc, (60 * params[3, ]), col = col[3])
lines(1:n.mcmc, (60 * params[4, ]), col = col[4])
abline(h = known[3] * 60)
abline(h = known[4] * 60)

#X params

#Single X
X = X.param[sample(1:Time, 1), ]
plot(0, 0, xlab = "MCMC Runs", ylab = "Single X", ylim = c(0,max(X)), 
     xlim = c(0,n.mcmc), type = "n", cex.lab = 1)
lines(1:n.mcmc, X, col = col[4])

#States over time
plot(X.est, type = "l", lwd = 3, cex.lab = 1, col = col[1])


#M
plot(0,0,xlab="MCMC Runs", ylab = "M", ylim=c(0, max(M.param)), xlim=c(0,n.mcmc), 
     type="n", cex.lab = 1)
for(i in 1:(4)){
  lines(1:n.mcmc, M.param[i, ], col = col[i])
}

dev.off()

