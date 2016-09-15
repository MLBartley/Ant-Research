############################
##
## 13 June 2016
## Continuous Aspects in MCMC
## Want to consider both number and
## duration of Trophallaxis interactions
##
###########################
#####
# Outline 
#####

  # code (turn into function?) to get data in order 
    # now N_t is number of ants (or pairs of ants) 
    # in trophallaxis at time t

  # code/function to run new mcmc model

  # Simulated data to ensure model works
  
  # apply to ant data



## Get Data in New format

high4 <- read.csv("./Data/Colony1_trophallaxis_high_density_4hr.csv")
low4 <- read.csv("./Data/Colony1_trophallaxis_low_density_4hr.csv")

low4.1 <- low4[which(low4$Location == 1), ]
low4.4 <- low4[which(low4$Location == 4), ]

prep.high = prep.troph.pairs(high4)
prep.low = prep.troph.pairs(low4)
prep.low.1 = prep.troph.pairs(low4.1)
prep.low.4 = prep.troph.pairs(low4.4)

N.high = prep.high$pairs2
N.low = prep.low$pairs2
N.low.1 = prep.low.1$pairs2
N.low.4 = prep.low.4$pairs2


##visualize data
##
##

col = c("#120d08", "#bc5356", "#538bbc", "#53bc84")

par(mfrow = c(2,2))

plot(N.high, main = "High Density", ylab = "# Interactions", xlab = "", type = "l", col = col[2])
plot(N.low, main = "Low Density", ylab = "# Interactions", xlab = "Time (Seconds)", type = "l", col = col[3])
plot(N.low.1, main = "Low Density, Chamber 1", ylab = "# Interactions", xlab = "Time (Seconds)", type = "l", col = col[3])
plot(N.low.4, main = "Low Density, Chamber 2", ylab = "# Interactions", xlab = "Time (Seconds)", type = "l", col = col[3])

par(mfrow = c(1,1))

Time = prep.low$hours * 60 * 60


###
# Simulated Data
###

#see terminal.R


###
# Inference
###

#Propose high and low gamma and lambda, 
  #first simplest case
  # then with additional state based conditions
#calculate R* and then P* matrices
#accept/reject

data = N.low.1[1:(4*60*60)] #only using first two hours for time

#hyperparameters
a = .01
b = .8
c = .04
d =  .8
r = .024
q = .5

#tuning parameter
tau = matrix( c(.05, 0, 0, 0,
                0, .05, 0, 0,
                0, 0, .05, 0, 
                0, 0, 0, .05), nrow = 4, ncol = 4)

n.mcmc = 5000

theta = matrix(c(70000, 1, 1, 70000), 2, 2)

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

params[1:2, 1] = c(.01, .04)/2 #recall, known[1] is gamma.high^tilde
params[3:4, 1] = (.24)/2

X.param[, 1] = sample(c(1, 2), size = Time, replace = T) #1 - high state, 2 - low state


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


for(l in restart:n.mcmc){
  
  # print out every 10 iterations completed
  if( l %% 100 == 0 ) cat(paste("iteration", l, "complete\n")) 
  
  #
  #update
  
  # #adaptive tuning parameter
  # if(l < n.mcmc/2 & l %% 100 == 0){
  #   
  #   sigma = 2.38 ^ 2 / 4 * var(log(t(params[, 1:(l - 1)])))
  #   tau = sigma
  # }
  
  proposal = rmvnorm(n = 1, mean = log(params[, l - 1]), sigma = tau)
  
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
      R_H.star[i, i - 2] = (i - 1)/2 * theta.star[3] #lambda_high
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
      R_L.star[i, i - 2] = (i - 1)/2 * theta.star[4] #lambda_low
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
      R_H[i, i - 2] = (i - 1)/2 * params[3, l - 1] #lambda_high
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
      R_L[i, i - 2] = (i - 1)/2 * params[4, l-1] #lambda_low
    }
    
  }
  
  P_L = Pctmc(Q = R_L, t = 1)
  rownames(P_L) <- 0:max(data)
  colnames(P_L) <- 0:max(data)
  
  
  
  #calculate probability
  MHprob = exp(log.fullcond(theta.star, P_L.star, P_H.star, data, X.param) -
                 log.fullcond(params[, l - 1], P_L, P_H, data, X.param))
  
  if(is.finite(MHprob) == FALSE){MHprob = 0}
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
  
  #
  
  
  
  ##X Parameters
  
  
  m = matrix(data = 0, nrow = 2, ncol = 2)
  rownames(m) <- c("high", "low")
  colnames(m) <- c("high", "low")
  # number states going from i to j, refreshes every run
  
  
  prob_H = P_H[data[1] + 1, data[2] + 1] * .5 *
    M[1, X.param[2, l - 1]]
  
  
  prob_L = P_L[data[1] + 1, data[2] + 1] * .5 *
    M[2, X.param[2, l - 1]]
  
  X.param[1, l] = sample(x = (1:2), size = 1, prob = c(prob_H, prob_L))
  
  m[X.param[1, l], X.param[1, l]] = m[X.param[1, l], X.param[1, l]] + 1
  # 
  # 
  for(t in 2:(Time - 1)){
    
    prob_H = P_H[data[t] + 1, data[t + 1] + 1]  *
      M[1, X.param[t + 1, l - 1]] *
      M[X.param[t - 1, l - 1], 1]
    
    prob_L = P_L[data[t] + 1, data[t + 1] + 1]  *
      M[2, X.param[t + 1, l - 1]] *
      M[X.param[t - 1, l - 1], 2]
    
    X.param[t, l] = sample(x = (1:2), 1,  prob = c(prob_H, prob_L))
    
    m[X.param[t - 1, l], X.param[t, l]] = m[X.param[t - 1, l],
                                            X.param[t, l]] + 1
  }
  
  prob_H = P_H[data[Time - 1] + 1, data[Time] + 1]  *
    M[X.param[1, l - 1], X.param[Time - 1, l - 1]]
  
  prob_L = P_L[data[X.param[Time - 1]] + 1, data[X.param[Time]] + 1]  *
    M[X.param[2, l - 1], X.param[Time - 1, l - 1]]
  
  X.param[Time, l] = sample(x = (1:2), 1,  prob = c(prob_H, prob_L))
  
  m[X.param[Time - 1, l], X.param[Time, l]] = m[X.param[Time - 1, l],
                                                X.param[Time, l]] + 1
  # 
  #M matrix parameter
  
  M[1, ] = (rdirichlet(n = 1 , alpha = theta[1, ] + m[1, ]))
  M[2, ] = (rdirichlet(n = 1 , alpha = theta[2, ] + m[2, ]))
  
  M.param[, l] = as.vector(t(M))
  
}


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


lines(1:n.mcmc, (60 * params[3, ]), col = col[3])
lines(1:n.mcmc, (60 * params[4, ]), col = col[4])

#X params

#Single X
X = X.param[sample(1:Time, 1), ]
plot(0, 0, xlab = "MCMC Runs", ylab = "Single X", ylim = c(0,max(X)), 
     xlim = c(0,n.mcmc), type = "n", cex.lab = 1)
lines(1:n.mcmc, X, col = col[4])

#States over time
plot(X.est, type = "l")
plot(round(X.est), type = 'l')

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
     xlim=c(0,l), 
     type="n",
     cex.lab = 1)
lines(1:(l-1), 60 * (params[1, 1:(l-1)] + params[2, 1:(l-1)]), col = col[1])
lines(1:(l - 1), 60 * params[2, 1:(l-1)], col = col[2])

lines(1:(l-1), (60 * params[3, 1:(l-1)]), col = col[3])
lines(1:(l-1), (60 * params[4, 1:(l-1)]), col = col[3])


#States over time
#
X.est = matrix(data = rep(NA, Time), nrow = Time, ncol = 1)

for(t in 1:Time ){
  X.est[t, 1] = mean(X.param[t, 1:(l-1)])  
}


plot(X.est, type = "l")
plot(round(X.est), type = "l")

#M
plot(0,0,xlab="MCMC Runs", ylab = "M", ylim = c(0, max(M.param[1:(l-1)])), xlim=c(0,n.mcmc), 
     type="n", cex.lab = 1)
for(i in 1:(4)){
  lines(1:(l-1), M.param[i, 1:(l-1)], col = col[i])
}

