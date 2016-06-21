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

  # code (function?) to get data in order 
    # now N_t is number of ants (or pairs of ants) 
    # in trophallaxis at time t

  # function to run new mcmc model

  # Simulated data to ensure model works
  
  # apply to ant data



## Get Data in New format

high4 <- read.csv("./Data/Colony1_trophallaxis_high_density_4hr.csv")

prep.high = prep.troph.pairs(high4)

N = prep.high$pairs2
Time = prep.high$hours * 60 * 60



###
# Inference
###

#Propose gamma, lambda, X_t
  #first simplest case
  # then with additional state based conditions
#calculate R, Q and then P matrices
#accept/reject

n.mcmc = 100
a = 0.03
b = 0.03
r = 0.03
q = 0.03
tau = c(0.02, 0.02)

data = N

# Need probability matrix for 2 state discrete time Markov Chain

DTMC.prob = matrix(c(0.99, 0.01, 0.99, 0.01), 2, 2)

#homes

params = matrix(NA, nrow = 2, ncol = n.mcmc)
rownames(params) <- c("gamma", "lambda")
colnames(params) <- 1:n.mcmc

X.param = matrix(NA, nrow = Time, 
                 ncol = n.mcmc, byrow = T)


#initialize

params[1, 1] = rgamma(n = 1, shape = a, rate = b)
params[2, 1] = rgamma(n = 1, shape = r, rate = q)

X.param[, 1] = rep(1, Time)


# #posterior likelihood
# 
# log.fullcond = function(params, P, data){
#   
#   sumP = 0 
#   
#   for(t in 2:length(data)){
#     sumP = sumP + log(P[data[t - 1] + 1, data[t] + 1])
#   }
#   
#   loglike = sumP + dgamma(params[1], a, b, log = T) + dgamma(params[2], r, q, log = T)
#   
#   return(loglike)
# }

#MCMC

for(l in 2:n.mcmc){
  
  #
  # Calculate indicator values for parameters
  #
  
  N.increase = rep(0, Time)

  for(t in 2:Time){
    if(N[t] !< N[t-1]){
    N.increase[t] = 1}
  }
  
  N.decrease = rep(0, Time)
  
  for(t in 2:Time){
    if(N[t] !> N[t-1]){
      N.decrease[t] = 1}
  }
  
    m = matrix(data = rep(0, n * n), nrow = n, ncol = n) 
  # number states going from i to j, refreshes every run
  
 I_lambda = sum(N.decrease)
 
 I_gamma = c(0, 0)
 I_gamma[1] = sum(N.increase[which(X.param[, l] == 1)])
   
 I_gamma[2] = sum(N.increase[which(X.param[, l] == 2)])   
 
  #
  #update
  #
  
  # proposal = rexp(n = 2, rate = 1 / params[, l - 1])

  params[1, l] = dgamma( (I_gamma[X.param[l-1]] + a), b)
  params[2, l] = dgamma( (I_lambda + r), q)
  

  ##X Parameters
  for(k in 1:n){
    gam[1, k] = lambda.param[k, l - 1] ^ data[1] * exp(-lambda.param[k, l - 1]) * 
      delta[k] * P.matrix[k, X.param[2, l - 1]]
  }
  
  X.param[1, l] = sample(x = (1:n), size = 1, prob = gam[1,])
  

  for(t in 2:(Time - 1)){
    
    for(k in 1:n){
      gam[t, k] = (lambda.param[k, l - 1] ^ data[t]) * 
        exp(-lambda.param[k, l - 1]) * P.matrix[X.param[t - 1, l - 1], k] *
        P.matrix[k, X.param[t + 1, l - 1]]
    }
    
    X.param[t, l] = sample(x = (1:n), 1,  prob = gam[t, ]) 
    
    m[X.param[t - 1, l], X.param[t, l]] = m[X.param[t - 1, l], 
                                            X.param[t,l]] + 1
  }
  
  for(k in 1:n){
    gam[Time, k] = lambda.param[k, l - 1] ^ data[Time] * 
      exp(-lambda.param[k, l - 1]) * P.matrix[X.param[(Time - 1), l - 1], k]
  }
  
  X.param[Time, l] = sample(x = 1:n, 1,  prob = gam[Time, ])
  
  m[X.param[Time - 1, l], X.param[Time, l]] = m[X.param[Time - 1, l], 
                                                X.param[Time, l]] + 1
  
  
#calculate P matrix
  
  R = matrix(0, nrow = max(N) + 1, ncol = max(N) + 1)
    rownames(R) <- 0:max(N)
    colnames(R) <- 0:max(N)
    
  for( i in 1:nrow(R)){
    
    if( i %% 2 != 0 & i != nrow(R) & i != (nrow(R) - 1)){
       R[i, i + 2] = proposal[1]  
    }
    
    if( i %% 2 != 0 & i != 1 & i != 2){
        R[i, i - 2] = (i - 1) * proposal[2] 
    }
    
  }
    
  P = Pctmc(Q = R, t = 1)
    rownames(P) <- 0:max(N)
    colnames(P) <- 0:max(N)

#   #calculate likelihood
# 
#   prob.params = exp(log.fullcond(proposal, P, data) - log.fullcond(params[, l - 1], P, data))
#     
#   #accept/reject
# 
#   if(runif(1) < prob.params){
#     params[, l] = proposal
#   }else{
#    params[, l] = params[, l - 1]
#   }
#   
#   
# }

  #compile estimates

gamma.est = mean(params[1, ])
lambda.est = mean(params[2, ])


  #plot the estimation runs.

col = c("#120d08", "#bc5356", "#538bbc", "#53bc84")

#gamma
plot(0,0,xlab="MCMC Runs",
     ylab="Gamma (per minute)",
     ylim=c(0,max(params)), 
     xlim=c(0,n.mcmc), 
     type="n",
     cex.lab = 1)
  lines(1:n.mcmc, 60 * params[1, ], col = col[1])



#lambda
plot(0,0,xlab = "MCMC Runs", ylab = "Lambda (per minute)", 
     ylim = c(0, max(params)), xlim = c(0, n.mcmc), 
     type = "n", cex.lab = 1)

  lines(1:n.mcmc, (60 * params[2, ]), col = col[2])




