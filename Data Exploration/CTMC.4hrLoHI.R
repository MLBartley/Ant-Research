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

prep.high = prep.troph.pairs(high4)
prep.low = prep.troph.pairs(low4)

N = prep.high$pairs2
N.low = prep.low$pairs2

Time = prep.high$hours * 60 * 60


###
# Simulated Data
###




###
# Inference
###

#Propose high and low gamma and lambda, 
  #first simplest case
  # then with additional state based conditions
#calculate R* and then P* matrices
#accept/reject

#gibbs updates on X_t and M

n.mcmc = 3000
a = 0.03 #why is this changing to 4
b = 0.03
r = 0.03
q = 0.03 
theta = matrix(c(99, 1, 1, 99), 2, 2) #prior on state transition probabilities

#might need new a and b for gamma_h^tilda 
tau = c(0.05, 0.05, 0.05, 0.05) #tuning



data = N.low[1:(2*60*60)] #only using first two hours for time


#homes

params = matrix(NA, nrow = 4, ncol = n.mcmc)
rownames(params) <- c("gamma_high", "gamma_low", "lambda_high_tilda", "lambda_low")
colnames(params) <- 1:n.mcmc

X.param = matrix(NA, nrow = Time, 
                 ncol = n.mcmc, byrow = T)


# Ne probability matrix for 2 state discrete time Markov Chain
M.param = matrix(data = NA, nrow = 4, 
                 ncol = n.mcmc, byrow = T)
  rownames(M.param) <- c("HH", "HL", "LH", "LL")
  colnames(M.param) <- 1:n.mcmc
M = matrix(c(0.9999, 0.0001, 0.0001, 0.9999),byrow = T,  2, 2)
  rownames(M) <- c("High", "Low")
  colnames(M) <- c("High", "Low")


#initialize

params[1:2, 1] = rgamma(n = 2, shape = a, rate = b)
params[3:4, 1] = rgamma(n = 2, shape = r, rate = q)

X.param[, 1] = rep(1, Time) #1 - high state, 2 - low state

M.param[, 1] = as.vector(t(M))


#posterior likelihood

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


  loglike = sumP_L + sumP_H +
                   dgamma(params[1], a, b, log = T) + 
                   dgamma(params[2], a, b, log = T) + 
                   dgamma(params[3], r, q, log = T) +
                   dgamma(params[4], r, q, log = T)

  return(loglike)
}


accept = 0

#MCMC

for(l in 2:n.mcmc){

  if( l %% 10 == 0 ) cat(paste("iteration", l, "complete\n")) 
  #
  #update
  #
  
  if(l < n.mcmc/2 & l %% 100 == 0){
    
    sigma = c(0, 0, 0, 0)
    for(a in 1:4){
      sigma[a] = ((2.38 ^ 2) / 3) * var(params[a, 1:(l - 1)])
    }
    
    if((sigma[1] != 0) & (sigma[2] != 0) & (sigma[3] != 0) & sigma[4] != 0){
      tau = sigma
    }

  }
  
  proposal = rnorm(4, mean = log(params[, l - 1]), sd = tau)
  
  theta.star = exp(proposal)
  
  theta.star[1] = theta.star[1] + theta.star[2] #calculate correct gamma_high
 
  
 #calculate P matrices - for high/low states
  
  R_H = matrix(0, nrow = max(data) + 1, ncol = max(data) + 1)
  rownames(R_H) <- 0:max(data)
  colnames(R_H) <- 0:max(data)
  
  for(i in 1:nrow(R_H)){
    
    if( i %% 2 != 0 & i != nrow(R_H) & i != (nrow(R_H) - 1)){
      R_H[i, i + 2] = theta.star[1]  #gamma_high
    }
    
    if( i %% 2 != 0 & i != 1 & i != 2){
      R_H[i, i - 2] = (i - 1) * theta.star[3] #lambda_high 
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
      R_L[i, i + 2] = theta.star[2]  #gamma_low
    }
    
    if( i %% 2 != 0 & i != 1 & i != 2){
      R_L[i, i - 2] = (i - 1) * theta.star[4] #lambda_low 
    }
    
  }
  
  P_L = Pctmc(Q = R_L, t = 1)
  rownames(P_L) <- 0:max(data)
  colnames(P_L) <- 0:max(data)
  
  
  #calculate probability
  MHprob = exp(log.fullcond(theta.star, P_L, P_H, data, X.param) - log.fullcond(params[, l - 1], P_L, P_H, data, X.param))
 
  #accept/reject

  
  if(runif(1) < MHprob){
    accept = accept + 1
    params[, l] = theta.star
  }else{
    params[, l] = params[, l - 1]
  }
  
   
  ##X Parameters
  
  m = matrix(data = 0, nrow = 2, ncol = 2) 
  # number states going from i to j, refreshes every run
  
  
  prob_H = P_H[data[1] + 1, data[2] + 1] * .5 *
    M[X.param[1, l - 1], X.param[2, l - 1]] 
    
  
  prob_L = P_L[data[1] + 1, data[2] + 1] * .5 *
   M[X.param[1, l - 1], X.param[2, l - 1]]
  
  X.param[1, l] = sample(x = (1:2), size = 1, prob = c(prob_H, prob_L))
  
  m[X.param[1, l], X.param[1, l]] = m[X.param[1, l], X.param[1, l]] + 1
  

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

  #M matrix parameter
    
  M[1, ] = (rdirichlet(n = 1 , alpha = theta[1, ] + m[1, ]))
  M[2, ] = (rdirichlet(n = 1 , alpha = theta[2, ] + m[2, ]))
 
  M.param[, l] = as.vector(t(M))
   
}



  #compile estimates
X.est = matrix(data = rep(NA, Time), nrow = Time, ncol = 1)


gamma.high.est = mean(params[1, ])
gamma.low.est = mean(params[2, ])
lambda.high.est = mean(params[3, ])
lambda.low.est = lambda.high.est + mean(params[4, ])

for(t in 1:Time ){
  X.est[t, 1] = mean(X.param[t, ])  
}

  #plot the estimation runs.

col = c("#120d08", "#bc5356", "#538bbc", "#53bc84")

#gamma
plot(0,0,xlab="MCMC Runs",
     ylab="Gamma (per minute)",
     ylim=c(0,max(params[1:2, ]) * 60), 
     xlim=c(0,n.mcmc), 
     type="n",
     cex.lab = 1)
  lines(1:n.mcmc, 60 * params[1, ], col = col[1])
  lines(1:n.mcmc, 60 * params[2, ], col = col[2])


#lambda
plot(0,0,xlab = "MCMC Runs", ylab = "Lambda (per minute)", 
     ylim = c(0, max(params) * 60), xlim = c(0, n.mcmc), 
     type = "n", cex.lab = 1)

  lines(1:n.mcmc, (60 * params[3, ]), col = col[3])
  lines(1:n.mcmc, (60 * params[4, ]), col = col[4])

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

  
  #########################################################
  ##
  ## Fancy Plots with Background Colors
  ##
  #########################################################
  start = low4$start_time
  location = low4$Location
  
  par(mfrow = c(1, 1))
  
  if(length(unique(location)) == 1){
    
    ##High Density - 4 Hours
    plot(start, 1:int.num, main="High", 
         xlab="Seconds", ylab = "Cumulative Interaction Count", 
         xlim = c(0, maxtime))
    states = X.est #from code above
    rr = rle(states[,1])
    rr$values = round(rr$values, digits = 0)
    embedded.chain = rr$values
    cs = c(0,cumsum(rr$lengths))*delta.t - delta.t
    cols=c('#bc535644','#538bbc44')
    for(j in 1:length(embedded.chain)){
      rect(cs[j],0,cs[j + 1],int.num, 
           col = cols[embedded.chain[j]], density = NA)
    }
  }
  
  
  #else{
    
    start = sort(low4$start_time)
    location = low4$Location
    int.num = length(start)
    maxtime = Time
    delta.t = 1 
    
    #Low Density - 4 Hours
    
    plot(start, 1:int.num, main="Low", xlab="Seconds",
         ylab = "Cumulative Interaction Count", 
         xlim=c(0,maxtime))
    states = X.est
    rr=rle(states[,1])
    rr$values = round(rr$values, digits = 0)
    embedded.chain=rr$values
    cs=c(0,cumsum(rr$lengths))*delta.t - delta.t
    cols=c('#bc535644','#538bbc44')
    for(j in 1:length(embedded.chain)){
      rect(cs[j],0,cs[j+1],int.num, col=cols[embedded.chain[j]] , density=NA)
    }
    