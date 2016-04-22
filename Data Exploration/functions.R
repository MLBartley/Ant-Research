#####
#####
##### Functions useful for Trophallaxis Data Analysis
#####
#####

###########################################################

##
## INDEX
##

# 1 - sim.mmpp 
# 2 - mcmc.troph
# 3 - mcmc.troph.cov 

###########################################################
#purpose of function: Simulate data
    ## X_t ~ Markov Chain ( X_t-1 , P ) [two states]
    ## Y_t ~ Pois (lambda_{X_t})

#Inputs: tmax, delta.t, start.state, P, lambda

#Output: (1) - x: (1, 2) unobserved two state process 
#        (2) - y: (0, ...) "observed" number of interactions at time t
#        (3) - N: (0, ...) cumulative count of interactions
#        (4) - delta.t
#        (5) - t: (0, ...) time (in seconds)
#        (6) - 1x1 visual of cumlative counts over time, colored by state

sim.mmpp <- function(tmax,delta.t,start.state=1,P,lambda){
  T=floor(tmax/delta.t)+1
  x=rep(NA,T)
  y=rep(NA,T)
  x[1]=start.state
  y[1]=0
  for(t in 2:T){
    ## sample latent state
    x[t]=sample(1:2,1,prob=P[x[t-1],])
    ## sample observed events
    y[t]=rpois(1,lambda=lambda[x[t]]*delta.t)
  }
  list(y=y,x=x,N=cumsum(y),delta.t=delta.t,t=(0:(T-1))*delta.t)
  
  par(mfrow = c(1, 1))
  plot((0:(T-1))*delta.t, cumsum(y), 
       type = "p", pch = ".", cex = 2, col = x,
       xlab = "Time", ylab = "Interactions")
}




#purpose of function: to find MCMC generated estimates of 
# (1) - the state (X_t = high/low troph rates) of the colony at time t
# (2) - the specific rates of interaction (lambda) of each state 
# (3) - the probability of moving from one state to another (P matrix)

#Output: (1) - estimates of X, lambda, P
#        (2) - 2x2 visual of estimates over time (runs)


mcmc.troph = function(data, title, a, b, theta, states, n.mcmc){
  
  Time = length(data)
  n = states
  delta = rep(1/n, n)
  
  library(gtools)
  
  #homes
  ## Build Homes for X(1:T), lambda(1:n), and P(nXn) and gam vectors
  
  X.param = matrix(data = rep(NA, Time * n.mcmc), nrow = Time, 
                   ncol = n.mcmc, byrow = T)
  
  lambda.param = matrix(data = rep(NA, n * n.mcmc), nrow = n, 
                        ncol = n.mcmc, byrow = T)
  
  P.param = matrix(data = rep(NA, n * n * n.mcmc), nrow = n * n, 
                   ncol = n.mcmc, byrow = T)
  
  gam = matrix(NA, nrow = Time, ncol = n, 
                 byrow = T)
  

 ## Initialize parameters

  X.param[,1] = rep(1, Time)

  lambda.param[, 1] = rgamma(n = n, shape = a, rate = b)
  
  P.matrix = matrix(data = theta/100, nrow=n, ncol =n, byrow = T) 
  
  P.param[,1] = as.vector(t(P.matrix))
  #holds all P.parameter values over runs
  
  ## Gibbs Updates
  
  for(l in 2:n.mcmc) {
    
    m = matrix(data = rep(0, n * n), nrow = n, ncol = n) 
    # number states going from i to j, refreshes every run
    
    P.matrix = matrix(data = c(P.param[, l - 1]), nrow = n, ncol = n, byrow = T)
    
    
    
    ##X Parameters
    for(k in 1:n){
      gam[1, k] = lambda.param[k, l - 1] ^ data[1] * exp(-lambda.param[k, l - 1]) * 
        delta[k] * P.matrix[k, X.param[2, l - 1]]
    }
    
    X.param[1, l] = sample(x = (1:n), size = 1, prob = gam[1,])
    
    m[X.param[1, l], X.param[1, l]] = m[X.param[1, l], X.param[1, l]] + 1
    
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
    
    #Lambda and P parameters
    for(h in 1:n) {
      lambda.param[h, l] = rgamma(n = 1, shape = 
                                    sum(data[which(X.param[, l] == h)]) + a,
                                  rate = sum(m[h, ]) + b )
      
      P.matrix[h, ] = (rdirichlet(n = 1 , alpha = theta[h, ] + m[h, ]))
      
    }
    P.param[, l] = as.vector(t(P.matrix))
  }
  
  ## Compile the Estimates
  
  ## X1:XT, Lambda, Pmatrix
  
  #homes
  X.est = matrix(data = rep(NA, Time), nrow = Time, ncol = 1)
  lambda.est = matrix(data = rep(NA, n), nrow = n, ncol = 1)
  P.est = matrix(data = rep(NA, n * n), nrow = n * n, ncol = 1)
  
  
  for(t in 1:Time ){
    X.est[t, 1] = mean(X.param[t, ])  
  }
  
  for(i in 1:n ){
    lambda.est[i, 1] = mean(lambda.param[i, ])
  }  
  
  for(l in 1:(n * n) ){
    P.est[l, 1] = mean(P.param[l, ])
  }
  
  P.est.matrix = matrix(data = c(P.est[, 1]), nrow = n, ncol = n, byrow = T)
  
  #plot the estimation runs.
  
  
  #lambda
  par(mfrow = c(2,2),
      oma = c(0,0,2,0) + 1,
      mar = c(1,1,1,1) + 3)
    plot(0,0,xlab="MCMC Runs", ylab="Lambda", ylim=c(0,max(lambda.param)), 
       xlim=c(0,n.mcmc), type="n", cex.lab = 1)
  for(i in 1:n){
     lines(1:n.mcmc, lambda.param[i, ], col = i)
  }
 
  
  #P
  plot(0,0,xlab="MCMC Runs", ylab="P", ylim=c(0, max(P.param)), xlim=c(0,n.mcmc), 
       type="n", cex.lab = 1)
  for(i in 1:(n*n)){
    lines(1:n.mcmc, P.param[i, ], col = i)
  }
  
  #Single X
  X = X.param[sample(1:Time, 1), ]
  plot(0, 0, xlab="MCMC Runs", ylab="Single X", ylim=c(0,max(X)), 
       xlim=c(0,n.mcmc), type="n", cex.lab = 1)
  lines(1:n.mcmc, X, col = "red")
  
  #States over time
  plot(X.est,type="l",lwd=3, cex.lab = 1)
  
    title(main=title, outer=T)

  
  list(X.est = X.est, lambda.est = lambda.est, P.est = P.est.matrix, P.run = P.param)
  
}



#purpose of function: to find MCMC generated estimates of 
# (1) - the state (X_t = high/low troph rates) of the colony at time t
# (2) - the specific rates of interaction (lambda) of each state 
# (3) - the probability of moving from one state to another (P matrix)

#Output: (1) - estimates of X, lambda, P
#        (2) - 2x2 visual of estimates over time (runs)

#New aspect: now we will incorporate a multinomial probit model 
# that handles the covariate information and uses it to inform
#  the P matrix estimates over time. 

mcmc.troph.cov = function(data, covariates, title, a, b, theta, states, n.mcmc){
  
  Time = length(data)
  n = states
  delta = rep(1 / n, n)
  
  library(gtools)
  
  #homes
  ## Build Homes for X(1:T), lambda(1:n), and P(nXn) and gam vectors
  
  X.param = matrix(data = rep(NA, Time * n.mcmc), nrow = Time, 
                   ncol = n.mcmc, byrow = T)
  
  lambda.param = matrix(data = rep(NA, n * n.mcmc), nrow = n, 
                        ncol = n.mcmc, byrow = T)
  
  P.param = matrix(data = rep(NA, n * n * n.mcmc), nrow = n * n, 
                   ncol = n.mcmc, byrow = T)
  
  gam = matrix(NA, nrow = Time, ncol = n, 
               byrow = T)
  
  
  ## Initialize parameters
  
  X.param[,1] = rep(1, Time)
  
  lambda.param[, 1] = rgamma(n = n, shape = a, rate = b)
  
  P.matrix = matrix(data = theta/100, nrow=n, ncol =n, byrow = T) 
  
  P.param[,1] = as.vector(t(P.matrix))
  #holds all P.parameter values over runs
  
  ## Gibbs Updates
  
  for(l in 2:n.mcmc) {
    
    m = matrix(data = rep(0, n * n), nrow = n, ncol = n) 
    # number states going from i to j, refreshes every run
    
    P.matrix = matrix(data = c(P.param[, l - 1]), nrow = n, ncol = n, byrow = T)
    
    
    
    ## X Parameters
    for(k in 1:n){
      gam[1, k] = lambda.param[k, l - 1] ^ data[1] * exp(-lambda.param[k, l - 1]) * 
        delta[k] * P.matrix[k, X.param[2, l - 1]]
    }
    
    X.param[1, l] = sample(x = (1:n), size = 1, prob = gam[1,])
    
    m[X.param[1, l], X.param[1, l]] = m[X.param[1, l], X.param[1, l]] + 1
    
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
    
    ## Lambda and P parameters
    for(h in 1:n) {
      lambda.param[h, l] = rgamma(n = 1, shape = 
                                    sum(data[which(X.param[, l] == h)]) + a,
                                  rate = sum(m[h, ]) + b )
      
      P.matrix[h, ] = (rdirichlet(n = 1 , alpha = theta[h, ] + m[h, ]))
      
    }
    P.param[, l] = as.vector(t(P.matrix))
  }
  
  ## Compile the Estimates
  
  ## X1:XT, Lambda, Pmatrix
  
  #homes
  X.est = matrix(data = rep(NA, Time), nrow = Time, ncol = 1)
  lambda.est = matrix(data = rep(NA, n), nrow = n, ncol = 1)
  P.est = matrix(data = rep(NA, n * n), nrow = n * n, ncol = 1)
  
  
  for(t in 1:Time ){
    X.est[t, 1] = mean(X.param[t, ])  
  }
  
  for(i in 1:n ){
    lambda.est[i, 1] = mean(lambda.param[i, ])
  }  
  
  for(l in 1:(n * n) ){
    P.est[l, 1] = mean(P.param[l, ])
  }
  
  P.est.matrix = matrix(data = c(P.est[, 1]), nrow = n, ncol = n, byrow = T)
  
  #plot the estimation runs.
  
  
  #lambda
  par(mfrow = c(2,2),
      oma = c(0,0,2,0) + 1,
      mar = c(1,1,1,1) + 3)
  plot(0,0,xlab="MCMC Runs", ylab="Lambda", ylim=c(0,max(lambda.param)), 
       xlim=c(0,n.mcmc), type="n", cex.lab = 1)
  for(i in 1:n){
    lines(1:n.mcmc, lambda.param[i, ], col = i)
  }
  
  
  #P
  plot(0,0,xlab="MCMC Runs", ylab="P", ylim=c(0, max(P.param)), xlim=c(0,n.mcmc), 
       type="n", cex.lab = 1)
  for(i in 1:(n*n)){
    lines(1:n.mcmc, P.param[i, ], col = i)
  }
  
  #Single X
  X = X.param[sample(1:Time, 1), ]
  plot(0, 0, xlab="MCMC Runs", ylab="Single X", ylim=c(0,max(X)), 
       xlim=c(0,n.mcmc), type="n", cex.lab = 1)
  lines(1:n.mcmc, X, col = "red")
  
  #States over time
  plot(X.est,type="l",lwd=3, cex.lab = 1)
  
  title(main=title, outer=T)
  
  
  list(X.est = X.est, lambda.est = lambda.est, P.est = P.est.matrix, P.run = P.param)
  
}
