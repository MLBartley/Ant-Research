#######
##
## 7 APRIL 2016 
## 
## Want to do a simulation of my model for ant interactions in a simpler
## manner so that I may incoperate auxilary variables for coefficients. 
##
##
########################################################################

##  Model needs to have the following
      # X_i - state of interactions. Currently high/low, but want to keep general
      # lambda - state based rate of interactions
      # P matrix - probability of switching between states
      
  ## And Then...
      # w_i - coefficients
      # beta 
      # z_i - latent auxillary variable
      # alpha_i 

########################################################################


## First Step - need to model the process with fixed parameters so we may
##               simulate data

#########################################################
##
## Simulation Function
##
## X_t ~ Markov Chain ( X_t-1 , P )
## Y_t ~ Pois (lambda_{X_t})
#########################################################

##
## simulate from the model
##

P = matrix(c(.99, .01, .01, .99), nrow = 2, byrow = T)
lambda = k = c(1, 10)

delta.t = 1 #needs to be 1, else observations dependent on time

sim = sim.mmpp(1000, delta.t, start.state = 1, P, lambda)
max(sim$y)



## Second Step - Load the covariates and other inputs
##

######## covariates 

#Foraging Movement Data
cov.data = read.csv("~/Google Drive/PSU/Projects/Ant-Research/Data/Colony1_high_foraging_2hr.csv")
               
foragers = unique(cov.data$Ant_ID)

leave = NA

for(i in foragers){
  leave = cbind(leave, max(cov.data$end[which(cov.data$Ant_ID == i)]))
}

leave = leave[-1]
leave = c(leave, 3000, 5000)
leave = sort(leave)

cov = c(rep(1000, leave[1] - 1), 
        0:(leave[2] - leave[1] - 1), 
        0:(leave[3] - leave[2] - 1), 
        0:(leave[4] - leave[3] - 1), 
        0:(leave[5] - leave[4] - 1), 
        0:(leave[6] - leave[5] - 1), 
        0:(7200 - leave[6])) 


#In/Out Chamber 4 Movement Data

in.out =  read.csv("Data/Colony_1_in&out_high_density_4hrs.csv")
in.out = in.out[which(in.out$Action == "enter"),]
in.out = in.out[order(in.out$time), ]


#Create vector of covariates - time since last ant entered
cov = rep(300, in.out$time[1])
  
  for(i in 2:nrow(in.out)){
    cov = c(cov, 0:(in.out$time[i] - in.out$time[i - 1] ))
  }

cov = c(cov, 0:(14400 - in.out$time[nrow(in.out)]))

# inputs
data = sim$y
states = 2
n.mcmc = 1000
theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 
a = 5
b = 2


mu.all = c(2, -1, -0.000004)
sig.all = matrix(data = c(0.2, 0, 0, 
                          0, 0.2, 0, 
                          0, 0, 0.0002), nrow = 3, ncol = 3, byrow = T)
#mu.a = 2
#mu.b = c(-1, -0.004)
#tau.a = 0.2
#tau.b = matrix(data = c(0.2, 0, 0, 0.2), nrow = 2, ncol = 2, byrow = T)

## Third Step - need to estimate X, lambda, P
## OPTION 1 - Latent Auxillary Variable

Time = length(data)
n = states
delta = rep(1 / n, n)

#library(gtools)

#homes
## Build Homes for X(1:T), lambda(1:n), and P(nXn) and gam vectors

X.param = matrix(data = rep(NA, Time * n.mcmc), nrow = Time, 
                 ncol = n.mcmc, byrow = T)

lambda.param = matrix(data = rep(NA, n * n.mcmc), nrow = n, 
                      ncol = n.mcmc, byrow = T)

# P.param = matrix(data = rep(NA, n * n * n.mcmc), nrow = n * n, 
#                  ncol = n.mcmc, byrow = T)

## P matrix now varies over time, need new homes

P.11.param = matrix(NA, nrow = Time, ncol = n.mcmc)
P.12.param = matrix(NA, nrow = Time, ncol = n.mcmc)
P.21.param = matrix(NA, nrow = Time, ncol = n.mcmc)
P.22.param = matrix(NA, nrow = Time, ncol = n.mcmc)

gam = matrix(NA, nrow = Time, ncol = n)

alph.beta.params = matrix(NA, nrow = 3, ncol = n.mcmc)

# alpha.param = matrix(NA, nrow = 1, ncol = n.mcmc)
# betas.param = matrix(NA, nrow = 2, ncol = n.mcmc)


## Initialize parameters

X.param[, 1] = rep(1, Time)

lambda.param[, 1] = rgamma(n = n, shape = a, rate = b)

P.matrix = matrix(data = theta / 100, nrow=n, ncol =n, byrow = T) 

#P.param[, 1] = as.vector(t(P.matrix)) #holds all P.parameter values over runs

P.11.param[, 1] = P.matrix[1, 1]
P.12.param[, 1] = P.matrix[1, 2]
P.21.param[, 1] = P.matrix[2, 1]
P.22.param[, 1] = P.matrix[2, 2]

alph.beta.params[, 1] = rnorm(3, mean = mu.all, sd = sqrt(diag(sig.all)))
 
# alpha.param[1] = rnorm(1, mean = mu.a, sd = tau.a)
# 
# betas.param[ , 1] = rmvnorm(1, mean = mu.b, sigma = tau.b)
#   

## Metropolis Hastings/Gibbs Updates

for(l in 2:n.mcmc) {
  
  m = matrix(data = rep(0, n * n), nrow = n, ncol = n) 
  # number states going from i to j, refreshes every run
  

  
  P.matrix = matrix(data = c(P.11.param[1, l - 1], P.12.param[1, l - 1], 
                             P.21.param[1, l - 1], P.22.param[1, l - 1]), 
                            nrow = n, ncol = n, byrow = T)
  
  
          #P matrix w/ alpha, betas - MH Algorithm
          
            log.fullcond = function(param){
              
              P.mult = 0
                for(i in 2:Time){
                P.mult = P.mult + log(P.matrix[X.param[i - 1, l - 1], X.param[i, l - 1]])
              }
              out = P.mult + 
                dnorm(param[1], mean = mu.all[1], sd = sqrt(diag(sig.all)[1]), log = T) + 
                dnorm(param[2], mean = mu.all[2], sd = sqrt(diag(sig.all)[2]), log = T) +
                dnorm(param[3], mean = mu.all[3], sd = sqrt(diag(sig.all)[3]), log = T)
              
              return(out)
            }
  
  
#           log.fullcond.alpha = function(param){
#             return(m[2,1] * (param - log(1 + exp(param))) - 
#                      (1 / (2 * tau.a ^ 2)) * (param - mu.a) ^ 2 )
#           }
#           
#           log.fullcond.betas = function(param, data.w){
#             return(m[1,2] * (t(param) %*% c(1, data.w) - 
#                                log(1 + exp(t(param) %*% c(1, data.w)))) + (dmvnorm(param, mean = mu.b, sigma = tau.b, log = T))
#             )
#           }
          
          #proposal
          
            proposal = rnorm(3, mean = alph.beta.params[, l - 1], sd =  c(.2, .2, .002))
            
#           proposal.alpha = rnorm(1, mean = alpha.param[l - 1], sd = tau.a)
#           
#           proposal.betas = rnorm(2, mean = betas.param[, l - 1], sd = diag(tau.b))
#           

            #accept/reject all params - Block update

            prob.all = exp(log.fullcond(proposal) - log.fullcond(alph.beta.params[, l - 1]))

            if(runif(1) < prob.all){
              alph.beta.params[, l] = proposal
            }else{
                alph.beta.params[, l] = alph.beta.params[, l-1]
              }
            
#           #accept/reject alpha 
#           prob.a = exp(log.fullcond.alpha(proposal.alpha) - log.fullcond.alpha(alpha.param[l-1]))
#           
#           if (runif(1) < prob.a){
#             alpha.param[l] = proposal.alpha
#           }
#           else{
#             alpha.param[l] = alpha.param[l-1]
#           }
#           
#           #accept/reject beta
#           prob.b = exp(log.fullcond.betas(proposal.betas, cov[l]) - 
#                          log.fullcond.betas(betas.param[, l-1], cov[l]))
#           
#           if (runif(1) < prob.b){
#             betas.param[, l] = proposal.betas
#           }else{
#             betas.param[, l] = betas.param[, l - 1]
#           }
          
          #Now use priors to calculate P matrix values 
          alpha = alph.beta.params[1, l]
          beta.0 = alph.beta.params[2, l]
          beta.1 = 0 
          #beta.1 = alph.beta.params[2, l]
          
          
          for(i in 1:Time){
            
          P.matrix[1, 2] = (exp(beta.0 + beta.1 * cov[i])) / (1 + exp(beta.0 + beta.1 * cov[i]))
          P.matrix[1, 1] = 1 - P.matrix[1, 2]
          P.matrix[2, 2] = (exp(alpha)) / (1 + exp(alpha))
          P.matrix[2, 1] = 1 - P.matrix[2, 2]  
          
          P.12.param[i, l] = P.matrix[1, 2]
          P.11.param[i, l] = P.matrix[1, 1]
          P.22.param[i, l] = P.matrix[2, 2]
          P.21.param[i, l] = P.matrix[2, 1]
          
          }
        
  ##X Parameters
     
    ## For X at t = 1
          
P.matrix = matrix(data = c(P.11.param[1, l - 1], P.12.param[1, l - 1], 
                          P.21.param[1, l - 1], P.22.param[1, l - 1]), 
                          nrow = n, ncol = n, byrow = T)        
  for(k in 1:n){
    gam[1, k] = lambda.param[k, l - 1] ^ data[1] * exp(-lambda.param[k, l - 1]) * 
      delta[k] * P.matrix[k, X.param[2, l - 1]]
  }
  
  X.param[1, l] = sample(x = (1:n), size = 1, prob = gam[1,])
  
  m[X.param[1, l], X.param[1, l]] = m[X.param[1, l], X.param[1, l]] + 1
  
  
     ##For X from t = 2 to t = Time - 1
  
  for(t in 2:(Time - 1)){
    
    P.matrix = matrix(data = c(P.11.param[t, l - 1], P.12.param[t, l - 1], 
                               P.21.param[t, l - 1], P.22.param[t, l - 1]), 
                                nrow = n, ncol = n, byrow = T)
    
    for(k in 1:n){
      gam[t, k] = (lambda.param[k, l - 1] ^ data[t]) * 
        exp(-lambda.param[k, l - 1]) * P.matrix[X.param[t - 1, l - 1], k] *
        P.matrix[k, X.param[t + 1, l - 1]]
    }
    
    X.param[t, l] = sample(x = (1:n), 1,  prob = gam[t, ]) 
    
    m[X.param[t - 1, l], X.param[t, l]] = m[X.param[t - 1, l], 
                                            X.param[t,l]] + 1
  }
  
  ## For X at t = Time
  
  P.matrix = matrix(data = c(P.11.param[Time, l - 1], P.12.param[Time, l - 1], 
                             P.21.param[Time, l - 1], P.22.param[Time, l - 1]), 
                    nrow = n, ncol = n, byrow = T)
  
  for(k in 1:n){
    gam[Time, k] = lambda.param[k, l - 1] ^ data[Time] * 
      exp(-lambda.param[k, l - 1]) * P.matrix[X.param[(Time - 1), l - 1], k]
  }
  
  X.param[Time, l] = sample(x = 1:n, 1,  prob = gam[Time, ])
  
  m[X.param[Time - 1, l], X.param[Time, l]] = m[X.param[Time - 1, l], 
                                                X.param[Time, l]] + 1
  
  
  #Lambda 
  for(h in 1:n) {
    lambda.param[h, l] = rgamma(n = 1, shape = 
                                  sum(data[which(X.param[, l] == h)]) + a,
                                rate = sum(m[h, ]) + b )
    
    # P.matrix[h, ] = (rdirichlet(n = 1 , alpha = theta[h, ] + m[h, ]))

  }

  
 
}

## Compile the Estimates

## X1:XT, Lambda, Pmatrix

#homes
X.est = matrix(NA, nrow = Time, ncol = 1)
lambda.est = matrix(data = rep(NA, n), nrow = n, ncol = 1)
#P.est = matrix(data = rep(NA, n * n), nrow = n * n, ncol = 1)

P.11.est = matrix(NA, nrow = Time, ncol = 1)
P.12.est = matrix(NA, nrow = Time, ncol = 1)
P.21.est = matrix(NA, nrow = Time, ncol = 1)
P.22.est = matrix(NA, nrow = Time, ncol = 1)

for(t in 1:Time ){
  X.est[t, 1] = mean(X.param[t, ])  
}

for(i in 1:n ){
  lambda.est[i, 1] = mean(lambda.param[i, ])
}  

for(t in 1:Time ){
  #P.est[l, 1] = mean(P.param[l, ])
  P.11.est[t, 1] = mean(P.11.param[t, ]) 
  P.12.est[t, 1] = mean(P.12.param[t, ]) 
  P.21.est[t, 1] = mean(P.21.param[t, ]) 
  P.22.est[t, 1] = mean(P.22.param[t, ])
  
}


#P.est.matrix = matrix(data = c(P.est[, 1]), nrow = n, ncol = n, byrow = T)
P.est.matrix =  matrix(data = c(mean(P.11.est), mean(P.12.est), 
                                mean(P.21.est), mean(P.22.est), 
                       nrow = n, ncol = n, byrow = T)


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
t = sample(1:Time, 1)
plot(0, 0, xlab="MCMC Runs", ylab="Single P", ylim=c(0, 1), xlim=c(0,n.mcmc), 
     type="n", cex.lab = 1)
#for(i in 1:(n * n)){
  lines(1:n.mcmc, P.11.param[t, ], col = "red")
  lines(1:n.mcmc, P.12.param[t, ], col = "blue")
  lines(1:n.mcmc, P.21.param[t, ], col = "black")
  lines(1:n.mcmc, P.22.param[t, ], col = "green")
#}

#Single X
X = X.param[sample(1:Time, 1), ]
plot(0, 0, xlab="MCMC Runs", ylab="Single X", ylim=c(0,max(X)), 
     xlim=c(0,n.mcmc), type="n", cex.lab = 1)
lines(1:n.mcmc, X, col = "red")

#States over time
plot(X.est,type="l",lwd=3, cex.lab = 1)

title(main=title, outer=T)

  ## Second Step - need to estimate X, lambda, P
## OPTION 2 - New P equations with alpha/beta priors















