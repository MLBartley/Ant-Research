##
## MCMC Algorithm - Trophallaxis Model 
## 
######################################

## Information Needed

#data 

#data = ylist[[1]]

sim.data = sim$y
data = sim.data

library("gtools")
  
Time = length(data)

n = 2 #number of proposed states

n.mcmc = 1000 #how many mcmc runs to do


##hyperparameters
a = 11 
b = 2
theta = matrix(data = c(9, 1, 1, 9), nrow = 2, ncol = 2, byrow = T) 
##would need two priors if want to weight the diagonals more heavily


delta = c(.5, .5) 

## Build Homes for X(1:T), lambda(1:n), and P(nXn) and gam vectors

X.param = matrix(data = rep(NA, Time * n.mcmc), nrow = Time, 
                 ncol = n.mcmc, byrow = T)

lambda.param = matrix(data = rep(NA, n * n.mcmc), nrow = n, 
                      ncol = n.mcmc, byrow = T)

P.param = matrix(data = rep(NA, n * n * n.mcmc), nrow = n * n, 
                 ncol = n.mcmc, byrow = T)

gam.1 = matrix(data = rep(NA, n), nrow = 1, ncol = n, 
               byrow = T)
  
gam.t = matrix(data = rep(NA, n * (Time - 2)), nrow = Time - 2, 
               ncol = n, byrow = T)
  
gam.T = matrix(data = rep(NA, n), nrow = 1, ncol = n, byrow = T)


## Initialize parameters


  X.param[,1] = rep(1, Time)
#X.param[,1:n.mcmc] = sim$x

lambda.param[, 1] = rgamma(n = 2, shape = a, rate = b)


P.param[,1] = c(.99, .01, .01, .99) #stay in same states
 
P.matrix = matrix(data = c(P.param[, 1]), nrow=n, ncol =n, byrow = T) 
#holds all P.parameter values over runs


## Gibbs Updates

for(l in 2:n.mcmc) {
  
  m = matrix(data = rep(0, n * n), nrow = n, ncol = n) 
  # number states going from i to j, refreshes every run

  P.matrix = matrix(data = c(P.param[, l - 1]), nrow = n, ncol = n, byrow = T)
 
  
  
  ##X Parameters
   for(k in 1:n){
    gam.1[1, k] = lambda.param[k, l - 1] ^ data[1] * exp(-lambda.param[k, l - 1]) * 
      delta[k] * P.matrix[k, X.param[2, l - 1]]
   }

   X.param[1, l] = sample(x = (1:n), size = 1, prob = gam.1)
  
          m[X.param[1, l], X.param[1, l]] = m[X.param[1, l], X.param[1, l]] + 1
   
  for(t in 2:(Time - 1)){
    
    for(k in 1:n){
      gam.t[t - 1, k] = (lambda.param[k, l - 1] ^ data[t]) * 
        exp(-lambda.param[k, l - 1]) * P.matrix[X.param[t - 1, l - 1], k] *
        P.matrix[k, X.param[t + 1, l - 1]]
    }
    
     X.param[t, l] = sample(x = (1:n), 1,  prob = gam.t[t - 1, ]) 
     
        m[X.param[t - 1, l], X.param[t, l]] = m[X.param[t - 1, l], 
                                                X.param[t,l]] + 1
  }
  
  for(k in 1:n){
      gam.T[k] = lambda.param[k, l - 1] ^ data[Time] * 
        exp(-lambda.param[k, l - 1]) * P.matrix[X.param[(Time - 1), l - 1], k]
  }
  
 X.param[Time, l] = sample(x = 1:n, 1,  prob = gam.T[1, ])
      
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

# X.est
# lambda.est
# P.est.matrix


#plot the estimation runs.


#lambda
par(mfrow = c(1, 1))
plot(1:n.mcmc, lambda.param[2, ])
lines(1:n.mcmc, lambda.param[1, ], col = "red")

##
## examine state recovery
##


true.states=sim$x

estimated.states=X.est
plot(true.states,type="l",lwd=3)
points(estimated.states,type="l",col="red",lwd=2,lty=2)




}

########################
##
## And now with a little umph:
## incorperating the ant data
##
########################
source("functions.R")

mcmc.troph(data = data, a = a, b = b, theta = theta, states = 2, n.mcmc = 1000)

data = ylist[[8]]



