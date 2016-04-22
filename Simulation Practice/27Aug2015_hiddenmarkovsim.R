##
##
##  Simulation for week of 27 August 2015
##  Hidden Markov Chains - Animal movement
##
########################################

##
## Step 1 - set initial values
##
T = 1000 #this will be the number in our vectors and the number of iterations (minus 1) to run
X = rep(NA, T) #build the vector to start to save on computing effort later
X[1] = 1           #this will be our two step markov chain, state one relates to 
                   # short distances (small mean), and state two the opposite

# lambda.1 = 2
# r.1 = 10
# 
# lambda.2 = 200
# r.2 = 10

lambda = c(2,20)
r = c(10,10)

P = matrix(
  c(.99, .01, .01, .99),
  nrow = 2,
  ncol = 2,
  byrow = TRUE)
  
d = rep(NA, T) #build the vector of distances

##
## Step 2 - draw d_0 ~ Gamma(lambda.1, lambda.1)
## 

d[1] = rgamma(1, lambda[1], r[1]) #set the first distance by drawing from a gamma dist
                               # with the parameters for X = 1 (our starting X)

##
## Step 3 - loop through changing states in X and then draw new d
##

for(t in 2:T){
  X[t] = sample(x = 1:2, size = 1, 
                replace = TRUE, prob = P[X[t-1], ]) #using sample function 
                                                    #to draw new state based on 
                                                    #previous state and P
  d[t] = rgamma(1, lambda[X[t]], r[X[t]] )
} 
#instead, make lambda and r vectors and use r[X(t)], etc
#this also help when extending to n number of states

##
## Step 4 - Visualize it!
##

plot(x = 1:T, y = d, col=X, type = "l")

##
## Step 5 - likelihood and forward equations
##

#use simulated distances (d) as my x here.
x = d

#T= length(x)
likelihood = (NA)
ones = rep(1, 2)
delta = c(.5,.5) #based on Probability matrix above

neg.log.like = function(theta, x){
  lambda.1 = theta[1]
  lambda.2 = theta[2]
  r.1 = theta[3]
  r.2 = theta[4]
  
  ones = rep(1, 2)
  
  L.mat = matrix(
    c(dgamma(x, lambda.1, r.1), dgamma(x, lambda.2, r.2)),
    nrow = length(x),
    ncol = 2
  ) 
   
  alpha.vec = matrix(
    c(rep(NA, length(x)*2)),
    nrow = length(x),
    ncol = 2
  )
  
  P_t = diag(L.mat[1,]) 
  alpha.vec[1,] = delta %*% P_t
  
  for(i in 2:T){
    P_t = diag(L.mat[i,]) 
    alpha.vec[i,] = alpha.vec[i-1,] %*% P %*%  P_t
  }
  like = (alpha.vec[T,] %*% (ones))
  return(-log(like))
} 
theta = c(2, 200, 10, 10)
neg.log.like(theta, x)

#then get neg log likelihood so we can minimize with optim()
# neg.l.like = function(theta, x=d){
#   -log(likelihood(theta, x))
#   }
# neg.l.like(theta, x)

#use optim() L-BFGS-B
optim(theta <- c(1,100,5,8),lower = c(0, 0, 0, 0),
      fn = neg.log.like, x = d, method = "L-BFGS-B")

#check out scaling
#Main idea: before the product of probabilities
#were quickly going to 0 and then trouble happens 
#when we try and do log(0)

#Instead: Algarithm with computational trick of multiplying
#the values by another value (then later canceling it out)
#There is an example in HMM Appendex A.1.3

#We need to redo our likelihood functions in this mannor

scale.neg.log.likelihood = function(theta, x=d){
  lambda.1 = theta[1]
  lambda.2 = theta[2]
  r.1 = theta[3]
  r.2 = theta[4] 
  
  phi = matrix(rep(NA, T*2),
               nrow = T,
               ncol =2)
  
  phi[1,] = delta
  
  lscaled = 0
  
  L.mat = matrix(
    c(dgamma(x, lambda.1, r.1), dgamma(x, lambda.2, r.2)),
    nrow = length(x),
    ncol = 2
  ) 

  ones = rep(1, 2)
  
  for(t in 2:T){
    P_t = diag(L.mat[t,]) 
    v.vec = phi[t-1,] %*% P %*% P_t
    u = v.vec %*% ones
    lscaled = lscaled + log(u)
    phi[t,] = v.vec/as.numeric(u)
    return(-lscaled)
  }
}


theta = c(2, 200, 10, 10)
neg.l.like(theta, x = d)
scale.neg.log.likelihood(theta = theta, x = d)

#should this give the same value as neg.log.likelihood

optim(theta <- c(5, 100, 10, 10), fn = scale.neg.log.likelihood, 
      method = "L-BFGS-B", lower = c(1,1, 1,1), 
      upper = c(400, 400, 999,400))

save(d, file="d.Rdata")
