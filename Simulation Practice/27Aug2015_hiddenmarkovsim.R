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

lambda = c(2,200)
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

plot(x = 1:T, y = d, col=X)

##
## Step 5 - likelihood and forward equations
##

#use simulated distances (d) as my x here.
x = d
T= length(x)
likelihood = (NA)
ones = rep(1, 2)
delta = c(.5,.5) #based on Probability matrix above

#still need a matrix akin to P(x) (different from ptm) as in the text
# I believe in this 2 state case, P(x) will be a 2x2 matrix?

# pgamma(1, lambda[1], r[1])
# pgamma(1, lambda[2], r[2])

L.mat = matrix(
  c(dgamma(x, lambda[1], r[1]), dgamma(x, lambda[2], r[2])),
  nrow = length(x),
  ncol = 2
)

# cond.prob = function(x,c){
#   dgamma(x, lambda[c], r[c]) #don't for get ddist is pdf and pdist is cdf
# }
# 
# cond.prob.matrix = function(x){
#   cp.matrix = matrix(
#     c(cond.prob(x,1), 0, 0, cond.prob(x,2)),
#     nrow = 2,
#     ncol = 2
#   )
#     return(cp.matrix)
# }


alpha.vec = matrix(
    c(rep(NA, length(x)*2)),
      nrow = length(x),
      ncol = 2
      )

  for(i in 2:T){
  P_t = diag(L.mat[i,]) 
  alpha.vec[1,] = delta %*% P_t #i think this should be a 1x2 matrix? but does that make sense?
  alpha.vec[i,] = alpha.vec[i-1,] %*% P %*%  P_t
}


likelihood = function(x, theta){
  lambda.1 = theta[1]
  lambda.2 = theta[2]
  r.1 = theta[3]
  r.2 = theta[4]
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
  
  for(i in 2:T){
    P_t = diag(L.mat[i,]) 
    alpha.vec[1,] = delta %*% P_t #i think this should be a 1x2 matrix? but does that make sense?
    alpha.vec[i,] = alpha.vec[i-1,] %*% P %*%  P_t
  }
  return(alpha.vec[T,] %*% (ones))
} 
theta = c(2, 200, 10, 10)
likelihood(x, theta)

#then get neg log likelihood so we can minimize with optim()
neg.l.like = function(x, theta){
  -log(likelihood(x, theta))
  }
neg.l.like(x, theta)

#use optim() L-BFGS-B
optim(theta <- c(1, 100, 5, 5), fn = neg.l.like, x=d, 
          method = "BFGS")

#check out scaling
neg.l.like(x,theta)
