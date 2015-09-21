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

lamda.1 = 2
r.1 = 10

lamda.2 = 200
r.2 = 10

P = matrix(
  c(.99, .01, .01, .99),
  nrow = 2,
  ncol = 2,
  byrow = TRUE)
  
d = rep(NA, T) #build the vector of distances

##
## Step 2 - draw d_0 ~ Gamma(lamda.1, lamda.1)
## 

d[1] = rgamma(1, lamda.1, r.1) #set the first distance by drawing from a gamma dist
                               # with the parameters for X = 1 (our starting X)

##
## Step 3 - loop through changing states in X and then draw new d
##

for(t in 2:T){
  X[t] = sample(x = 1:2, size = 1, 
                replace = TRUE, prob = P[X[t-1], ]) #using sample function 
                                                    #to draw new state based on 
                                                    #previous state and P
  if(X[t] == 1){
  d[t] = rgamma(1, lamda.1, r.1 )
  }
  else {
  d[t] = rgamma(1, lamda.2, r.2)
  }
} 
#instead, make lambda and r vectors and use r[X(t)], etc
#this also help when extending to n number of states
##
## Step 4 - Visualize it!
##

plot(x = 1:T, y = d, col=X)
