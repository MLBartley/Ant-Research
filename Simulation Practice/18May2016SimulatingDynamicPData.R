#######
##
## 18 May 2016 
## 
## Want to do a simulation of data (both high and low density) where 
## the probability of engaging in trophallaxis increases when a Forager
## enters the nest. 
##
## Issue #1 on GitHub
########################################################################

#Outline

# estimate lambda
# Arrivals
# time since arrivals
# P matrix based on time since last arrivals
# use dynamic P matrix to simulate interactions

########
##
## Obtain Lambda estimate for simulations
##
########

cov.data = read.csv("./Data/Colony1_high_foraging_2hr.csv")
lambda = nrow(cov.data)/max(cov.data$start)
lambda

########
##
## Forager arrival times
##
########
time = 4 * 60 * 60
forager.number = rpois(n = 1, lambda = lambda*time)
forager.arrivals = rep(NA, forager.number)
  for (i in 1:forager.number){
    forager.arrivals[i] = rgamma(n = 1, shape = i,
                                 rate = lambda)
  }
  
forager.arrivals = round(forager.arrivals)
forager.arrivals = sort(forager.arrivals)
forager.arrivals = forager.arrivals[which(duplicated(forager.arrivals) == FALSE)]
forager.arrivals = c(forager.arrivals, time+1)
########
##
## Time since Arrivals
##
########

covariate = rep(NA, time)
covariate[1] = 300 #five minutes since arrival

arrival = 1

for(i in 2:time){
  
  if(forager.arrivals[arrival] == i){
    covariate[i] = 0
    arrival = arrival + 1
  }else{
    covariate[i] = covariate[i - 1] + 1
  }
  
}

########
##
## Calculate P matrix based on arrival data (covariate)
##
########
Time = time
theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 
n.mcmc = 1
n = 2

alpha = -2.2
beta.0 = -2.2
beta.1 = 0.0004

## P matrix now varies over time, needs homes

P.11.param = matrix(NA, nrow = Time, ncol = n.mcmc)
P.12.param = matrix(NA, nrow = Time, ncol = n.mcmc)
P.21.param = matrix(NA, nrow = Time, ncol = n.mcmc)
P.22.param = matrix(NA, nrow = Time, ncol = n.mcmc)

P.matrix = matrix(data = theta / 100, nrow = n, ncol = n, byrow = T) 

P.11.param[1, 1] = P.matrix[1, 1]
P.12.param[1, 1] = P.matrix[1, 2]
P.21.param[1, 1] = P.matrix[2, 1]
P.22.param[1, 1] = P.matrix[2, 2]

for(i in 2:Time){ 
  
  P.matrix[1, 2] = (exp(beta.0 + beta.1 * covariate[i])) / 
    (1 + exp(beta.0 + beta.1 * covariate[i]))
  P.matrix[1, 1] = 1 - P.matrix[1, 2]
  P.matrix[2, 1] = (exp(alpha)) / (1 + exp(alpha))
  P.matrix[2, 2] = 1 - P.matrix[2, 1]  
  
  P.12.param[i, 1] = P.matrix[1, 2]
  P.11.param[i, 1] = P.matrix[1, 1]
  P.22.param[i, 1] = P.matrix[2, 2]
  P.21.param[i, 1] = P.matrix[2, 1]
  
}


#############
##
## Use new P matrix to simulate data
##
#############