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

cov.data = read.csv("./Data/Colony1_in&out_high_density_4hr.csv")
lambda.f = nrow(cov.data)/max(cov.data$time)
lambda.f

########
##
## Forager arrival times
##
########

#may just use actual covariate data

prep.cov = prep.inout.data(data = cov.data, delta.t = 1, hours = 4)

# hours = 4
# time = hours * 60 * 60
# forager.number = rpois(n = 1, lambda = lambda.f*time)
# forager.arrivals = rep(NA, forager.number)
#   for (i in 1:forager.number){
#     forager.arrivals[i] = rgamma(n = 1, shape = i,
#                                  rate = lambda.f)
#   }
#   
# forager.arrivals = round(forager.arrivals)
# forager.arrivals = sort(forager.arrivals)
# forager.arrivals = forager.arrivals[which(duplicated(forager.arrivals) == FALSE)]
# forager.arrivals = c(forager.arrivals, time + 1)

########
##
## Time since Arrivals
##
########

# covariate = rep(NA, time)
# covariate[1] = 300 #five minutes since arrival
# 
# arrival = 1
# 
# for(i in 2:time){
#   
#   if(forager.arrivals[arrival] == i){
#     covariate[i] = 0
#     arrival = arrival + 1
#   }else{
#     covariate[i] = covariate[i - 1] + 1
#   }
#   
# }

########
##
## Calculate P matrix based on arrival data (covariate)
##
########
covariate = prep.cov$cov
Time = length(covariate)
theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 
n = 2

alpha = -4.6
beta.0 = -4
beta.1 = -0.001

## P matrix now varies over time, needs homes

P.11.param = matrix(NA, nrow = Time, ncol = 1)
P.12.param = matrix(NA, nrow = Time, ncol = 1)
P.21.param = matrix(NA, nrow = Time, ncol = 1)
P.22.param = matrix(NA, nrow = Time, ncol = 1)

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


#visualized P matrices

#P
t = sample(1:Time, 1)
plot(0, 0, xlab="Time", ylab="Single P", ylim=c(0, 1), xlim=c(0,Time), 
     type="n", cex.lab = 1)

lines(1:Time, P.11.param, col = "red")
lines(1:Time, P.12.param, col = "blue")
lines(1:Time, P.21.param, col = "black")
lines(1:Time, P.22.param, col = "green")


#############
##
## Use new P matrix to simulate data
##
#############

lambda = k = c(.006, .1)

sim = sim.mcmc.dynamP(tmax = length(P.22.param), start.state = 1, 
                P11 = P.11.param, P12 = P.12.param, P21 = P.21.param, 
                P22 = P.22.param, lambda = lambda)
points(cov.data$time, rep(0, length(cov.data$time)), col="#53bc84")

P = matrix(c(.999, .001, .001, .999), nrow = 2, byrow = T)
sim.2 = sim.mmpp(tmax = length(P.21.param), delta.t = 1, start.state = 1, P = P, lambda = lambda)

####
##
## Thoughts - 
##
####

# could combine two simulation functions to be more general
# could add code to dynamic simulation function to add arrival times to plot
# 

################################
##
## 20 May 2016 
## 
## Want to fit simulated data to model and retrieve "true" values before
## using model on full data
##
## Issue #2 on GitHub
########################################################################

### should we simulate 'arrival' times for simulated interactions 
    # useful for mcmc.troph.cov where ant.file is needed to
    # plot the interactions and recovered states

sum(sim$y)

sim$start_time = sim$t[which(sim$y >= 1)]

sim$Location = rep(1, length(sim$start_time))


#we've got simulated data and we know the truth:
    # lambda = c(1, 15) for low/high rates
    # n = 2
    # alpha = -4.6
    # beta.0 = -4
    # beta.1 = 0.0001

mu.all = c(-4, -4, -0.001)
sig.all = matrix(data = c(0.2, 0, 0, 
                          0, 0.2, 0, 
                          0, 0, 0.0002), nrow = 3, ncol = 3, byrow = T)
tau = c(0.2, 0.2, 0.002)

recov = mcmc.troph.cov(y.data = sim$y, ant.file = sim, 
                       inout.file = cov.data, title = "Test", 
                       a = 5, b = 2, theta = theta, states = 2, 
                       n.mcmc = 2000, cov = covariate,
                       mu.cov = mu.all , sig.cov = sig.all, tau = tau,
                       delta.t = 1)

