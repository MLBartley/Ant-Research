##
## 
## Simulation Test for DT.pen.cov.mcmc.torph_function
## 
## 


# What goes in:
  #Ant data file 
  #Ant Trophallaxis data - # in troph at each second
  #Covariance data - time since forager entered
  #Number of States - 2
  #Hours of observations - 4
  #Number of Runs 
  #Hyperpriors
  #Penalty prior
  #Proposal density tuning
  #Seconds (are the data binned?)
  #Whether or not to save figures
  #Starting values of chains


# What comes out: 
  #
  #
  #
  #
  


#Outline:
  # Simulate ant data
  # Simulate covariance data
  #
  
alpha.beta = c(.005, .001, .005, .001)
lambda = c(.5, 2)

simulated = sim.DT.cov.troph(tmax = 3600, delta.t = 1, 
                             start.state = 1, alpha.beta = alpha.beta, 
                             lambda = lambda, num.locations = 1)

tau = matrix( c(.0001, 0, 0 , 0,
                0, .0001, 0, 0,
                0, 0, .0001, 0, 
                0, 0, 0, .00001), nrow = 4, ncol = 4)
#penalty = exp(seq(-5, 5, by =  1)) 
penalty = exp(-20)
X = simulated$state
lambda = lambda
alpha.beta = alpha.beta
start = list(X = X, lambda = lambda, alpha.beta = alpha.beta)

covariate = simulated$covariate
y.data = simulated$inter.persec
states = 2
ant.file = simulated
hours = 1
a = .08
b = .005
c = .08
d = .005
n.mcmc = 500
seconds = 1
fig.save = T
start = start

model = DT.pen.cov.mcmc.troph(penalty = penalty, 
                              covariate = simulated$covariate,
                              y.data = simulated$inter.persec,
                              states = 2,
                              ant.file = simulated, hours = 1,
                              a = .08, b = .005, c = .08, d = .005,
                              tau = tau, n.mcmc = 1000,
                              seconds = 1, fig.save = FALSE,
                              start = start)
  
   