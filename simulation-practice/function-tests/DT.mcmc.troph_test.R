##
##
## Test code for Ants function: sim.mmpp
## 28 August 2016
##
########################################################

# OUTLINE
# Check discription 
# simulate data for 1 and 30 second intervals
# fit simulated data using funcition
# check given variables vs estimates
# repeat with different starting parameters


#check description
?DT.mcmc.troph

#simulate data using sim.mmpp function
delta.t = 5
lambda = c(.05, .3)

P30 = matrix(c(.997, .003, .003, .997), nrow = 2, byrow = T)
hours = 1
sim30 = sim.DT.troph(tmax = hours * 60 * 60, delta.t = delta.t, 
               start.state = 1, P = P30, lambda = lambda,
               num.locations = 1)


#useful for working through function 
#
# 
y.data = sim30$inter.persec
y.data = sim30$bin.inter
ant.file = sim30
title = "test"
a = .08
b = .005
c = .08
d = .005
theta = theta
states = 2
n.mcmc = 500
# delta.t = 1
delta.t = 5
hours = hours


#fit simulated date using function
#
X = sim30$state
X = sample(x = c(1, 2), size = hours*60*60, replace = T)
X5 = sample(c(1,2), hours*60*60/delta.t, replace = T)
lambda = lambda
P = P30
start = list(X = X, lambda = lambda, P = P)

start5 = list(X = X5, lambda = lambda, P = P)

theta = matrix(data = c(5000, 1, 1, 5000), nrow = 2, ncol = 2, byrow = T) 

test.bin1a = DT.mcmc.troph(y.data = sim30$inter.persec, ant.file = sim30,
                  title = "5000,1", a = .08, b = .005, c = .08, d = .005, 
                  theta = theta, states = 2, n.mcmc = 5000,
                  delta.t = 1, hours = hours, param.start = start)

theta = matrix(data = c(6000, 1, 1, 6000), nrow = 2, ncol = 2, byrow = T) 
test.bin1b = DT.mcmc.troph(y.data = sim30$inter.persec, ant.file = sim30,
                          title = "6000, 1", a = .1, b = .5, c = .1, d = .5, 
                          theta = theta, states = 2, n.mcmc = 3000,
                          delta.t = 1, hours = hours, param.start = start)

    theta = matrix(data = c(1, 1, 1, 1), nrow = 2, ncol = 2, byrow = T) 
test.bin1c = DT.mcmc.troph(y.data = sim30$inter.persec, ant.file = sim30,
                          title = "1, 1", a = .1, b = .5, c = .1, d = .5, 
                          theta = theta, states = 2, n.mcmc = 1000,
                          delta.t = 1, hours = hours, param.start = start)


############################################################

  
theta30 = matrix(data = c(1, 1, 1, 1), nrow = 2, ncol = 2, byrow = T) 
test.bin30 = DT.mcmc.troph(y.data = sim30$bin.inter, ant.file = sim30,
                      title = "1,1", a = .1, b = .5, c = .1, d = .5,
                      theta = theta30, states = 2, n.mcmc = 5000,
                      delta.t = delta.t, hours = hours, param.start = start5)

theta30 = matrix(data = c(80, 20, 20, 80), nrow = 2, ncol = 2, byrow = T) 
test.bin30a = DT.mcmc.troph(y.data = sim30$bin.inter, ant.file = sim30,
                           title = "80,20", a = .1, b = .5, c = .1, d = .5,
                           theta = theta30, states = 2, n.mcmc = 5000,
                           delta.t = 30, hours = hours)

theta30 = matrix(data = c(110, 20, 20, 110), nrow = 2, ncol = 2, byrow = T) 
test.bin30b = DT.mcmc.troph(y.data = sim30$bin.inter, ant.file = sim30,
                           title = "110,20", a = .1, b = .5, c = .1, d = .5,
                           theta = theta30, states = 2, n.mcmc = 5000,
                           delta.t = 30, hours = 2)

theta30 = matrix(data = c(110, 80, 80, 110), nrow = 2, ncol = 2, byrow = T) 
test.bin30c = DT.mcmc.troph(y.data = sim30$bin.inter, ant.file = sim30,
                           title = "test", a = .1, b = .5, c = .1, d = .5,
                           theta = theta30, states = 2, n.mcmc = 5000,
                           delta.t = 30, hours = 2)

theta30 = matrix(data = c(200, 1, 1, 200), nrow = 2, ncol = 2, byrow = T) 
test.bin30d = DT.mcmc.troph(y.data = sim30$bin.inter, ant.file = sim30,
                           title = "test", a = .1, b = .5, c = .1, d = .5,
                           theta = theta30, states = 2, n.mcmc = 5000,
                           delta.t = 30, hours = 2)
