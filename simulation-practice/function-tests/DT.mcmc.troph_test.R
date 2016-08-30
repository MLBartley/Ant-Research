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

lambda = c(.05, 1)

P30 = matrix(c(.8, .2, .2, .8), nrow = 2, byrow = T)
sim30 = sim.DT.troph(tmax = 1 * 60 * 60, delta.t = 30, 
               start.state = 1, P = P30, lambda = lambda,
               num.locations = 1)

P1 = matrix(c(.99, .01, .01, .99), nrow = 2, byrow = T)
sim1 = sim.DT.troph(tmax = 1 * 60 * 60, delta.t = 1, 
               start.state = 1, P = P1, lambda = lambda, num.locations = 1)

#fit simulated date using function

theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 

test.bin1 = DT.mcmc.troph(y.data = sim30$y, ant.file = sim30,
                  title = "test", a = .1, b = 2, c = .1, d = .5, 
                  theta = theta, states = 2, n.mcmc = 2000,
                  delta.t = 1, hours = 1)

# test.bin5 =  mcmc.troph(y.data = sim$intbin, ant.file = sim,
#                         title = "test", a = .1, b = 2, 
#                         theta = theta, states = 2, n.mcmc = 2000,
#                         delta.t = 5, hours = 1)

test.bin30 = mcmc.troph(y.data = sim$intbin, ant.file = sim,
                      title = "test", a = .1, b = 2,
                      theta = theta, states = 2, n.mcmc = 2000,
                      delta.t = 30, hours = 1)
# 
#  test.bin60 = mcmc.troph(y.data = sim$intbin, ant.file = sim,
#                        title = "test", a = .1, b = 2,
#                        theta = theta, states = 2, n.mcmc = 2000,
#                        delta.t = 60, hours = 1)
