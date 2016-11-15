##
##
## Test code for Ants function: DT.pen.mcmc.troph
## 11 October 2016
##
########################################################

# OUTLINE
# Check discription 
# simulate data for 1 and 30 second intervals
# fit simulated data using funcition
# check given variables vs estimates
# repeat with different starting parameters


#check description
?DT.pen.mcmc.troph


#simulate data using sim.mmpp function

lambda = c(.01, .08)

P30 = matrix(c(.995, .005, .005, .995), nrow = 2, byrow = T)
gamma = c(0.005, 0.005)

sim30 = sim.DT.troph(tmax = 4 * 60 * 60, delta.t = 30, 
                     start.state = 1, P = P30, lambda = lambda, 
                     gamma = gamma,
                     num.locations = 1)

 tau = matrix( c(.01, 0, 
                0, .01), nrow = 2, ncol = 2)
 
tau.pen = .01 

X = sim30$state
X.30 = sim30$bin.state
lambda = lambda
gamma = gamma

start = list(X = X, lambda = lambda, gamma = gamma)
start.30 = list(X = X.30, lambda = lambda, gamma = gamma)




#fit simulated date using function

test.bin1a = DT.pen.mcmc.troph(y.data = sim30$inter.persec, states = 2,
                               ant.file = sim30, hours = 4, tau = tau, tau.pen = tau.pen,
                               a = .08, b = .005, c = .08, d = .005,  
                               penalty = 2, n.mcmc = 5000, seconds = 1, fig.save = T)


##############


test.bin30a = DT.pen.mcmc.troph(y.data = sim30$bin.inter, states = 2, 
                            ant.file = sim30, hours = 4, tau = tau, tau.pen = .01, 
                            a = .1, b = .05, c = .1, d = .05,
                          penalty = 100, n.mcmc = 1000,
                           seconds = 30, fig.save = T, start = start.30)


### Compare to non-penalized functions
### 

theta = matrix(data = c(8000, 1, 1, 8000), nrow = 2, ncol = 2, byrow = T) 
test.bin1COMPARE = DT.mcmc.troph(y.data = sim30$inter.persec, ant.file = sim30,
                           title = "8000,1", a = .1, b = .5, c = .1, d = .5, 
                           theta = theta, states = 2, n.mcmc = 1000,
                           delta.t = 1, hours = 4)


theta30 = matrix(data = c(1, 1, 1, 1), nrow = 2, ncol = 2, byrow = T) 
test.bin30COMPARE = DT.mcmc.troph(y.data = sim30$bin.inter, ant.file = sim30,
                           title = "1,1", a = .1, b = .05, c = .1, d = .05,
                           theta = theta30, states = 2, n.mcmc = 5000,
                           delta.t = 30, hours = 4)



