##
## 
## Simulation Test for sim.DT.cov.torph_function
## 
## 


# What goes in:
  #


# What comes out: 
  #

alpha.beta = c(.005, .1, .005, .0001)
lambda = c(.5, 2)

test = sim.DT.cov.troph(tmax = 1*60*60, delta.t = 1, start.state = 1, 
                        alpha.beta = alpha.beta, lambda = lambda, 
                        num.locations = 1)


