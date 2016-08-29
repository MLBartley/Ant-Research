##
##
## Test code for Ants function: sim.DT.troph
## 28 August 2016
##
########################################################

#check description
?sim.DT.troph

P = matrix(c(.995, .005, .005, .995), nrow = 2, byrow = T)
lambda = k = c(.1, .4)

simulate = sim.DT.troph(tmax = (4*60*60), delta.t = 1,
                        start.state = 1, P = P, lambda = lambda, num.locations = 1)
