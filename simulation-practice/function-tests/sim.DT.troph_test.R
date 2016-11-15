##
##
## Test code for Ants function: sim.DT.troph
## 28 August 2016
##
########################################################

#check description
?sim.DT.troph


#CASE 1
P = matrix(c(.995, .005, .005, .995), nrow = 2, byrow = T)
lambda = k = c(.1, .4)

simulate = sim.DT.troph(tmax = (1*60*60), delta.t = 1,
                        start.state = 1, gamma = 0
                        P = P, lambda = lambda, num.locations = 1)


#check per second data
summary(simulate$inter.persec)


##
## Outcome - too many interactions per second (max = 5)
## Solution - lessen lambda parameters
##
##

P = matrix(c(.997, .003, .003, .997), nrow = 2, byrow = T)
lambda = k = c(.01, .08)

simulate = sim.DT.troph(tmax = (4*60*60), delta.t = 30,
                        start.state = 1,
                        P = P, lambda = lambda, num.locations = 1)


#check per second data
summary(simulate$inter.persec)
summary(simulate$bin.inter)


## Case 2 - Using rate matrix to calculate P matrix (with gammas)
gamma = c(0.005, 0.005)
simulate = sim.DT.troph(tmax = (4*60*60), delta.t = 30,
                        start.state = 1, gamma = gamma,
                        P = P, lambda = lambda, num.locations = 1)
