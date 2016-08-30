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
                        start.state = 1, P = P, lambda = lambda, num.locations = 1)


#check per second data
summary(simulate$inter.persec)


##
## Outcome - too many interactions per second (max = 5)
## Solution - lessen lambda parameters
##
##

P = matrix(c(.995, .005, .005, .995), nrow = 2, byrow = T)
lambda = k = c(.01, .04)

simulate = sim.DT.troph(tmax = (1*60*60), delta.t = 30,
                        start.state = 1, P = P, lambda = lambda, num.locations = 1)


#check per second data
summary(simulate$inter.persec)
summary(simulate$bin.inter)

