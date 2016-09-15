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

lambda = c(.001, .008)

P30 = matrix(c(.995, .005, .005, .995), nrow = 2, byrow = T)
sim30 = sim.DT.troph(tmax = 4 * 60 * 60, delta.t = 30, 
               start.state = 1, P = P30, lambda = lambda,
               num.locations = 1)


#useful for working through function 
#
# 
# y.data = sim30$inter.persec
y.data = sim30$bin.inter
ant.file = sim30
title = "test"
a = .1
b = .5
c = .1
d = .5
theta = theta
states = 2
n.mcmc = 2000
# delta.t = 1
delta.t = 30
hours = 1

data = y.data
Time = length(data)
n = 2
delta = rep(1/n, n)

#needed for final graphic
location = sim30$location
start = sim30$start_time
start = sort(start)
int.num = length(start)
maxtime = hours * 60 * 60



#fit simulated date using function

theta = matrix(data = c(1000, 1, 1, 1000), nrow = 2, ncol = 2, byrow = T) 
test.bin1a = DT.mcmc.troph(y.data = sim30$inter.persec, ant.file = sim30,
                  title = "1000,1", a = .1, b = .5, c = .1, d = .5, 
                  theta = theta, states = 2, n.mcmc = 3000,
                  delta.t = 1, hours = 2)

theta = matrix(data = c(1500, 1, 1, 1500), nrow = 2, ncol = 2, byrow = T) 
test.bin1b = DT.mcmc.troph(y.data = sim30$inter.persec, ant.file = sim30,
                          title = "1500, 1", a = .1, b = .5, c = .1, d = .5, 
                          theta = theta, states = 2, n.mcmc = 3000,
                          delta.t = 1, hours = 2)

theta = matrix(data = c(1700, 1, 1, 1700), nrow = 2, ncol = 2, byrow = T) 
test.bin1c = DT.mcmc.troph(y.data = sim30$inter.persec, ant.file = sim30,
                          title = "test", a = .1, b = .5, c = .1, d = .5, 
                          theta = theta, states = 2, n.mcmc = 3000,
                          delta.t = 1, hours = 2)

theta = matrix(data = c(1900, 1, 1, 1900), nrow = 2, ncol = 2, byrow = T) 
test.bin1d = DT.mcmc.troph(y.data = sim30$inter.persec, ant.file = sim30,
                          title = "1900,1", a = .1, b = .5, c = .1, d = .5, 
                          theta = theta, states = 2, n.mcmc = 3000,
                          delta.t = 1, hours = 2)

theta = matrix(data = c(2000, 1, 1, 2000), nrow = 2, ncol = 2, byrow = T) 
test.bin1e = DT.mcmc.troph(y.data = sim30$inter.persec, ant.file = sim30,
                          title = "2000,1", a = .1, b = .5, c = .1, d = .5, 
                          theta = theta, states = 2, n.mcmc = 3000,
                          delta.t = 1, hours = 2)


theta = matrix(data = c(3000, 1, 1, 3000), nrow = 2, ncol = 2, byrow = T) 
test.bin1e = DT.mcmc.troph(y.data = sim30$inter.persec, ant.file = sim30,
                           title = "3000,1", a = .1, b = .5, c = .1, d = .5, 
                           theta = theta, states = 2, n.mcmc = 3000,
                           delta.t = 1, hours = 2)
############################################################

  
theta30 = matrix(data = c(1, 1, 1, 1), nrow = 2, ncol = 2, byrow = T) 
test.bin30 = DT.mcmc.troph(y.data = sim30$bin.inter, ant.file = sim30,
                      title = "80,1", a = .1, b = .5, c = .1, d = .5,
                      theta = theta30, states = 2, n.mcmc = 5000,
                      delta.t = 30, hours = 2)

theta30 = matrix(data = c(80, 20, 20, 80), nrow = 2, ncol = 2, byrow = T) 
test.bin30a = DT.mcmc.troph(y.data = sim30$bin.inter, ant.file = sim30,
                           title = "80,20", a = .1, b = .5, c = .1, d = .5,
                           theta = theta30, states = 2, n.mcmc = 5000,
                           delta.t = 30, hours = 2)

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
