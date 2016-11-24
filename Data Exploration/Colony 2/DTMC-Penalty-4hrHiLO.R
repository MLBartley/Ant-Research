#####
##
## 26 September 2016 
## 
## Colony 2 
## 
## Penalized Number of Jumps Model
##
#####################

##Needed Packages
##
##  fdrtool
##  mvtnorm

#Outline
#   write up MCMC for simple model
#   using half normal prior on gamma
#   tuning on sigma^2 tunes our penalization
#   
#   

lambda = c(.01, .08)

P30 = matrix(c(.995, .005, .005, .995), nrow = 2, byrow = T)
sim = sim.DT.troph(tmax = 4 * 60 * 60, delta.t = 30, 
                     start.state = 1, P = P30, lambda = lambda,
                     num.locations = 1)


col2.high <- read.delim("~/Google Drive/PSU/Projects/Ant-Research/Data/Colony2_foraging_high_formatted.txt")
col2.low <- read.delim("~/Google Drive/PSU/Projects/Ant-Research/Data/Colony2_foraging_low_formatted.txt")

high30 = prep.troph.data(col2.high, 30)
low30 = prep.troph.data(col2.low, 30)

states = 2
ant.file = col2.high
hours = 4
a = .08
b = .05
c = .08
d = .05
tau = matrix( c(.01, 0, 
                0, .01), nrow = 2, ncol = 2)
tau.pen = .01
penalty = seq(0.00001, 10, length.out = 10)
n.mcmc = 1000
seconds = 1


bin1.1 = DT.pen.mcmc.troph(y.data = sim$inter.persec, states = 2,
                          ant.file = sim, hours = 4, tau = tau, tau.pen = .01,
                          a = .08, b = .005, c = .08, d = .005,  
                          penalty = penalty[1], n.mcmc = 1000, seconds = 1)

  bin1.1$MSPE


bin1.2 = DT.pen.mcmc.troph(y.data = sim$inter.persec, states = 2,
                           ant.file = sim, hours = 4, tau = tau, tau.pen = .01,
                           a = .08, b = .005, c = .08, d = .005,  
                           penalty = penalty[2], n.mcmc = 1000, seconds = 1)

bin1.2$MSPE



bin1.3 = DT.pen.mcmc.troph(y.data = sim$inter.persec, states = 2,
                           ant.file = sim, hours = 4, tau = tau, tau.pen = .01,
                           a = .08, b = .005, c = .08, d = .005,  
                           penalty = penalty[3], n.mcmc = 1000, seconds = 1)

bin1.3$MSPE


bin1.4 = DT.pen.mcmc.troph(y.data = sim$inter.persec, states = 2,
                           ant.file = sim, hours = 4, tau = tau, tau.pen = .01,
                           a = .08, b = .005, c = .08, d = .005,  
                           penalty = penalty[4], n.mcmc = 1000, seconds = 1)

bin1.4$MSPE
match4 = sim$state == round(bin1.4$X.est)
summary(match4)

bin1.5 = DT.pen.mcmc.troph(y.data = sim$inter.persec, states = 2,
                           ant.file = sim, hours = 4, tau = tau, tau.pen = .01,
                           a = .08, b = .005, c = .08, d = .005,  
                           penalty = penalty[5], n.mcmc = 1000, seconds = 1)

bin1.5$MSPE

match5 = sim$state == round(bin1.5$X.est)
summary(match5)

bin1.6 = DT.pen.mcmc.troph(y.data = sim$inter.persec, states = 2,
                           ant.file = sim, hours = 4, tau = tau, tau.pen = tau,
                           a = .08, b = .005, c = .08, d = .005,  
                           penalty = penalty[6], n.mcmc = 1000, seconds = 1)

match6 = sim$state == round(bin1.6$X.est)
summary(match6)

bin1.7 = DT.pen.mcmc.troph(y.data = sim$inter.persec, states = 2,
                           ant.file = sim, hours = 4, tau = tau, tau.pen = tau,
                           a = .08, b = .005, c = .08, d = .005,  
                           penalty = penalty[7], n.mcmc = 1000, seconds = 1)

match7 = sim$state == round(bin1.7$X.est)
summary(match7)



bin1.8 = DT.pen.mcmc.troph(y.data = sim$inter.persec, states = 2,
                           ant.file = sim, hours = 4, tau = tau, tau.pen = tau,
                           a = .08, b = .005, c = .08, d = .005,  
                           penalty = penalty[8], n.mcmc = 1000, seconds = 1)

match8 = sim$state == round(bin1.8$X.est)
summary(match8)

bin1.9 = DT.pen.mcmc.troph(y.data = sim$inter.persec, states = 2,
                           ant.file = sim, hours = 4, tau = tau, tau.pen = 0.1,
                           a = .08, b = .005, c = .08, d = .005,  
                           penalty = penalty[9], n.mcmc = 1000, seconds = 1, fig.save = T)

bin1.9$MSPE 

match9 = sim$state == round(bin1.9$X.est)
summary(match9)


bin1.10 = DT.pen.mcmc.troph(y.data = sim$inter.persec, states = 2,
                           ant.file = sim, hours = 4, tau = tau, tau.pen = tau,
                           a = .08, b = .005, c = .08, d = .005,  
                           penalty = penalty[10], n.mcmc = 1000, seconds = 1)


match10 = sim$state == round(bin1.10$X.est)
summary(match10)





