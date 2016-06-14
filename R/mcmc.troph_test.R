lambda = c(.05, 1)

P = matrix(c(.8, .2, .2, .8), nrow = 2, byrow = T)
sim = sim.mmpp(tmax = 2 * 60 * 60, delta.t = 30, 
               start.state = 1, P = P, lambda = lambda)

P = matrix(c(.99, .01, .01, .99), nrow = 2, byrow = T)
sim = sim.mmpp(tmax = 1 * 60 * 60, delta.t = 1, 
               start.state = 1, P = P, lambda = lambda)
####

delta.t = 60
sum(sim$y)

sim$start_time = sim$t[which(sim$y >= 1)] #only works with per second data

sim$Location = rep(1, length(sim$start_time))

sim$intbin = rep(NA, length(sim$y)/delta.t)

tint = 1

for(i in 1:length(sim$intbin)){
  
  sim$intbin[i] = sum(sim$y[tint:(tint + delta.t - 1)])
  tint = tint + delta.t
}

#######
# what if instead of binning by totaling interactions, we instead average?

sim$average = sim$intbin / delta.t



####
theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 

test = mcmc.troph(y.data = sim$y, ant.file = sim,
                  title = "test", a = .1, b = 2, 
                  theta = theta, states = 2, n.mcmc = 2000,
                  delta.t = 1, hours = 1)

test.bin5 =  mcmc.troph(y.data = sim$intbin, ant.file = sim,
                        title = "test", a = .1, b = 2, 
                        theta = theta, states = 2, n.mcmc = 2000,
                        delta.t = 5, hours = 1)

# test.bin30 = mcmc.troph(y.data = sim$intbin, ant.file = sim,
#                       title = "test", a = .1, b = 2, 
#                       theta = theta, states = 2, n.mcmc = 2000,
#                       delta.t = 30, hours = 1)
 
 test.bin60 = mcmc.troph(y.data = sim$intbin, ant.file = sim,
                       title = "test", a = .1, b = 2,
                       theta = theta, states = 2, n.mcmc = 2000,
                       delta.t = 60, hours = 1)

# test.avg = mcmc.troph(y.data = sim$average, ant.file = sim,
#                       title = "test", a = .1, b = 2, 
#                       theta = theta, states = 2, n.mcmc = 2000,
#                       delta.t = 5, hours = 1)
