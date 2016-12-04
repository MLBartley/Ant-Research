install.packages("devtools",  repos='http://cran.us.r-project.org')
library("devtools")

install_github("mlbartley/Ant-Research", subdir = "Ants")

# library("ctmcmove")
# # library("fdrtool")
# library("ggplot2")
 library("gtools")
 library("mvtnorm")
# library("roxygen2")
 library("Ants")

lambda = c(.01, .12)

P30 = matrix(c(.995, .005, .005, .995), nrow = 2, byrow = T)
gamma = c(0.005, 0.005)

pdf(file = paste("./output/", Sys.time(), ".pdf", sep = ""))

sim = sim.DT.troph(tmax = 4 * 60 * 60, delta.t = 5, gamma = gamma, 
                   start.state = 1, P = P30, lambda = lambda,
                   num.locations = 1)

dev.off()


#function parameters

tau = matrix( c(.001, 0, 
                0, .001), nrow = 2, ncol = 2)
penalty = seq(0.000000001, .0001, length.out = 10)
# penalty = .00001

X = sim$state
X.30 = sim$bin.state
lambda = lambda
gamma = gamma
start = list(X = X, lambda = lambda, gamma = gamma)

#apply funciton to penalty parameters
#
n.mcmc = 50
seconds = 1

results = lapply(penalty, FUN = DT.pen.mcmc.troph, y.data = sim$inter.persec, states = 2,
                 ant.file = sim, hours = 4, tau = tau, tau.pen = .01,
                 a = .08, b = .005, c = .08, d = .005,  
                 n.mcmc = n.mcmc, seconds = 1, fig.save = T, start = start)





         #what do I want to pull out for table?
# penalty
# lambda estimates - poisson rates
# P estimates
# gamma esitmates?
# MSPE
#

lambda.low.est = lambda[1] * seconds

for(i in 1:length(penalty)){
    lambda.low.est = c(lambda.low.est, results[[i]]$lambda.est[[1]]$est)
}

lambda.high.est = lambda[2] * seconds

for(i in 1:length(penalty)){
  lambda.high.est = c(lambda.high.est, results[[i]]$lambda.est[[2]]$est)
}

gamma.low.est = 0.005

for(i in 1:length(penalty)){
  gamma.low.est = c(gamma.low.est, results[[i]]$gamma.est[[1]]$est)
}

gamma.high.est = 0.005

for(i in 1:length(penalty)){
  gamma.high.est = c(gamma.high.est, results[[i]]$gamma.est[[2]]$est)
}



P.11.est = P30[1, 1]

for(i in 1:length(penalty)){
  P.11.est = c(P.11.est, results[[i]]$P.est[[1]]$est)
}

P.12.est = P30[1, 2]

for(i in 1:length(penalty)){
  P.12.est = c(P.12.est, results[[i]]$P.est[[2]]$est)
}

P.21.est = P30[2, 1]

for(i in 1:length(penalty)){
  P.21.est = c(P.21.est, results[[i]]$P.est[[3]]$est)
}

P.22.est = P30[2, 2]

for(i in 1:length(penalty)){
  P.22.est = c(P.22.est, results[[i]]$P.est[[4]]$est)
}


MSPE.est = 0

for(i in 1:length(penalty)){
 MSPE.est = c(MSPE.est, results[[i]]$MSPE)
}


accept = n.mcmc

for(i in 1:length(penalty)){
  accept = c(accept, results[[i]]$accept)
}
  


table = data.frame(c(0,penalty), lambda.low.est, lambda.high.est, 
                   gamma.low.est, gamma.high.est,
                   P.11.est, P.12.est, P.21.est, P.22.est, 
                   MSPE.est, accept)
  
  rownames(table) <- c("Truth", "Model 1", "Model 2", "Model 3", 
                       "Model 4", "Model 5", "Model 6", "Model 7", 
                       "Model 8", "Model 9", "Model 10")
write.csv(x = table, file = paste("./output/", Sys.time(), ".csv", sep = "") )






start.30 = list(X = X.30, lambda = lambda, gamma = gamma)
seconds = 5

results.30 = lapply(penalty, FUN = DT.pen.mcmc.troph, y.data = sim$bin.inter, states = 2,
                 ant.file = sim, hours = 4, tau = tau, tau.pen = .01,
                 a = .08, b = .005, c = .08, d = .005,
                 n.mcmc = n.mcmc, seconds = seconds, fig.save = T, start = start.30)



lambda.low.est = lambda[1] * seconds

for(i in 1:length(penalty)){
  lambda.low.est = c(lambda.low.est, results.30[[i]]$lambda.est[[1]]$est)
}

lambda.high.est = lambda[2] * seconds

for(i in 1:length(penalty)){
  lambda.high.est = c(lambda.high.est, results.30[[i]]$lambda.est[[2]]$est)
}

gamma.low.est = 0.005

for(i in 1:length(penalty)){
  gamma.low.est = c(gamma.low.est, results.30[[i]]$gamma.est[[1]]$est)
}

gamma.high.est = 0.005

for(i in 1:length(penalty)){
  gamma.high.est = c(gamma.high.est, results.30[[i]]$gamma.est[[2]]$est)
}



P.11.est = P30[1, 1]

for(i in 1:length(penalty)){
  P.11.est = c(P.11.est, results.30[[i]]$P.est[[1]]$est)
}

P.12.est = P30[1, 2]

for(i in 1:length(penalty)){
  P.12.est = c(P.12.est, results.30[[i]]$P.est[[2]]$est)
}

P.21.est = P30[2, 1]

for(i in 1:length(penalty)){
  P.21.est = c(P.21.est, results.30[[i]]$P.est[[3]]$est)
}

P.22.est = P30[2, 2]

for(i in 1:length(penalty)){
  P.22.est = c(P.22.est, results.30[[i]]$P.est[[4]]$est)
}


MSPE.est = 0

for(i in 1:length(penalty)){
  MSPE.est = c(MSPE.est, results.30[[i]]$MSPE)
}


accept = n.mcmc

for(i in 1:length(penalty)){
  accept = c(accept, results.30[[i]]$accept)
}

table = data.frame(c(0,penalty), lambda.low.est, lambda.high.est, 
                   gamma.low.est, gamma.high.est,
                   P.11.est, P.12.est, P.21.est, P.22.est, 
                   MSPE.est, accept)

rownames(table) <- c("Truth", "Model 1", "Model 2", "Model 3", 
                     "Model 4", "Model 5", "Model 6", "Model 7", 
                     "Model 8", "Model 9", "Model 10")
write.csv(x = table, file = paste("./output/", Sys.time(), ".csv", sep = "") )


