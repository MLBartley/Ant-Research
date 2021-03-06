remove.packages("Ants")

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

# setwd("~/Google Drive/PSU/Projects/Ant-Research/simulation-practice")
# run above if running code on personal computer

lambda = c(.01, .08)

P30 = matrix(c(.997, .003, .003, .997), nrow = 2, byrow = T)
gamma = c(0.005, 0.005)

pdf(file = paste("./output/", Sys.time(), ".pdf", sep = ""))

sim = sim_DT_troph(time_max = 1 * 60 * 60, delta_t = 5, 
                   switch_rate = gamma, 
                   start_state = 1, state_tpm = P30, int_rate = lambda,
                   num_locations = 1)



#function parameters

tau = matrix( c(.001, 0, 
                0, .001), nrow = 2, ncol = 2)
penalty = exp(seq(-5, 5, by =  1)) 

penalty = exp(seq(-20, -10, by = 1))
# penalty = exp(c(-20, -14, 7))

# penalty = .00001

X = sim$state
X.30 = sim$bin.state
lambda = lambda
gamma = gamma
start = list(X = X, lambda = lambda, gamma = gamma)

#apply funciton to penalty parameters
#
n.mcmc = 5000
seconds = 1

results = lapply(penalty, FUN = DT_pen_mcmc, starts_data = sim$inter_persec, states = 2, ant_file = sim, hours = 1, a = .08, b = .005, c = .08, d = .005, tau = tau, tau.pen = 0, n_mcmc = n.mcmc, delta_t = 1, start = start, fig_save = T, fig_path = "./Comprehensive-Exam-Prep/output_images/", fig_name = "simu_penalty")

sumtable_model(results = results, compare = penalty, 
  file_path = "./Comprehensive-Exam-Prep/output_tables/", 
  file_name = "sim_penalty", model = "penalized")



# , y.data = sim$inter.persec, states = 2,
#                  ant.file = sim, hours = 4, tau = tau, tau.pen = .01,
#                  a = .08, b = .005, c = .08, d = .005,
#                  n.mcmc = n.mcmc, seconds = 1, fig.save = T, start = start)

results.exp = lapply(penalty, FUN = DT.pen.mcmc.troph.expprior, y.data = sim$inter.persec, states = 2,
                     ant.file = sim, hours = 4, tau = tau, tau.pen = .01,
                     a = .08, b = .005, c = .08, d = .005,  
                     n.mcmc = n.mcmc, seconds = 1, fig.save = T, start = start)


results.exp = results


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
  
  # rownames(table) <- c("Truth", "Model 1", "Model 2", "Model 3", 
  #                      "Model 4", "Model 5", "Model 6", "Model 7", 
  #                      "Model 8", "Model 9", "Model 10")
write.csv(x = table, file = paste("./output/", Sys.time(), ".csv", sep = "") )




# 
# 
# start.30 = list(X = X.30, lambda = lambda, gamma = gamma)
# seconds = 5
# 
# results.30 = lapply(penalty, FUN = DT.pen.mcmc.troph, y.data = sim$bin.inter, states = 2,
#                  ant.file = sim, hours = 4, tau = tau, tau.pen = .01,
#                  a = .08, b = .005, c = .08, d = .005,
#                  n.mcmc = n.mcmc, seconds = seconds, fig.save = T, start = start.30)
# 
# 
# 
# lambda.low.est = lambda[1] * seconds
# 
# for(i in 1:length(penalty)){
#   lambda.low.est = c(lambda.low.est, results.30[[i]]$lambda.est[[1]]$est)
# }
# 
# lambda.high.est = lambda[2] * seconds
# 
# for(i in 1:length(penalty)){
#   lambda.high.est = c(lambda.high.est, results.30[[i]]$lambda.est[[2]]$est)
# }
# 
# gamma.low.est = 0.005
# 
# for(i in 1:length(penalty)){
#   gamma.low.est = c(gamma.low.est, results.30[[i]]$gamma.est[[1]]$est)
# }
# 
# gamma.high.est = 0.005
# 
# for(i in 1:length(penalty)){
#   gamma.high.est = c(gamma.high.est, results.30[[i]]$gamma.est[[2]]$est)
# }
# 
# 
# 
# P.11.est = P30[1, 1]
# 
# for(i in 1:length(penalty)){
#   P.11.est = c(P.11.est, results.30[[i]]$P.est[[1]]$est)
# }
# 
# P.12.est = P30[1, 2]
# 
# for(i in 1:length(penalty)){
#   P.12.est = c(P.12.est, results.30[[i]]$P.est[[2]]$est)
# }
# 
# P.21.est = P30[2, 1]
# 
# for(i in 1:length(penalty)){
#   P.21.est = c(P.21.est, results.30[[i]]$P.est[[3]]$est)
# }
# 
# P.22.est = P30[2, 2]
# 
# for(i in 1:length(penalty)){
#   P.22.est = c(P.22.est, results.30[[i]]$P.est[[4]]$est)
# }
# 
# 
# MSPE.est = 0
# 
# for(i in 1:length(penalty)){
#   MSPE.est = c(MSPE.est, results.30[[i]]$MSPE)
# }
# 
# 
# accept = n.mcmc
# 
# for(i in 1:length(penalty)){
#   accept = c(accept, results.30[[i]]$accept)
# }
# 
# table = data.frame(c(0,penalty), lambda.low.est, lambda.high.est, 
#                    gamma.low.est, gamma.high.est,
#                    P.11.est, P.12.est, P.21.est, P.22.est, 
#                    MSPE.est, accept)
# # 
# # rownames(table) <- c("Truth", "Model 1", "Model 2", "Model 3", 
# #                      "Model 4", "Model 5", "Model 6", "Model 7", 
# #                      "Model 8", "Model 9", "Model 10")
# write.csv(x = table, file = paste("./output/", Sys.time(), ".csv", sep = "") )
# 
# 
