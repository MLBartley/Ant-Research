################################
##
## 20 May 2016 
## 
## Want to fit 4 hours low/high density 
## data to model. 
##
## Issue #3 on GitHub
################################

#Outline

#load data
#load covariates
#visualize data and covariates together
#prep data using 
#run code and covariates through mcmc 
#how can it be improved?


#high density trophallaxis data

troph.high.4 = read.csv("./Data/Colony1_trophallaxis_high_density_4hr.csv")
troph.low.4 = read.csv("./Data/Colony1_trophallaxis_low_density_4hr.csv")


#in and out data 

inout.high.4 = read.csv("./Data/Colony1_in&out_high_density_4hr.csv")
inout.low.4 = read.csv("./Data/Colony1_in&out_low_density_4hr.csv")

#visualize high data

sumhigh = sumvis.troph.data(data = troph.high.4,
                        entrance = inout.high.4, hours = 4, 
                        density = "high")

sumlow = sumvis.troph.data(data = troph.low.4, 
                           entrance = inout.low.4, 
                           hours = 4, density = "low")

#########
##
## prep trophallaxis data 
##
###########

high.prep = prep.troph.data(data = troph.high.4, delta.t = 60)
low.prep = prep.troph.data(data = troph.low.4, delta.t = 60)

###################
##
##  prep covariate data - could be made into function
##
###################

high.in.prep = prep.inout.data(data = inout.high.4, delta.t = 60, hours = 4)
low.in.prep = prep.inout.data(data = inout.low.4, delta.t = 60, hours = 4)

##################
##
## bring in simulated data
##
##################







###################
##
## use model on data
##
###################

theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 
mu.all = c(-2.2, -1,  -0.003)
sig.all = matrix(data = c(0.2, 0, 0, 
                          0, 0.2, 0, 
                          0, 0, 0.002), nrow = 3, ncol = 3, byrow = T)
tau = c(0.2, 0.2, 0.002)

run.high = mcmc.troph(y.data = high.prep$high.y, ant.file = troph.high.4, title = "Test", a = 2, b = 2, 
                  theta = theta, states = 2, n.mcmc = 5000, delta.t = 60, hours = 4)

run.high.cov = mcmc.troph.cov(y.data = high.prep$high.y, ant.file = troph.high.4,
                    inout.file = inout.high.4, title = "Test", a = 5, b = 2,
                    theta = theta, states = 2, n.mcmc = 5000, 
                    cov = high.in.prep$cov, mu.cov = mu.all, 
                    sig.cov = sig.all, tau, delta.t = 60)


run.low = mcmc.troph(y.data = low.prep$low.y, ant.file = troph.low.4,
                      title = "test", a = 2, b = 2, theta = theta, 
                      states = 2, n.mcmc = 5000, delta.t = 60, hours = 4)

  run.low.cov = mcmc.troph.cov(y.data = low.prep$low.y, ant.file = troph.low.4, 
                          inout.file = inout.low.4, title = "test", a = 5, 
                          b = 2, theta = theta, states = 2, n.mcmc = 5000, 
                          cov = low.in.prep$cov, mu.cov = mu.all, 
                          sig.cov = sig.all, tau, delta.t = 60)

