#######################
##
## Script for Combined Results to be used
## for June 7th Meeting for EEID Ant Grant
##
#######################

#Outline

#load data
  # 2 hours
  # 4 hours
  # 12 hours

#load covariates
  # 4 hours

#visualize data and covariates together
#prep data using 
#run code and covariates through mcmc 
#how can it be improved?


# trophallaxis data

troph.high.2 = read.csv("./Data/Colony1_trophallaxis_high_density_2hr.csv")
troph.low.2 = read.csv("./Data/Colony1_trophallaxis_low_density_2hr.csv")

troph.high.4 = read.csv("./Data/Colony1_trophallaxis_high_density_4hr.csv")
troph.low.4 = read.csv("./Data/Colony1_trophallaxis_low_density_4hr.csv")

troph.high.12 = read.csv("./Data/Colony1_trophallaxis_high_density_12hr.csv")

#in and out data 

inout.high.4 = read.csv("./Data/Colony1_in&out_high_density_4hr.csv")
inout.low.4 = read.csv("./Data/Colony1_in&out_low_density_4hr.csv")

#visualize high data

sumhigh2 = sumvis.troph.data(data = troph.high.2,
                             entrance = inout.high.4, hours = 2, 
                             density = "high")

sumlow2 = sumvis.troph.data(data = troph.low.2,
                             entrance = inout.low.4, hours = 2, 
                             density = "low")

sumhigh4 = sumvis.troph.data(data = troph.high.4,
                            entrance = inout.high.4, hours = 4, 
                            density = "high")

sumlow4 = sumvis.troph.data(data = troph.low.4, 
                           entrance = inout.low.4, 
                           hours = 4, density = "low")

sumhigh12 = sumvis.troph.data(data = troph.high.12,
                             entrance = inout.high.4, hours = 12, 
                             density = "high")

#########
##
## prep trophallaxis data 
##
###########

high.prep.2 = prep.troph.data(data = troph.high.2, delta = 60)
low.prep.2 = prep.troph.data(data = troph.low.2, delta = 60)

high.prep.4 = prep.troph.data(data = troph.high.4, delta.t = 60)
low.prep.4 = prep.troph.data(data = troph.low.4, delta.t = 60)

high.prep.12 = prep.troph.data(data = troph.high.12, delta.t = 60)

###################
##
##  prep covariate data - could be made into function
##
###################

high.in.prep = prep.inout.data(data = inout.high.4, delta.t = 60, hours = 4)
low.in.prep = prep.inout.data(data = inout.low.4, delta.t = 60, hours = 4)



###################
##
## use model on data
##
###################

theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 

run.high.2 = mcmc.troph(y.data = high.prep.2$high.y, ant.file = troph.high.2, title = "Test", a = 5, b = 2, 
                      theta = theta, states = 2, n.mcmc = 5000, delta.t = 60, hours = 2)


run.low.2 = mcmc.troph(y.data = low.prep.2$low.y, ant.file = troph.low.2,
                     title = "test", a = 5, b = 2, theta = theta, 
                     states = 2, n.mcmc = 5000, delta.t = 60)

run.high.4 = mcmc.troph(y.data = high.prep$high.y, ant.file = troph.high.4, title = "Test", a = 5, b = 2, 
                        theta = theta, states = 2, n.mcmc = 5000, delta.t = 60)


run.low.4 = mcmc.troph(y.data = low.prep$low.y, ant.file = troph.low.4,
                       title = "test", a = 5, b = 2, theta = theta, 
                       states = 2, n.mcmc = 5000, delta.t = 60)

run.high.12 = mcmc.troph(y.data = high.prep$high.y, ant.file = troph.high.4, title = "Test", a = 5, b = 2, 
                        theta = theta, states = 2, n.mcmc = 5000, delta.t = 60)

