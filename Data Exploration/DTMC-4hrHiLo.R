#####
##
## 16 August 2016 
## Need to reconcile differences between discrete and continuous
## time models for 4 hours of ant trophallaxis data.
##
#####################

#Outline
#
##Read in data
##Prep data
##Visualize data
##call two model functions
##run dt model for 30 second intervals
##run dt model for 1 second intervals


#Read in Data

high4 <- read.csv("./Data/Colony1_trophallaxis_high_density_4hr.csv")
low4 <- read.csv("./Data/Colony1_trophallaxis_low_density_4hr.csv")

#Prep Data

high1 = prep.troph.data(high4, 1)
low1 = prep.troph.data(low4, 1)

high30 = prep.troph.data(high4, 30)
low30 = prep.troph.data(low4, 30)

#Visualize Data


#dt model 1 second
theta1 = matrix(data = c(90, 1, 1, 90), nrow = 2, ncol = 2, byrow = T)
theta30 = matrix(data = c(100, 1, 1, 100), nrow = 2, ncol = 2, byrow = T) 

run.high1 = mcmc.troph(y.data = high1$high.y, ant.file = high4, 
                       title = "Test", a = 2, b = .5,
                       theta = theta, states = 2, n.mcmc = 1000, 
                       delta.t = 1, hours = 4)

run.high30 = mcmc.troph(y.data = high30$high.y, ant.file = high4, 
                        title = "Test", a = 2, b = .5, 
                        theta = theta30, states = 2, n.mcmc = 5000,
                        delta.t = 30, hours = 4)

run.low1 = mcmc.troph(y.data = low1$low.y, ant.file = low4,
                      title = "Test", a = 2, b = .5, 
                      theta = theta1, states = 2, n.mcmc = 5000, 
                      delta.t = 1, hours = 4)

run.low30 = mcmc.troph(y.data = low30$low.y, ant.file = low4,
                       title = "Test", a = 2, b = .5,
                       theta = theta30, states = 2, n.mcmc = 5000,
                       delta.t = 30, hours = 4)
                            
