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
theta1 = matrix(data = c(140, 1, 1, 140), nrow = 2, ncol = 2, byrow = T)

theta1L = matrix(data = c(74, 1, 1, 74), nrow = 2, ncol = 2, byrow = T)



run.high1 = mcmc.troph(y.data = high1$high.y, ant.file = high4, 
                       title = "Test", a = 2, b = .5,
                       theta = theta1, states = 2, n.mcmc = 5000, 
                       delta.t = 1, hours = 4)

  run.low1 = mcmc.troph(y.data = low1$low.y, ant.file = low4,
                      title = "Test", a = 2, b = .5, 
                      theta = theta1L, states = 2, n.mcmc = 5000, 
                      delta.t = 1, hours = 4)

run.low1.entrance = mcmc.troph(y.data = low1$low1.y, ant.file = low4,
                               title = "Test", a = 2, b = .5, 
                               theta = theta1L, states = 2, n.mcmc = 5000, 
                               delta.t = 1, hours = 4)

run.low1.queen = mcmc.troph(y.data = low1$low4.y, ant.file = low4,
                            title = "Test", a = 2, b = .5, 
                            theta = theta1L, states = 2, n.mcmc = 5000, 
                            delta.t = 1, hours = 4)

#dt model 30 seconds


theta30 = matrix(data = c(200, 1, 1, 200), nrow = 2, ncol = 2, byrow = T) 

run.high30 = mcmc.troph(y.data = high30$high.y, ant.file = high4, 
                        title = "Test", a = 2, b = .5, 
                        theta = theta30, states = 2, n.mcmc = 5000,
                        delta.t = 30, hours = 4)


run.low30 = mcmc.troph(y.data = low30$low.y, ant.file = low4,
                       title = "Test", a = 2, b = .5,
                       theta = theta30, states = 2, n.mcmc = 5000,
                       delta.t = 30, hours = 4)

run.low30.entrance = mcmc.troph(y.data = low30$low1.y, ant.file = low4,
                                title = "Test", a = 2, b = .5,
                                theta = theta30, states = 2, n.mcmc = 5000,
                                delta.t = 30, hours = 4)

run.low30.queen = mcmc.troph(y.data = low30$low4.y, ant.file = low4,
                             title = "Test", a = 2, b = .5,
                             theta = theta30, states = 2, n.mcmc = 5000,
                             delta.t = 30, hours = 4)
                            

#cd model 
