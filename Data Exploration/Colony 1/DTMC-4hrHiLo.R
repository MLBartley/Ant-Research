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

#Prep Data - bins interaction data by n seconds, separates by location if applicable

high1 = prep.troph.data(high4, 1)
low1 = prep.troph.data(low4, 1)

high5 = prep.troph.data(high4, 5)
low5 = prep.troph.data(low4, 5)

high30 = prep.troph.data(high4, 30)
low30 = prep.troph.data(low4, 30)

#Visualize Data

par(mfrow = c(1, 1))
plot(high4$start_time,1:nrow(high4),main="High Density Trophallaxis, 4 Hours",
     xlab = "Start Time", 
     ylab = "Number of Interactions")


plot(low4$start_time, 1:nrow(low4), main="Low Density Trophallaxis",
     xlab = "Start Time", 
     ylab = "Number of Interactions", 
     col=low4$Location)
legend(5000, 100, c("Loc1", "Loc4"), lty = c(1,1), col = c("black", "blue"))


#dt model 1 second
theta1 = matrix(data = c(7000, 1, 1, 7000), nrow = 2, ncol = 2, byrow = T)

theta5 = matrix(data = c(100, 1, 1, 100), nrow = 2, ncol = 2, byrow = T)

theta30 = matrix(data = c(70, 1, 1, 70), nrow = 2, ncol = 2, byrow = T)

theta1L = matrix(data = c(74, 1, 1, 74), nrow = 2, ncol = 2, byrow = T)


run.high1.DT = DT.mcmc.troph(y.data = high1$high.y, ant.file = high4,
                             title = "DT Test, 1 sec bin", 
                             a = 2, b = .5, c = 2, d = 5, 
                             theta = theta1, states = 2, n.mcmc = 2000, 
                             delta.t = 1, hours = 4 )

  run.low1 = DT.mcmc.troph(y.data = low1$low.y, ant.file = low4,
                      title = "DT Test, 1 sec bin",
                      a = 2, b = .5, c = 2, d = .5, 
                      theta = theta1L, states = 2, n.mcmc = 2000, 
                      delta.t = 1, hours = 4)

run.low1.entrance = DT.mcmc.troph(y.data = low1$low1.y, ant.file = low4,
                               title = "DT test, 1 sec bin", 
                               a = 2, b = .5, c = 2, d = .5, 
                               theta = theta1L, states = 2, n.mcmc = 5000, 
                               delta.t = 1, hours = 4)

run.low1.queen = mcmc.troph(y.data = low1$low4.y, ant.file = low4,
                            title = "Test", a = 2, b = .5, 
                            theta = theta1L, states = 2, n.mcmc = 5000, 
                            delta.t = 1, hours = 4)


#dt model 5 seconds

run.high5.DT = DT.mcmc.troph(y.data = high5$high.y, ant.file = high4,
                             title = "DT Test, 5 sec bin", 
                             a = 2, b = .5, c = 2, d = .5, 
                             theta = theta5, states = 2, n.mcmc = 2000, 
                             delta.t = 5, hours = 4 )


#dt model 30 seconds


theta30 = matrix(data = c(40, 1, 1, 40), nrow = 2, ncol = 2, byrow = T) 

run.high30.DT = DT.mcmc.troph(y.data = high30$high.y, ant.file = high4,
                             title = "DT Test, 30 sec bin", 
                             a = 2, b = .5, c = 2, d = .5, 
                             theta = theta30, states = 2, n.mcmc = 2000, 
                             delta.t = 30, hours = 4 )

run.low30.DT = DT.mcmc.troph(y.data = low30$low.y, ant.file = low4,
                              title = "DT Test, 30 sec bin", 
                              a = 2, b = .5, c = 2, d = .5, 
                              theta = theta30, states = 2, n.mcmc = 2000, 
                              delta.t = 30, hours = 4 )




run.low30.entrance = DT.mcmc.troph(y.data = low30$low1.y, ant.file = low4,
                                title = "Test", a = 2, b = .5, c= 2, d = .5, 
                                theta = theta30, states = 2, n.mcmc = 5000,
                                delta.t = 30, hours = 4)

run.low30.queen = DT.mcmc.troph(y.data = low30$low4.y, ant.file = low4,
                             title = "Test", a = 1, b = .5,c = 1, d = .5, 
                             theta = theta30, states = 2, n.mcmc = 5000,
                             delta.t = 30, hours = 4)
                            

#cd model 
