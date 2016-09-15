#####
##
## 5 September 2016 
## 
## Colony 2 - Discrete Model
##
#####################

col2.high <- read.delim("~/Google Drive/PSU/Projects/Ant-Research/Data/Colony2_foraging_high_formatted.txt")
col2.low <- read.delim("~/Google Drive/PSU/Projects/Ant-Research/Data/Colony2_foraging_low_formatted.txt")

high1 = prep.troph.data(col2.high, 1)
low1 = prep.troph.data(col2.low, 1)

high30 = prep.troph.data(col2.high, 30)
low30 = prep.troph.data(col2.low, 30)


par(mfrow = c(1, 1))
plot(1:(length(cumsum(high30$high.y))), cumsum(high30$high.y), 
     type = "p", pch = ".", cex = 2,
     xlab = "Time", ylab = "Interactions")

plot(1:(length(cumsum(low30$low.y))), cumsum(low30$low.y), 
     type = "p", pch = ".", cex = 2,
     xlab = "Time", ylab = "Interactions")

plot(1:(length(cumsum(low30$low.y))), cumsum(low30$low1.y), 
     type = "p", pch = ".", cex = 2,
     xlab = "Time", ylab = "Interactions")

plot(1:(length(cumsum(low30$low.y))), cumsum(low30$low4.y), 
     type = "p", pch = ".", cex = 2,
     xlab = "Time", ylab = "Interactions")


#dt model 1 second
theta1 = matrix(data = c(200, 1, 1, 200), nrow = 2, ncol = 2, byrow = T)

theta1L = matrix(data = c(75, 1, 1, 75), nrow = 2, ncol = 2, byrow = T)



run.high1 = DT.mcmc.troph(y.data = high1$high.y, ant.file = col2.high, 
                       title = "Test", a = 2, b = .5, c = 1, d = .5,
                       theta = theta1, states = 2, n.mcmc = 5000, 
                       delta.t = 1, hours = 4)

run.low1 = mcmc.troph(y.data = low1$low.y, ant.file = low4,
                      title = "Test", a = 2, b = .5, 
                      theta = theta1L, states = 2, n.mcmc = 3000, 
                      delta.t = 1, hours = 4)

run.low1.entrance = mcmc.troph(y.data = low1$low1.y, ant.file = low4,
                               title = "Test", a = 2, b = .5, 
                               theta = theta1L, states = 2, n.mcmc = 3000, 
                               delta.t = 1, hours = 4)

run.low1.queen = mcmc.troph(y.data = low1$low4.y, ant.file = low4,
                            title = "Test", a = 2, b = .5, 
                            theta = theta1L, states = 2, n.mcmc = 3000, 
                            delta.t = 1, hours = 4)

#dt model 30 seconds


theta30 = matrix(data = c(200, 1, 1, 200), nrow = 2, ncol = 2, byrow = T) 

run.high30 = DT.mcmc.troph(y.data = high30$high.y, ant.file = col2.high, 
                        title = "Test", a = 2, b = .5, c = 1, d = .5,
                        theta = theta30, states = 2, n.mcmc = 5000,
                        delta.t = 30, hours = 4)


run.low30 = DT.mcmc.troph(y.data = low30$low.y, ant.file = col2.low,
                       title = "Test", a = 2, b = .5, c = 1, d = .5,
                       theta = theta30, states = 2, n.mcmc = 5000,
                       delta.t = 30, hours = 4)

theta30L = matrix(data = c(100, 1, 1, 100), nrow = 2, ncol = 2, byrow = T) 

run.low30.entrance = DT.mcmc.troph(y.data = low30$low1.y, ant.file = col2.low,
                                title = "Test", a = 2, b = .5, c = 1, d = .5,
                                theta = theta30L, states = 2, n.mcmc = 5000,
                                delta.t = 30, hours = 4)

run.low30.queen = DT.mcmc.troph(y.data = low30$low4.y, ant.file = col2.low,
                             title = "Test", a = 2, b = .5, c = 1, d = .5,
                             theta = theta30L, states = 2, n.mcmc = 5000,
                             delta.t = 30, hours = 4)

