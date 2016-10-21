##
## September 2016
## Colony 2 
## 
## 1, 5, 30 second of DT Models
##
###########################
col2.high <- read.delim("~/Google Drive/PSU/Projects/Ant-Research/Data/Colony2_foraging_high_formatted.txt")
col2.low <- read.delim("~/Google Drive/PSU/Projects/Ant-Research/Data/Colony2_foraging_low_formatted.txt")

high1 = prep.troph.data(col2.high, 1)
low1 = prep.troph.data(col2.low, 1)

high5 = prep.troph.data(col2.high, 5)
low5 = prep.troph.data(col2.low, 5)

high15 = prep.troph.data(col2.high, 15)
low15 = prep.troph.data(col2.low, 15)

high30 = prep.troph.data(col2.high, 30)
low30 = prep.troph.data(col2.low, 30)


par(mfrow = c(1, 1))

#Visualize full time series


plot(1:(length(cumsum(high1$high.y))), cumsum(high1$high.y), 
     type = "p", pch = ".", cex = 2,
     xlab = "Time", ylab = "Interactions", 
     main = "Colony 2: High Density, 1s Intervals")

plot(1:(length(cumsum(low1$low.y))), cumsum(low1$low.y), 
     type = "p", pch = ".", cex = 2,
     xlab = "Time", ylab = "Interactions", 
     main = "Colony 2: Low Density, 1s Intervals")

plot(1:(length(cumsum(low1$low.y))), cumsum(low1$low1.y), 
     type = "p", pch = ".", cex = 2,
     xlab = "Time", ylab = "Interactions",
     main = "Colony 2: Low Density, Entrance, 1s Intervals")

plot(1:(length(cumsum(low1$low.y))), cumsum(low1$low4.y), 
     type = "p", pch = ".", cex = 2,
     xlab = "Time", ylab = "Interactions",
     main = "Colony 2: Low Density, Queen's Chamber, 1s Intervals")




#Visualize binned intervals (30 s)

plot(1:(length(cumsum(high30$high.y))), cumsum(high30$high.y), 
     type = "p", pch = ".", cex = 2,
     xlab = "Time", ylab = "Interactions",
     main = "Colony 2: High Density, 30s Intervals")

plot(1:(length(cumsum(low30$low.y))), cumsum(low30$low.y), 
     type = "p", pch = ".", cex = 2,
     xlab = "Time", ylab = "Interactions",
     main = "Colony 2: Low Density, 30s Intervals")


plot(1:(length(cumsum(low30$low.y))), cumsum(low30$low1.y), 
     type = "p", pch = ".", cex = 2,
     xlab = "Time", ylab = "Interactions",
     main = "Colony 2: Low Density, Entrance, 30s Intervals")

plot(1:(length(cumsum(low30$low.y))), cumsum(low30$low4.y), 
     type = "p", pch = ".", cex = 2,
     xlab = "Time", ylab = "Interactions",
     main = "Colony 2: Low Density, Queen's Chamber, 30s Intervals")


### OBSERVATION: High Density, and Low Density Entrance have the 
### most observable switching behaviors.


#dt model 1 second
theta1 = matrix(data = c(8000, 1, 1, 8000), nrow = 2, ncol = 2, byrow = T)

theta1L = matrix(data = c(75, 1, 1, 75), nrow = 2, ncol = 2, byrow = T)



run.high1 = DT.mcmc.troph(y.data = high1$high.y, ant.file = col2.high, 
                       title = "High Density, 1s", a = 2, b = .5, c = 1, d = .5,
                       theta = theta1, states = 2, n.mcmc = 5000, 
                       delta.t = 1, hours = 4)

 run.low1 = DT.mcmc.troph(y.data = low1$low.y, ant.file = low4,
                      title = "Low Density, 1s", a = 2, b = .5, 
                      theta = theta1L, states = 2, n.mcmc = 5000, 
                      delta.t = 1, hours = 4)

run.low1.entrance = mcmc.troph(y.data = low1$low1.y, ant.file = low4,
                               title = "Low Density, Entrance, 1s", a = 2, b = .5, 
                               theta = theta1L, states = 2, n.mcmc = 5000, 
                               delta.t = 1, hours = 4)

run.low1.queen = mcmc.troph(y.data = low1$low4.y, ant.file = low4,
                            title = "Low Density, Queen's Chamber, 1s", a = 2, b = .5, 
                            theta = theta1L, states = 2, n.mcmc = 5000, 
                            delta.t = 1, hours = 4)



#dt model 5 second
theta1 = matrix(data = c(8000, 1, 1, 8000), nrow = 2, ncol = 2, byrow = T)

theta1L = matrix(data = c(75, 1, 1, 75), nrow = 2, ncol = 2, byrow = T)



run.high5 = DT.mcmc.troph(y.data = high5$high.y, ant.file = col2.high, 
                          title = "High Density, 1s", a = 2, b = .5, c = 1, d = .5,
                          theta = theta1, states = 2, n.mcmc = 5000, 
                          delta.t = 5, hours = 4)

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
                        title = "High Density, 30s", a = 2, b = .5, c = 1, d = .5,
                        theta = theta30, states = 2, n.mcmc = 5000,
                        delta.t = 30, hours = 4)


run.low30 = DT.mcmc.troph(y.data = low30$low.y, ant.file = col2.low,
                       title = "Low Density, 30s", a = 2, b = .5, c = 1, d = .5,
                       theta = theta30, states = 2, n.mcmc = 5000,
                       delta.t = 30, hours = 4)

theta30L = matrix(data = c(100, 1, 1, 100), nrow = 2, ncol = 2, byrow = T) 

run.low30.entrance = DT.mcmc.troph(y.data = low30$low1.y, ant.file = col2.low[which(col2.low$Location == 1),],
                                title = "Low Density, Entrance, 30s", a = 2, b = .5, c = 1, d = .5,
                                theta = theta30L, states = 2, n.mcmc = 5000,
                                delta.t = 30, hours = 4)

run.low30.queen = DT.mcmc.troph(y.data = low30$low4.y, ant.file = col2.low[which(col2.low$Location == 4),],
                             title = "Low Density, Queen's, 30s", a = 2, b = .5, c = 1, d = .5,
                             theta = theta30L, states = 2, n.mcmc = 5000,
                             delta.t = 30, hours = 4)

