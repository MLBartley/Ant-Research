#####
##
## 14 FEB 2016 
## Use MCMC approach to explore new ant trophallaxis data.
## Have both low and high density data.
##
#####################


####
#### Read in Data
####

#Location
#Ant_ID
#Ant_ID_.partner
#start_time
#end_time
#Note: q or Queen is queen ant. Long live the Queen.

high <- read.csv("~/Google Drive/PSU/Projects/Ant-Research/Data/Colony_1_trophallaxis_high_density.csv")
low <- read.csv("~/Google Drive/PSU/Projects/Ant-Research/Data/Colony_1_trophallaxis_low_density.csv")

head(high)
head(low)
str(high)
str(low)


#########################################################
##
## 
## Summary of Data
##
##
#########################################################

#How many ants interacting in each?

num.high = unique(high$Ant_ID)
length(num.high) 

table(high$Ant_ID)

par(mfrow = c(1, 1))
hist(table(high$Ant_ID), xlab = "Count",
     main = "Interactions per Ant: High Density", 
     breaks = 15)
abline(v=mean(table(high$Ant_ID)), lty = 3, col="red", lwd = 3)

num.low = unique(low$Ant_ID)
length(num.low) 

table(low$Ant_ID)
hist(table(low$Ant_ID), xlab = "Count",
     main = "Interactions per Ant:Low Density", 
     breaks = 15)


par(mfrow=c(1,2), oma = c(0, 0, 2, 0))

num.low1 = unique(low$Ant_ID[which(low$Location == 1)])
length(num.low1)

table(low$Ant_ID[which(low$Location == 1)])
hist(table(low$Ant_ID[which(low$Location == 1)]), xlab = "Count",
     main = "",
     #main = "Interactions per Ant:Low Density, Loc 1", 
     breaks = 15, col = "black",
     xlim = c(0, 35))


num.low4 = unique(low$Ant_ID[which(low$Location == 4)])
length(num.low4) 

table(low$Ant_ID[which(low$Location == 4)])
hist(table(low$Ant_ID[which(low$Location == 4)]), xlab = "Count",
     #main = "Interactions per Ant:Low Density, Loc 4", 
     main = "",
     breaks = 15, col = "blue", 
     xlim = c(0, 35))

mtext("Interaction per Ant, by Location", outer = TRUE, cex=1.5)

#########################################################
##
## Get rid of duplicate entries
## All entries recorded twice with each ant in main/partner position
##
#########################################################

high = high[seq(1, nrow(high), by = 2), ]
low = low[seq(1, nrow(low), by = 2), ]

#########################################################
##
## plot time of trophylaxis events like a counting process
##
#########################################################

#order data frame by start time so plot works better
high = high[order(high$start_time), ]
low = low[order(low$start_time), ]


## Separate Low Density by Location

low.1 = low[which(low$Location == 1), ]
  
low.4 = low[which(low$Location == 4), ]
  low.4 = low.4[order(low.4$start_time), ]


par(mfrow = c(1, 1))
plot(high$start_time,1:nrow(high),main="High Density Trophallaxis",
     xlab = "Start Time", 
     ylab = "Number of Interactions")

plot(low$start_time, 1:nrow(low), main="Low Density Trophallaxis",
     xlab = "Start Time", 
     ylab = "Number of Interactions", 
     col=low$Location)
     legend(5000, 100, c("Loc1", "Loc4"), lty = c(1,1), col = c("black", "blue"))

par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
    plot(low.1$start_time, 1:nrow(low.1), main="Location 1",
         xlim = c(0, 7200),
         xlab = "Start Time", 
         ylab = "Number of Interactions")
    
    plot(low.4$start_time, 1:nrow(low.4), main="Location 4",
         xlim = c(0, 7200),
         xlab = "Start Time", 
         ylab = "Number of Interactions",
         col = "blue")
      mtext("Low Density Trophallaxis", outer = TRUE, cex = 1.5 )

#########################################################
##
## Bin Interactions into time chuncks
##
#########################################################
delta.t = 60
  
#High Density Data
for(i in 1:nrow(high)){
  tmp = high$start_time 
  y = rep(0, max(high$end_time) / delta.t)
  mint = 0
  for(t in 1:length(y)){
    y[t] = length(which(tmp > mint & tmp <= mint + delta.t))
    mint = mint + delta.t
  }
}
high.y = y

#Low Density Data: Both Locations
for(i in 1:nrow(low)){
  tmp = low$start_time
  y = rep(0, max(high$end_time) / delta.t) #max time same ~7200
  mint = 0
  for(t in 1:length(y)){
    y[t] = length(which(tmp > mint & tmp <= mint + delta.t))
    mint = mint + delta.t
  }
}
low.y = y

#Low Density Data: Location 1
for(i in 1:nrow(low)){
  tmp = low.1$start_time
  y = rep(0, max(high$end_time) / delta.t)
  mint = 0
  for(t in 1:length(y)){
    y[t] = length(which(tmp > mint & tmp <= mint + delta.t))
    mint = mint + delta.t
  }
}
low1.y = y


#Low Density Data: Location 4
for(i in 1:nrow(low)){
  tmp = low.4$start_time
  y = rep(0, max(high$end_time) / delta.t)
  mint = 0
  for(t in 1:length(y)){
    y[t] = length(which(tmp > mint & tmp <= mint + delta.t))
    mint = mint + delta.t
  }
}
low4.y = y


#check all interactions captured
sum(high.y)
nrow(high)

sum(low.y)
nrow(low)


#########################################################
##
## Stop! It's MCMC Time
##
#########################################################

source("functions.R")

theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 


###
### Two States
###

out.high = mcmc.troph(data = high.y, title = "High Density", a = 5, b = 2, 
                      theta = theta, states = 2, n.mcmc = 3000)

out.low = mcmc.troph(data = low.y, title = "Low Density", a = 5, b = 2, 
                     theta = theta, states = 2, n.mcmc = 1000)


out.low1 = mcmc.troph(data = low1.y, title = "Low Density, Location 1", a = 5, b = 2, 
                      theta = theta, states = 2, n.mcmc = 5000)

out.low4 = mcmc.troph(data = low4.y, title = "Low Density, Location 4",  a = 5, b = 2, 
                      theta = theta, states = 2, n.mcmc = 5000)


#########################################################
##
## Fancy Plots with Background Colors
##
#########################################################
par(mfrow = c(1, 1))

plot(high$start_time, 1:nrow(high), main="High", xlab="Seconds", xlim=c(0,7200))
##    plot(one.day,1:length(one.day),main=day,xlab="Minutes")
states = out.high$X.est
rr=rle(states[,1])
rr$values = round(rr$values, digits = 0)
embedded.chain=rr$values
cs=c(0,cumsum(rr$lengths))*delta.t - delta.t
cols=c('#FF000022','#0000FF22')
for(j in 1:length(embedded.chain)){
  rect(cs[j],0,cs[j+1],nrow(high), col=cols[embedded.chain[j]] , density=NA)
}
##    axis(4,pretty(c(0,100)),col="green")




###
### Three States
###

theta = matrix(data = c(90, 5, 5, 
                        5, 90, 5,
                        5, 5, 90), nrow = 3, ncol = 3, byrow = T) 


out.high = mcmc.troph(data = high.y, title = "High Density", a = 5, b = 2, 
                      theta = theta, states = 3, n.mcmc = 5000)

out.low = mcmc.troph(data = low.y, title = "Low Density", a = 5, b = 2, 
                     theta = theta, states = 3, n.mcmc = 5000)


out.low1 = mcmc.troph(data = low1.y, title = "Low Density, Location 1", a = 5, b = 2, 
                      theta = theta, states = 3, n.mcmc = 5000)

out.low4 = mcmc.troph(data = low4.y, title = "Low Density, Location 4",  a = 5, b = 2, 
                      theta = theta, states = 3, n.mcmc = 5000)


#########################################################
##
## Fancy Plots with Background Colors
##
#########################################################
par(mfrow = c(1, 1))

plot(high$start_time, 1:nrow(high), main="High", xlab="Seconds", xlim=c(0,7200))
##    plot(one.day,1:length(one.day),main=day,xlab="Minutes")
states = out.high$X.est
rr=rle(states[,1])
rr$values = round(rr$values, digits = 0)
embedded.chain=rr$values
cs=c(0,cumsum(rr$lengths))*delta.t - delta.t
cols=c('#FF000022','#0000FF22')
for(j in 1:length(embedded.chain)){
  rect(cs[j],0,cs[j+1],nrow(high), col=cols[embedded.chain[j]] , density=NA)
}
##    axis(4,pretty(c(0,100)),col="green")


#########################################################
##
## Checking Dirichlet Prior vs Posterier
##
##
#########################################################
theta = matrix(data = c(.9, .1, .1, .9), nrow = 2, ncol = 2, byrow = T) 

out.high.01 = mcmc.troph(data = high.y, title = "High Density", a = 5, b = 2, 
                        theta = theta, states = 2, n.mcmc = 1000)
out.low.01 = mcmc.troph(data = low.y, title = "Low Density", a = 5, b = 2, 
                       theta = theta, states = 2, n.mcmc = 1000)


theta = matrix(data = c(9, 1, 1, 9), nrow = 2, ncol = 2, byrow = T) 

out.high.1 = mcmc.troph(data = high.y, title = "High Density", a = 5, b = 2, 
                      theta = theta, states = 2, n.mcmc = 1000)
out.low.1 = mcmc.troph(data = low.y, title = "Low Density", a = 5, b = 2, 
                      theta = theta, states = 2, n.mcmc = 1000)

theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 

out.high.10 = mcmc.troph(data = high.y, title = "High Density", a = 5, b = 2, 
                        theta = theta, states = 2, n.mcmc = 1000)
out.low.10 = mcmc.troph(data = low.y, title = "Low Density", a = 5, b = 2, 
                       theta = theta, states = 2, n.mcmc = 1000)

theta = matrix(data = c(900, 100, 100, 900), nrow = 2, ncol = 2, byrow = T) 

out.high.100 = mcmc.troph(data = high.y, title = "High Density", a = 5, b = 2, 
                        theta = theta, states = 2, n.mcmc = 1000)
out.low.100 = mcmc.troph(data = low.y, title = "Low Density", a = 5, b = 2, 
                       theta = theta, states = 2, n.mcmc = 1000)



par(mfrow = c(2, 4))
      plot(density(rdirichlet(n = 10000, alpha = c(.9, .1))), col = "purple", ylim = c(0, 4), main = "High: .9, .1")
      for(i in 1:4)  
        lines(density(out.high.01$P.run[i,]), col=i)

      plot(density(rdirichlet(n = 10000, alpha = c(9, 1))), col = "purple", ylim = c(0, 4), main = "High: 9, 1")
       for(i in 1:4)  
         lines(density(out.high.1$P.run[i,]), col=i)

      plot(density(rdirichlet(n = 10000, alpha = c(90, 10))), col = "purple", ylim = c(0, 4), main = "High: 90, 10")
      for(i in 1:4)  
        lines(density(out.high.10$P.run[i,]), col=i)
      #lines(density(rdirichlet(n = 10000, alpha = c(90, 10))), col = "purple")
      
      plot(density(rdirichlet(n = 10000, alpha = c(900, 100))), col = "purple", ylim = c(0, 4), xlim = c(0, 1),
           main = "High: 900, 100")
      for(i in 1:4)  
        lines(density(out.high.100$P.run[i,]), col=i)
      #lines(density(rdirichlet(n = 10000, alpha = c(90, 10))), col = "purple")
      
      
    plot(density(rdirichlet(n = 10000, alpha = c(.9, .1))), col = "purple", ylim = c(0, 4), main = "Low: .9, .1")
    for(i in 1:4)  
      lines(density(out.low.01$P.run[i,]), col=i)


      plot(density(rdirichlet(n = 10000, alpha = c(9, 1))), col = "purple", ylim = c(0, 4), main = "Low: 9, 1")
      for(i in 1:4)  
        lines(density(out.low.1$P.run[i,]), col=i)
      #lines(density(rdirichlet(n = 10000, alpha = c(90, 10))), col = "purple")
      
      plot(density(rdirichlet(n = 10000, alpha = c(90, 10))), col = "purple", ylim = c(0, 4), main = "Low: 90, 10")
      for(i in 1:4)  
        lines(density(out.low.10$P.run[i,]), col=i)
      #lines(density(rdirichlet(n = 10000, alpha = c(90, 10))), col = "purple")
      
      plot(density(rdirichlet(n = 10000, alpha = c(900, 100))), col = "purple", ylim = c(0, 4), xlim = c(0, 1),
           main = "Low, 900, 100")
      for(i in 1:4)  
        lines(density(out.low.100$P.run[i,]), col=i)
      #lines(density(rdirichlet(n = 10000, alpha = c(90, 10))), col = "purple")
