#####
##
## 27 March 2016 
## Use MCMC approach to explore 4 hour ant trophallaxis data.
## Have both low and high density data.
##
#####################


####
#### Read in Data
####


##Interaction Data
  #Location
  #Ant_ID
  #Ant_ID_.partner
  #start_time
  #end_time
  #Note: q or Queen is queen ant. Long live the Queen.

high4 <- read.csv("./Data/Colony1_trophallaxis_high_density_4hr.csv")
low4 <- read.csv("./Data/Colony1_trophallaxis_low_density_4hr.csv")

head(high4)
head(low4)
str(high4)
str(low4)

##In/Out Movement Data
  #Ant_ID
  #Time
  #Action (enter/exit Chamber 1)

in.out.high = read.csv("Data/Colony1_in&out_high_density_4hr.csv")
in.out.low = read.csv("Data/Colony1_in&out_low_density_4hr.csv")

c =  prep.inout.data(data = in.out.low, delta.t = 60, hours = 4)
#########################################################
##
## 
## Summary of Data
##
##
#########################################################

#### A look at in/out data 

num.inout.high = unique(in.out.high$Ant_ID)
length(num.inout.high) 

table(in.out.high$Ant_ID)

num.inout.low = unique(in.out.low$Ant_ID)
length(num.inout.low)

table(in.out.low$Ant_ID)
  

# only entrance times
cov.high = in.out.high[which(in.out.high$Action == "enter"), ]
cov.low = in.out.low[which(in.out.low$Action == "Enter"), ]



#### How many ants interacting in each?

num.high = unique(high4$Ant_ID)
length(num.high) 

table(high4$Ant_ID)

par(mfrow = c(1, 1))
hist(table(high4$Ant_ID), xlab = "Count",
     main = "Interactions per Ant: High Density, 4 Hours", 
     breaks = 20)
abline(v = mean(table(high4$Ant_ID)), lty = 3, col="red", lwd = 3)


num.low = unique(low4$Ant_ID)
length(num.low) 

table(low4$Ant_ID)

hist(table(low4$Ant_ID), xlab = "Count",
     main = "Interactions per Ant:Low Density", 
     breaks = 20)
abline(v=mean(table(low4$Ant_ID)), lty = 3, col="red", lwd = 3)


par(mfrow=c(1,2), oma = c(0, 0, 2, 0))

num.low1 = unique(low4$Ant_ID[which(low4$Location == 1)])
length(num.low1)

inter.num1 = table(low4$Ant_ID[which(low4$Location == 1)])
inter.num1 #table of interactions per ant by ID

num.low4 = unique(low4$Ant_ID[which(low4$Location == 4)])
length(num.low4) 

inter.num4 = table(low4$Ant_ID[which(low4$Location == 4)])
inter.num4

xlim = max(inter.num1, inter.num4) #useful so both histograms are on same xlim scale

hist(table(low4$Ant_ID[which(low4$Location == 1)]), xlab = "Count",
     main = "",
     #main = "Interactions per Ant:Low Density, Loc 1", 
     breaks = 20, col = "black",
     xlim = c(0, xlim))

hist(table(low4$Ant_ID[which(low4$Location == 4)]), xlab = "Count",
     #main = "Interactions per Ant:Low Density, Loc 4", 
     main = "",
     breaks = 15, col = "blue", 
     xlim = c(0, xlim))

mtext("Interaction per Ant, by Location", outer = TRUE, cex=1.5)

#########################################################
##
## Get rid of duplicate entries
## All entries recorded twice with each ant in main/partner position
##
#########################################################

high4 = high4[seq(1, nrow(high4), by = 2), ]
low4 = low4[seq(1, nrow(low4), by = 2), ]

#########################################################
##
## plot time of trophylaxis events like a counting process
##
#########################################################

#order data frame by start time so plot works better
high4 = high4[order(high4$start_time), ]
low4 = low4[order(low4$start_time), ]


## Separate Low Density by Location

low4.1 = low4[which(low4$Location == 1), ]
  
low4.4 = low4[which(low4$Location == 4), ]
  low4.4 = low4.4[order(low4.4$start_time), ]


par(mfrow = c(1, 1))
plot(high4$start_time,1:nrow(high4),main="High Density Trophallaxis, 4 Hours",
     xlab = "Start Time", 
     ylab = "Number of Interactions")
points(cov.high$time, rep(0, nrow(cov.high)), pch=8, col="red")

    
plot(low4$start_time, 1:nrow(low4), main="Low Density Trophallaxis",
     xlab = "Start Time", 
     ylab = "Number of Interactions", 
     col=low4$Location)
     legend(5000, 100, c("Loc1", "Loc4"), lty = c(1,1), col = c("black", "blue"))
points(cov.low$time, rep(0, nrow(cov.low)), pch=8, col="red")
  
  
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
    plot(low4.1$start_time, 1:nrow(low4.1), main="Location 1",
         xlim = c(0, max(low4$end_time)),
         xlab = "Start Time", 
         ylab = "Number of Interactions")
points(cov.low$time, rep(0, nrow(cov.low)), pch=8, col="red")
  
    plot(low4.4$start_time, 1:nrow(low4.4), main="Location 4",
         xlim = c(0, max(low4$end_time)),
         xlab = "Start Time", 
         ylab = "Number of Interactions",
         col = "blue")
  points(cov.low$time, rep(0, nrow(cov.low)), pch=8, col="red")
  
      mtext("Low Density Trophallaxis", outer = TRUE, cex = 1.5 )

#########################################################
##
## Bin Interactions into time chuncks
##
#########################################################
delta.t = 60
  
#High Density Data
for(i in 1:nrow(high4)){
  tmp = high4$start_time 
  y = rep(0, max(high4$end_time) / delta.t)
  mint = 0
  for(t in 1:length(y)){
    y[t] = length(which(tmp > mint & tmp <= mint + delta.t))
    mint = mint + delta.t
  }
}
high.y = y

#Low Density Data: Both Locations
for(i in 1:nrow(low4)){
  tmp = low4$start_time
  y = rep(0, max(high4$end_time) / delta.t) #max time same ~7200
  mint = 0
  for(t in 1:length(y)){
    y[t] = length(which(tmp > mint & tmp <= mint + delta.t))
    mint = mint + delta.t
  }
}
low.y = y

#Low Density Data: Location 1
for(i in 1:nrow(low4)){
  tmp = low4.1$start_time
  y = rep(0, max(high4$end_time) / delta.t)
  mint = 0
  for(t in 1:length(y)){
    y[t] = length(which(tmp > mint & tmp <= mint + delta.t))
    mint = mint + delta.t
  }
}
low1.y = y


#Low Density Data: Location 4
for(i in 1:nrow(low4)){
  tmp = low4.4$start_time
  y = rep(0, max(high4$end_time) / delta.t)
  mint = 0
  for(t in 1:length(y)){
    y[t] = length(which(tmp > mint & tmp <= mint + delta.t))
    mint = mint + delta.t
  }
}
low4.y = y


#check all interactions captured
sum(high.y)
nrow(high4)

sum(low.y)
nrow(low4)


#########################################################
##
## Stop! It's MCMC Time
##
#########################################################

#source("functions.R")

theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 

mu.all = c(2, -1, -0.000004)
sig.all = matrix(data = c(0.2, 0, 0, 
                          0, 0.2, 0, 
                          0, 0, 0.0002), nrow = 3, ncol = 3, byrow = T)
tau = c(0.2, 0.2, 0.002)

###
### Two States
###

out.high = mcmc.troph(data = high.y, title = "High Density, 4 Hrs", a = 5, b = 2, 
                      theta = theta, states = 2, n.mcmc = 3000)

out.high.cov = mcmc.troph.cov(data = high.y, title = "High Density w/ Cov, 4 Hrs", 
                              a = 5, b = 2, theta = theta, states = 2, n.mcmc = 1000,
                              cov = cov, mu.cov = mu.all, sig.cov = sig.all )

out.low = mcmc.troph(y.data = low.y, ant.file = low4, title = "Low Density", a = 5, b = 2, 
                     theta = theta, states = 2, n.mcmc = 3000, delta.t = 60, hours = 4)

out.low.cov = mcmc.troph.cov(y.data = low.y, ant.file = low4,
                             inout.file = in.out.low, title = "Test",
                             a = 5, b = 2, theta = theta, states = 2,
                             n.mcmc = 3000, cov = c$cov, mu.cov = mu.all, 
                             sig.cov = sig.all, tau = tau, delta.t = 60, 
                             hours = 4 )

out.low1 = mcmc.troph(data = low1.y, title = "Low Density, Location 1", a = 5, b = 2, 
                      theta = theta, states = 2, n.mcmc = 5000, hours = 4)

out.low4 = mcmc.troph(data = low4.y, title = "Low Density, Location 4",  a = 5, b = 2, 
                      theta = theta, states = 2, n.mcmc = 5000)


#########################################################
##
## Fancy Plots with Background Colors
##
#########################################################
par(mfrow = c(1, 1))

##High Density - 4 Hours
plot(high4$start_time, 1:nrow(high4), main="High, 4 Hours", 
     xlab="Seconds", xlim = c(0,max(high4$end_time)))
##    plot(one.day,1:length(one.day),main=day,xlab="Minutes")
states = out.high$X.est
rr = rle(states[,1])
rr$values = round(rr$values, digits = 0)
embedded.chain = rr$values
cs = c(0,cumsum(rr$lengths))*delta.t - delta.t
cols=c('#FF000022','#0000FF22')
for(j in 1:length(embedded.chain)){
  rect(cs[j],0,cs[j + 1],nrow(high4), 
       col=cols[embedded.chain[j]], density=NA)
}
points(cov.high$time, rep(0, nrow(cov.high)), 
       pch=8, col="forestgreen")

##    axis(4,pretty(c(0,100)),col="green")

#Low Density - 4 Hours

plot(low4$start_time, 1:nrow(low4), main="Low", xlab="Seconds", xlim=c(0,max(high4$end_time)))
##    plot(one.day,1:length(one.day),main=day,xlab="Minutes")
states = out.low$X.est
rr=rle(states[,1])
rr$values = round(rr$values, digits = 0)
embedded.chain=rr$values
cs=c(0,cumsum(rr$lengths))*delta.t - delta.t
cols=c('#FF000022','#0000FF22')
for(j in 1:length(embedded.chain)){
  rect(cs[j],0,cs[j+1],nrow(low4), col=cols[embedded.chain[j]] , density=NA)
}
points(cov.low$time, rep(0, nrow(cov.low)), 
       pch=8, col="forestgreen")
#Low Density - Location 1 

plot(low4.1$start_time, 1:nrow(low4.1), main="Low, Loc 1", xlab="Seconds", xlim=c(0,max(high4$end_time)))
##    plot(one.day,1:length(one.day),main=day,xlab="Minutes")
states = out.low1$X.est
rr=rle(states[,1])
rr$values = round(rr$values, digits = 0)
embedded.chain=rr$values
cs=c(0,cumsum(rr$lengths))*delta.t - delta.t
cols=c('#FF000022','#0000FF22')
for(j in 1:length(embedded.chain)){
  rect(cs[j],0,cs[j+1],nrow(low4.1), col=cols[embedded.chain[j]] , density=NA)
}
points(cov.low$time, rep(0, nrow(cov.low)), 
       pch=8, col="forestgreen")

#Low Density - Location 4
plot(low4.4$start_time, 1:nrow(low4.4), main="Low, Loc 4", xlab="Seconds", xlim=c(0,max(high4$end_time)))
##    plot(one.day,1:length(one.day),main=day,xlab="Minutes")
states = out.low4$X.est
rr=rle(states[,1])
rr$values = round(rr$values, digits = 0)
embedded.chain=rr$values
cs=c(0,cumsum(rr$lengths))*delta.t - delta.t
cols=c('#FF000022','#0000FF22')
for(j in 1:length(embedded.chain)){
  rect(cs[j],0,cs[j+1],nrow(low4.4), col=cols[embedded.chain[j]] , density=NA)
}
points(cov.low$time, rep(0, nrow(cov.low)), 
       pch=8, col="forestgreen")

#########################################################
##
## Checking for Location Lagged Correlation
##
#########################################################
k = 6 #lag

no.delay = lm(low1.y ~ low4.y)
summary(no.delay)

delay = lm(low1.y[(1+k):length(low1.y)]  ~ low4.y[1:length(low4.y)-k])
summary(delay)


