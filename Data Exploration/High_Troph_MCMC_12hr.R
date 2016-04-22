#####
##
## 20 April 2016 
## Use MCMC approach to explore 12 hour ant trophallaxis data.
## Have only high density data.
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

high12 <- read.csv("Data/Colony_1_trophallaxis_high_density_12hr.csv")

head(high12)
str(high12)


#########################################################
##
## 
## Summary of Data
##
##
#########################################################

#How many ants interacting in each?

num.high = unique(high12$Ant_ID)
length(num.high) 

table(high12$Ant_ID)

par(mfrow = c(1, 1))
hist(table(high12$Ant_ID), xlab = "Count",
     main = "Interactions per Ant: High Density, 12 Hrs", 
     breaks = 20)
abline(v = mean(table(high12$Ant_ID)), lty = 3, col = "red", lwd = 3)


#########################################################
##
## Get rid of duplicate entries
## All entries recorded twice with each ant in main/partner position
##
#########################################################

high12 = high12[seq(1, nrow(high12), by = 2), ]

#########################################################
##
## plot time of trophylaxis events like a counting process
##
#########################################################

#order data frame by start time so plot works better
high12 = high12[order(high12$start_time), ]


par(mfrow = c(1, 1))
plot(high12$start_time,1:nrow(high12),main="High Density Trophallaxis, 12 Hrs",
     xlab = "Start Time", 
     ylab = "Number of Interactions")

###Add ant entrance times to first four hours?? MLB 20/Apr


#########################################################
##
## Bin Interactions into time chuncks
##
#########################################################
delta.t = 60

#High Density Data
for(i in 1:nrow(high12)){
  tmp = high12$start_time 
  y = rep(0, max(high12$end_time) / delta.t)
  mint = 0
  for(t in 1:length(y)){
    y[t] = length(which(tmp > mint & tmp <= mint + delta.t))
    mint = mint + delta.t
  }
}
high.y = y


#check all interactions captured
sum(high.y)
nrow(high12)


#########################################################
##
## Stop! It's MCMC Time
##
#########################################################

#source("functions.R")
#above now in Ants Package 

theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 


###
### Two States
###

out.high = mcmc.troph(data = high.y, title = "High Density, 12 Hrs", a = 5, b = 2, 
                      theta = theta, states = 2, n.mcmc = 3000)


#########################################################
##
## Fancy Plots with Background Colors
##
#########################################################
par(mfrow = c(1, 1))

##High Density - 4 Hours
plot(high12$start_time, 1:nrow(high12), main="High, 12 Hours", xlab="Seconds", xlim=c(0,max(high12$end_time)))
##    plot(one.day,1:length(one.day),main=day,xlab="Minutes")
states = out.high$X.est
rr=rle(states[,1])
rr$values = round(rr$values, digits = 0)
embedded.chain=rr$values
cs=c(0,cumsum(rr$lengths))*delta.t - delta.t
cols=c('#FF000022','#0000FF22')
for(j in 1:length(embedded.chain)){
  rect(cs[j],0,cs[j+1],nrow(high12), col=cols[embedded.chain[j]] , density=NA)
}
##    axis(4,pretty(c(0,100)),col="green")

####
## Save the Plots
####

jpeg(file="Presentations/High12MCMC.jpeg")

plot(high12$start_time, 1:nrow(high12), main="High, 12 Hours", xlab="Seconds", xlim=c(0,max(high12$end_time)))
##    plot(one.day,1:length(one.day),main=day,xlab="Minutes")
states = out.high$X.est
rr = rle(states[,1])
rr$values = round(rr$values, digits = 0)
embedded.chain = rr$values
cs = c(0,cumsum(rr$lengths))*delta.t - delta.t
cols = c('#FF000022','#0000FF22')
for(j in 1:length(embedded.chain)){
  rect(cs[j],0,cs[j+1],nrow(high12), col=cols[embedded.chain[j]] , density=NA)
}

dev.off()
