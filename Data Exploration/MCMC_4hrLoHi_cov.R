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
#run code and covariates through mcmc 
#how can it be improved?


#high density trophallaxis data

troph.high.4 = read.csv("./Data/Colony1_trophallaxis_high_density_4hr.csv")
troph.low.4 = read.csv("./Data/Colony1_trophallaxis_low_density_4hr.csv")


#in and out data 

inout.high.4 = read.csv("./Data/Colony1_in&out_high_density_4hr.csv")
#only want entrances
inout.high.4 = inout.high.4[which(inout.high.4$Action == "enter"),]

inout.low.4 = read.csv("./Data/Colony1_in&out_low_density_4hr.csv")
#onlly want entrances
inout.low.4 = inout.low.4[which(inout.low.4$Action == "Enter"),]


#visualize high data

sum = sumvis.troph.data(data = troph.high.4,
                        entrance = inout.high.4, hours = 4, 
                        density = "high")

#########
##
## prep trophallaxis data - could be made into function
##
###########

high4 = troph.high.4[seq(1, nrow(troph.high.4), by = 2), ]

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

###################
##
##  prep covariate data - could be made into function
##
###################

forager.arrivals = sort(inout.high.4$time)
forager.arrivals = forager.arrivals[which(duplicated(forager.arrivals) == FALSE)]
forager.arrivals = c(forager.arrivals, 999999)

## Time since Arrivals

time = 4 * 60 * 60
covariate = rep(NA, time)
covariate[1] = 300 #five minutes since arrival

arrival = 1

for(i in 2:time){
  
  if(forager.arrivals[arrival] == i){
    covariate[i] = 0
    arrival = arrival + 1
  }else{
    covariate[i] = covariate[i - 1] + 1
  }
  
}

#bin covariates into similar chunks as the data

max.time = max(troph.high.4$end_time)
delta.t = 60

for(i in 1:length(covariate)){
  c = rep(0, max.time / delta.t)
  mint = 1
  for(t in 1:length(c)){
    c[t] = min(covariate[mint:(mint + delta.t - 1)])
    mint = mint + delta.t
  }
}
high.cov = c


#use model on data
theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 
mu.all = c(2, -1, -0.000004)
sig.all = matrix(data = c(0.2, 0, 0, 
                          0, 0.2, 0, 
                          0, 0, 0.0002), nrow = 3, ncol = 3, byrow = T)


run = mcmc.troph.cov(data = high.y, title = "Test", a = 5, b = 2,
                    theta = theta, states = 2, n.mcmc = 3000, 
                    cov = covariate, mu.cov = mu.all, sig.cov = sig.all)

run2 = mcmc.troph(y.data = high.y, ant.file = troph.high.4, title = "Test", a = 5, b = 2, 
                  theta = theta, states = 2, n.mcmc = 3000, delta.t = 60)

