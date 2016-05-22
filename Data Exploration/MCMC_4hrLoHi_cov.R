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

