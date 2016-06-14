######################################
##
##  CTMC Coding
##
######################################


#Goal: simulate interaction chain between ants in High Density Data

####
#### Read in Data
####

#Location
#Ant_ID
#Ant_ID_.partner
#start_time
#end_time
#Note: q or Queen is queen ant. Long live the Queen.

high <- read.csv("./Data/Colony1_trophallaxis_high_density_4hr.csv")
hours = 4



head(high)
str(high)


num.high = unique(high$Ant_ID)
ant = length(num.high)  #number of ants


#high = high[seq(1, nrow(high), by = 2), ]
Ti = nrow(high)/2
inter = 1 + #no interactions
  choose(ant, 2)  #two way interactions

#Embedded Chain and Rate Transitions
Master = matrix(0, ant, hours * 60 * 60)
n = 0

for(i in sort(num.high)){
  n = n+1
  high.i = high[which(high$Ant_ID == i), ]
  
  for(l in 1:nrow(high.i)){
    start.i = high.i$start_time[l]
  end.i = high.i$end_time[l]
  Master[n, start.i:end.i] = 1
  }
  
}

Master.sum = colSums(Master)

cont.time = rle(Master.sum)

EC = cont.time$values
RT = cont.time$lengths



# #transitions rates (not probabilities!)
# alpha = matrix(1, n, n, Ti) #(1->2, 2->1)
# 
# for(i in 1:n){
#     alpha[i, i] = 0
#   }
# 
# v= rep(NA, n) #rate for state i
# 
# for(i in 1:n){
#   v[i] = sum(alpha[i, ])
# }
# 
# P = matrix(NA, n, n) #transition probability matrix
# 
# for(i in 1:n){
#   for(k in 1:n){
#     P[i, k] = alpha[i, k] / v[i]
#   }
# }


##plot of interactions


t = cumsum(RT)  
X = EC

plot(0,0,xlab="t",ylab="State",ylim=c(1,max(X)),xlim=c(0,hours * 60 * 60),type="n")
lines(x = c(0, t[1]), y = c(X[1], X[1]))

for(i in 2:length(t)){
  lines(x = c(t[i-1], t[i]), y = c(X[i], X[i]))
}


##Estimating alphas
##lots of vectors are 13 units long bc the states range from 0 to 12

n.states = length(unique(EC))
Alpha = matrix(0, 13, 13)

M = matrix(0, 13, 13) #number of state changes from i to j
for(i in 1:(length(EC)-1)){
   M[EC[i]+1, EC[i +1]+1] =  M[EC[i]+1, EC[i + 1]+1] + 1
}

Time = rep(NA, 13)

for(i in 1:length(Time)){
  Time[i] = sum(RT[which(EC == i-1)])
}

for(i in 1:length(Time)){
  for(j in 1:length(Time)){
    Alpha[i,j] = M[i, j] / Time[i]
  }
}



