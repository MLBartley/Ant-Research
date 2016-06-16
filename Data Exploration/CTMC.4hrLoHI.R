############################
##
## 13 June 2016
## Continuous Aspects in MCMC
## Want to consider both number and
## duration of Trophallaxis interactions
##
###########################
#####
# Outline 
#####

  # code (function?) to get data in order 
    # now N_t is number of ants (or pairs of ants) 
    # in trophallaxis at time t

  # function to run new mcmc model

  # Simulated data to ensure model works
  
  # apply to ant data



## Get Data in New format

high4 <- read.csv("./Data/Colony1_trophallaxis_high_density_4hr.csv")

prep.high = prep.troph.pairs(high4)

N = prep.high$pairs2

# hours = 4
# 
# num.high = unique(high4$Ant_ID)
# ant = length(num.high)  #number of ants
# 
# #Embedded Chain and Rate Transitions
# Master = matrix(0, ant, hours * 60 * 60)
# n = 0
# 
# for(i in sort(num.high)){
#   n = n+1
#   high.i = high4[which(high4$Ant_ID == i), ]
# 
#   for(l in 1:nrow(high.i)){
#     start.i = high.i$start_time[l]
#     end.i = high.i$end_time[l]
#     Master[n, start.i:end.i] = 1
#   }
# 
# }
# 
# Master.sum = colSums(Master)
# 
# cont.time = rle(Master.sum)
# 
# EC = cont.time$values
# RT = cont.time$lengths
# 
# 
# ##plot of interactions
# 
# t = cumsum(RT)
# X = EC
# 
# plot(0,0,xlab="t",ylab="State",ylim=c(1,max(X)),xlim=c(0,hours * 60 * 60),type="n")
# lines(x = c(0, t[1]), y = c(X[1], X[1]))
# 
# for(i in 2:length(t)){
#   lines(x = c(t[i-1], t[i]), y = c(X[i], X[i]))
# }
# 
# N = c()
# 
# for(i in 1:length(EC)){
#   N = c(N, rep(EC[i], RT[i]))
# }


###
# Inference
###

#Propose gamma, lambda
  #first simplest case
  # then with additional state based conditions
#calculate R, Q and then P matrices
#accept/reject

n.mcmc = 1000
a = 0.03
b = 0.03
r = 0.03
q = 0.03
tau = c(0.02, 0.02)

data = N[1:1000]

#homes

params = matrix(NA, nrow = 2, ncol = n.mcmc)
rownames(params) <- c("gamma", "lambda")
colnames(params) <- 1:n.mcmc

#initialize

params[1, 1] = rgamma(n = 1, shape = a, rate = b)
params[2, 1] = rgamma(n = 1, shape = r, rate = q)

#posterior likelihood

log.fullcond = function(params, P, data){
  
  sumP = 0 
  
  for(t in 2:length(data)){
    sumP = sumP + log(P[data[t - 1] + 1, data[t] + 1])
  }
  
  loglike = sumP + dgamma(params[1], a, b, log = T) + dgamma(params[2], r, q, log = T)
  
  return(loglike)
}

#MCMC

for(l in 2:n.mcmc){
  

  #propose

  proposal = rexp(n = 2, rate = 1 / params[, l - 1])

  #calculate P matrix
  
  R = matrix(0, nrow = max(N) + 1, ncol = max(N) + 1)
    rownames(R) <- 0:max(N)
    colnames(R) <- 0:max(N)
    
  for( i in 1:nrow(R)){
    
    if( i %% 2 != 0 & i != nrow(R) & i != (nrow(R) - 1)){
       R[i, i + 2] = proposal[1]  
    }
    
    if( i %% 2 != 0 & i != 1 & i != 2){
        R[i, i - 2] = (i - 1) * proposal[2] 
    }
    
  }
    
  P = Pctmc(Q = R, t = 1)
    rownames(P) <- 0:max(N)
    colnames(P) <- 0:max(N)

  #calculate likelihood

  prob.params = exp(log.fullcond(proposal, P, data) - log.fullcond(params[, l - 1], P, data))
    
  #accept/reject

  if(runif(1) < prob.params){
    params[, l] = proposal
  }else{
   params[, l] = params[, l - 1]
  }
  
  
}

  #compile estimates

gamma.est = mean(params[1, ])
lambda.est = mean(params[2, ])


  #plot the estimation runs.

col = c("#120d08", "#bc5356", "#538bbc", "#53bc84")

#gamma
plot(0,0,xlab="MCMC Runs",
     ylab="Gamma (per second)",
     ylim=c(0,max(params)), 
     xlim=c(0,n.mcmc), 
     type="n",
     cex.lab = 1)
  lines(1:n.mcmc, params[1, ], col = col[i])



#lambda
plot(0,0,xlab="MCMC Runs", ylab="Lambda (per second)", 
     ylim=c(0, max(params)), xlim=c(0,n.mcmc), 
     type="n", cex.lab = 1)

  lines(1:n.mcmc, params[2, ], col = col[i])




