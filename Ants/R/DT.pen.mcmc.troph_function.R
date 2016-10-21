#' Penalized Discrete Time MCMC Estimation function for Trophallaxis data
#'
#' The purpose of this function is to find MCMC generated estimates of 
#' (1) - the state (X_t = high/low troph rates) of the colony at time t
#' (2) - the specific rates of interaction (lambda) of each state 
#' (3) - the specific rates of state switching (gamma) for the ants
#' (3a) - the probability of moving from one state to another (P matrix) 
#'        generated from gamma values. 
#'
#' @param y.data, ant.file, title, a, b, c, d, theta, states, n.mcmc, delta.t
#' @return  (1) - estimates of X, lambda, gamma, P
#'          (2) - 2x2 visual of estimates over time (runs)
#'          (3) - color block state switching graph
#' @export
#' @examples
#' theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 
#' out.high = DT.pen.mcmc.troph()
#' 
#' 


DT.pen.mcmc.troph = function(y.data,states, ant.file, 
                                  hours, a, b, c, d, tau,
                                  penalty, n.mcmc, seconds){

#any changes to data

data = y.data
Time = length(data)
n = states
delta = rep(1/n, n)

#needed for final graphic 
location = ant.file$Location 
start = ant.file$start_time
start = sort(start)
int.num = length(start)
maxtime = hours * 60 * 60

#build homes for chains

## Build Homes for X(1:T), lambda(1:n) vector , and P(nXn),  gamma (nxn) matrices

X.param = matrix(data = rep(NA, Time * n.mcmc), nrow = Time, 
                 ncol = n.mcmc, byrow = T)

lambda.param = matrix(data = rep(NA, n * n.mcmc), nrow = n, 
                      ncol = n.mcmc, byrow = T)

P.param = matrix(data = rep(NA, n * n * n.mcmc), nrow = n * n, 
                 ncol = n.mcmc, byrow = T)


#note that we're just collecting the off diagonal (non zero) values needed to calculate P
gamma.param = matrix(NA, nrow = n, ncol = n.mcmc, 
                     byrow = T)

#probability home for generating X.param, 
gam = matrix(NA, nrow = Time, ncol = n, 
             byrow = T)


penalty.param = matrix(NA, 1, n.mcmc, T)

## Initialize parameters

X.param[,1] = sample(1:2, replace = T, Time)

lambda.param[1, 1] = rgamma(n = 1, shape = a, rate = b) #lambda low
lambda.param[2, 1] = rgamma(n = 1, shape = c, rate = d) #change in lambda

lambda.high = lambda.param[1,1] + lambda.param[2, 1] #lambda high, not needed, just a reminder  

gamma.param[1, 1] = rhalfnorm(n = 1, theta = sqrt(pi/2)/penalty) #LH
gamma.param[2, 1] = rhalfnorm(n = 1, theta = sqrt(pi/2)/penalty) #HL 

P.matrix = matrix(NA, nrow = n, ncol = n, byrow = T) 

P.matrix[1, 2] = gamma.param[1, 1] * exp(- gamma.param[1, 1] * seconds)
P.matrix[1, 1] = 1 - P.matrix[1,2]
P.matrix[2, 1] = gamma.param[2, 1] * exp(- gamma.param[2, 1] * seconds)
P.matrix[2, 2] = 1 - P.matrix[1,2]

P.param[,1] = as.vector(t(P.matrix))
#holds all P.parameter values over runs
#

penalty.param[1, 1] = penalty

#log likelihood

log.fullcond = function(P.param, params, data, X.param, penalty){
  
  P.matrix = matrix(data = P.param, nrow = n, ncol = n, byrow = T)
  
  sumX = 0
  
  for(t in 2:length(data)){
    
    sumX = sumX + log(P.matrix[X.param[t-1, l-1], X.param[t, l-1]])
  }
  
  loglike = sumX -
    log(penalty) -
    1 / (2 * penalty) * (params[1]^2 + params[2]^2)  
  
  
  
  return(loglike)
}

accept = 0



for(l in 2:n.mcmc) {
  
  # print out every 10 iterations completed
  if( l %% 100 == 0 ) cat(paste("iteration", l, "complete\n")) 
  
  
  #MH updates - want to propose/accept/reject gammaLH and gammaHL
  
  # #adaptive tuning parameter
  if(l < n.mcmc/2 & l %% 100 == 0){
    
    sigma = 2.38 ^ 2 / 2 * var(log(t(gamma.param[, 1:(l - 1)])))
    tau = sigma
  }
  
  
  proposal = rmvnorm(n = 1, mean = log(gamma.param[, l - 1]), sigma = tau)
  
  
  #unlog 
  theta.star = exp(proposal)
  
  
  ## Need to take the gamma values 
  ## we've proposed and caluate the 
  ## P matrix for probability of 'jumping' 
  ## between states
  
  P.matrix[1, 2] = theta.star[1] * exp(- theta.star[1] * seconds)
  P.matrix[1, 1] = 1 - P.matrix[1,2]
  P.matrix[2, 1] = theta.star[2] * exp(- theta.star[2] * seconds)
  P.matrix[2, 2] = 1 - P.matrix[1,2]
  
  P.star = as.vector(t(P.matrix))
  
  #calculate probability
  MHprob = exp(log.fullcond(P.star, theta.star, data, X.param, penalty) -
                 log.fullcond(P.param[, l - 1], gamma.param[, l-1], data, X.param, penalty))
  
  if(is.finite(MHprob) == FALSE){MHprob = 0}
  
  
  #accept/reject
  
  if(runif(1) < MHprob){
    accept = accept + 1
    gamma.param[, l] = theta.star
    P.param[, l] = as.vector(t(P.star))
    
  }else{
    gamma.param[, l] = gamma.param[, l - 1]
    P.param[, l] = P.param[, l - 1]
  }
  
  
  #gibbs updates
  
  ## X Values over time
  
  ##X Parameters
  
  
  m = matrix(data = 0, nrow = 2, ncol = 2)
  rownames(m) <- c("low", "high")
  colnames(m) <- c("low", "high")
  # number states going from i to j, refreshes every run
  
  lambda.low = lambda.param[1, l - 1]
  lambda.high = lambda.low + lambda.param[2, l - 1]
  
  ##X Parameters
  
  gam[1, 1] = lambda.low ^ data[1] * exp(-lambda.low) * 
    delta[1] * P.matrix[1, X.param[2, l - 1]]
  
  gam[1, 2] = lambda.high ^ data[1] * exp(-lambda.high) * 
    delta[1] * P.matrix[1, X.param[2, l - 1]]
  
  
  
  X.param[1, l] = sample(x = (1:n), size = 1, prob = gam[1,])
  
  m[X.param[1, l], X.param[1, l]] = m[X.param[1, l], X.param[1, l]] + 1
  
  for(t in 2:(Time - 1)){
    
    gam[t, 1] = lambda.low ^ data[t] * exp(-lambda.low) * 
      P.matrix[X.param[t - 1, l - 1], 1] *
      P.matrix[1, X.param[t + 1, l - 1]]
    
    gam[t, 2] = lambda.high ^ data[t] * exp(-lambda.high) * 
      P.matrix[X.param[t - 1, l - 1], 2] *
      P.matrix[2, X.param[t + 1, l - 1]]
    
    
    
    X.param[t, l] = sample(x = (1:n), 1,  prob = gam[t, ]) 
    
    m[X.param[t - 1, l], X.param[t, l]] = m[X.param[t - 1, l], 
                                            X.param[t,l]] + 1
  }
  
  gam[Time, 1] = lambda.low ^ data[Time] * exp(-lambda.low) * 
    P.matrix[X.param[Time - 1, l - 1], 1] 
  
  gam[Time, 2] = lambda.high ^ data[Time] * exp(-lambda.high) * 
    P.matrix[X.param[Time - 1, l - 1], 2] 
  
  
  X.param[Time, l] = sample(x = 1:n, 1,  prob = gam[Time, ])
  
  m[X.param[Time - 1, l], X.param[Time, l]] = m[X.param[Time - 1, l], 
                                                X.param[Time, l]] + 1
  
  ## Lambda Parameters 
  
  lambda.param[1, l] = rgamma(n = 1, shape =
                                sum(data[which(X.param[, l] == 1)]) + a,
                              rate = sum(m[1, ]) + b )
  
  lambda.param[2, l] = rgamma(n = 1, shape =
                                sum(data[which(X.param[, l] == 2)]) + c,
                              rate = sum(m[2, ]) + d)
}


#estimation
source("http://www.stat.psu.edu/~mharan/batchmeans.R")

lambda.est = apply(lambda.param, 1, bm) 
lambda.var = apply(lambda.param, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 

gamma.est = apply(gamma.param, 1, bm) 
gamma.var = apply(gamma.param, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 


P.est = apply(P.param[,1:(l-1)], 1, bm)
P.var = apply(P.param, 1 , quantile, probs = c(0.025, 0.975, na.rm = T))


X.est = matrix(NA, Time, 1, T)

for(t in 1:Time ){
  X.est[t, 1] = mean(X.param[t, ])  
}



#visualization

#plot the estimation runs.

col = c("#120d08", "#bc5356", "#538bbc", "#53bc84")

# if(save == T){
#   pdf(file = paste("./output/", Sys.time(), ".pdf", sep = ""))
# }

par(mfrow = c(2,2),
    oma = c(0,0,2,0) + 1,
    mar = c(1,1,1,1) + 3)


#Rate Parameters
plot(0,0,xlab="MCMC Runs",
     ylab="Rates (per minute)",
     ylim=c(0, 60* max(gamma.param, lambda.param)), 
     xlim=c(0,n.mcmc), 
     type="n",
     cex.lab = 1)
lines(1:n.mcmc, 60 * (lambda.param[1, ] + lambda.param[2, ]), col = col[1])
lines(1:n.mcmc, 60 * lambda.param[1, ], col = col[2])

lines(1:n.mcmc, (60 * gamma.param[1, ]), col = col[3])
lines(1:n.mcmc, (60 * gamma.param[2, ]), col = col[4])

#X params

#Single X
X = X.param[sample(1:Time, 1), ]
plot(0, 0, xlab = "MCMC Runs", ylab = "Single X", ylim = c(0,max(X)), 
     xlim = c(0,n.mcmc), type = "n", cex.lab = 1)
lines(1:n.mcmc, X, col = col[4])

#States over time
#plot(X.est, type = "l")
plot(round(X.est), type = "l")

#P
plot(0,0,xlab="MCMC Runs", ylab = "Probability Matrix for State Switching", ylim = c(0, max(P.param)), xlim=c(0,n.mcmc), 
     type="n", cex.lab = 1)
for(i in 1:(4)){
  lines(1:n.mcmc, P.param[i, ], col = col[i])
}

title(main = "Diagnostic Plots", outer = T)

#########################################################
##
## Fancy Plots with Background Colors
##
#########################################################

par(mfrow = c(1, 1))


if(length(unique(location)) == 1){
  
  ##High Density - 4 Hours
  plot(start, 1:int.num, main="High", 
       xlab="Seconds", ylab = "Cumulative Interaction Count", 
       xlim = c(0, maxtime))
  states = X.est #from code above
  rr = rle(states[,1])
  rr$values = round(rr$values, digits = 0)
  embedded.chain = rr$values
  cs = c(0,cumsum(rr$lengths))*seconds - seconds
  cols=c('#bc535644','#538bbc44')
  for(j in 1:length(embedded.chain)){
    rect(cs[j],0,cs[j + 1],int.num, 
         col = cols[embedded.chain[j]], density = NA)
    
  }
  points(start, 1:int.num, main="Low", xlab="Seconds",
         ylab = "Cumulative Interaction Count", 
         xlim=c(0,maxtime))
}else{
  #Low Density - 4 Hours
  
  plot(start, 1:int.num, main="Low", xlab="Seconds",
       ylab = "Cumulative Interaction Count", 
       xlim=c(0,maxtime))
  states = X.est
  rr=rle(states[,1])
  rr$values = round(rr$values, digits = 0)
  embedded.chain=rr$values
  cs=c(0,cumsum(rr$lengths))*delta.t - delta.t
  cols=c('#bc535644','#538bbc44')
  for(j in 1:length(embedded.chain)){
    rect(cs[j],0,cs[j+1],int.num, col=cols[embedded.chain[j]] , density=NA)
  }
  points(start, 1:int.num, main="Low", xlab="Seconds",
         ylab = "Cumulative Interaction Count", 
         xlim=c(0,maxtime))
  }

list(X.est = X.est, lambda.est = lambda.est, gamma.est = gamma.est, P.est = P.est, P.run = P.param)

}


