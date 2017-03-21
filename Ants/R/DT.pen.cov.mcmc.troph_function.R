#' Penalized Discrete Time MCMC Estimation function for Trophallaxis data
#'
#' The purpose of this function is to find MCMC generated estimates of 
#' (1) - the state (X_t = high/low troph rates) of the colony at time t
#' (2) - the specific rates of interaction (lambda) of each state 
#' (3) - the specific rates of state switching (gamma) for the ants
#' (3a) - the probability of moving from one state to another (P matrix) 
#'        generated from gamma values. 
#'
#' @param   y.data, ant.file, title, a, b, c, d, theta, states, n.mcmc, delta.t
#' @return  (1) - estimates of betas, alpha, X, lambda, gamma, P
#'          (2) - 2x2 visual of estimates over time (runs)
#'          (3) - color block state switching graph
#' @export
#' @examples
#'
#' 
#' 


DT.pen.cov.mcmc.troph = function(penalty, covariate, y.data, states, ant.file, 
                                  hours, a, b, c, d, tau, #tau.pen,
                                   n.mcmc, seconds, fig.save, start){
  
#starting values - mostly to keep this all in one place to easily check
X.start = start$X #same number of values as number of seconds
lambda.start = start$lambda #two values
betas.start = start$alpha.beta #four values (for now)

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

Y.1SA.param = matrix(NA, Time, n.mcmc, T)

lambda.param = matrix(data = rep(NA, n * n.mcmc), nrow = n, 
                      ncol = n.mcmc, byrow = T)

## P matrix now varies over time, need new homes

P.11.param = matrix(NA, nrow = Time, ncol = n.mcmc)
P.12.param = matrix(NA, nrow = Time, ncol = n.mcmc)
P.21.param = matrix(NA, nrow = Time, ncol = n.mcmc)
P.22.param = matrix(NA, nrow = Time, ncol = n.mcmc)

alpha.betas.param = matrix(data = NA, 
                           nrow = 4, #NOT GENERALIZED 
                           ncol = n.mcmc)
#if one covariate this is two by two

#note that we're just collecting the off diagonal (non zero) values 
    #needed to calculate P
#Since gammas now vary second-to-second based on covariate(s) 
    #we've split this into two separate variables
    
gamma.param.LH = matrix(NA, nrow = Time, ncol = n.mcmc, 
                     byrow = T)

gamma.param.HL = matrix(NA, Time, n.mcmc)



#probability home for generating X.param, 
gam = matrix(NA, nrow = Time, ncol = n, 
             byrow = T)


penalty.param = matrix(NA, 1, n.mcmc, T) #only used when penalty is not fixed

## Initialize parameters

X.param[,1] = X.start

Y.1SA.param[, 1] = data
 
lambda.param[1, 1] = lambda.start[1] #lambda low
lambda.param[2, 1] = lambda.start[2]# - lambda.start[1] #change in lambda

#lambda.high = lambda.param[1,1] + lambda.param[2, 1] #lambda high, not needed, just a reminder  

alpha.betas.param[1, 1] = betas.start[1] #exp{beta_0}
alpha.betas.param[2, 1] = betas.start[2] #beta_1

alpha.betas.param[3, 1] = betas.start[3] #exp{alpha_0}
alpha.betas.param[4, 1] = betas.start[4] #alpha_1

## IN FUTURE NEED TO GENERALIZE ABOVE CODE FOR ANY NUMBER OF COVARIATES!
## Current Meridith apologizes to Future Meridith, but this is on Her.

#gammas now function of betas (and alpha?)
gamma.param.LH[1, 1] = alpha.betas.param[1,1] * exp(alpha.betas.param[2,1] * covariate[1]) #LH #covariate[1]??? CHECK
gamma.param.HL[1, 1] = alpha.betas.param[3,1] * exp(alpha.betas.param[4,1] * covariate[1]) #HL 

#Don't forget that now Pmatrix and Gamma values change over time
P.matrix = matrix(NA, nrow = n, ncol = n, byrow = T) 


P.matrix[1, 2] = gamma.param.LH[1, 1] * exp(-gamma.param.LH[1, 1] * seconds) / 
                  (gamma.param.LH[1, 1] / (-1 + exp(gamma.param.LH[1, 1])) )
P.matrix[1, 1] = 1 - P.matrix[1,2]
P.matrix[2, 1] = gamma.param.HL[1, 1] * exp(-gamma.param.HL[1, 1] * seconds) /
                    (gamma.param.HL[1, 1] / (-1 + exp(gamma.param.HL[1, 1])) )
P.matrix[2, 2] = 1 - P.matrix[2,1]

#need to send each start value to full P variables (same for all runs)

P.11.param[, 1] = P.matrix[1, 1]
P.12.param[, 1] = P.matrix[1, 2]
P.21.param[, 1] = P.matrix[2, 1]
P.22.param[, 1] = P.matrix[2, 2]
        

      #Uncomment to keep P matrix held constant - check gamma values
      #
      # P.param[ ,1] = c(.995, .005, .005, .995)
      #


#holds all P.parameter values over runs

penalty.param[1, 1] = penalty
#log likelihood

log.fullcond = function(P.param, params, X.param, penalty){
  
  #LogFullCond: [X|P(gamma)] + [e^beta_0] + [betas] + [e^alpha_b] + [alphas]
  #Note that the P.matrix varies over time
  # and is represented as a matrix now 
  # P.param = [P11, P12, P.21, P.22] are column heads

  sumX = 0
  
  for(t in 2:length(data)){
    
    #pull off P matrix for each time
    P.matrix = matrix(data = P.param[t - 1, ], nrow = n, ncol = n, byrow = T)
    
    sumX = sumX + log(P.matrix[X.param[t - 1, l - 1], X.param[t, l - 1]])
  }
  
  loglike = sumX -
    1/penalty * (params[1]^2)  + log(params[1]) -  #e^beta_0 
    1/penalty * (params[3]^2) + log(params[3]) -  #e^alpha_0
    1/penalty * (params[2]^2) - #betas
    1/penalty * (params[4]^2)  #alphas
  return(loglike)
}

accept = 0


for(l in 2:n.mcmc) {
  
  # print out every 100 iterations completed
  if( l %% 100 == 0 ) cat(paste("iteration", l, "complete\n")) 
  
  
  #MH updates - want to propose/accept/reject gammaLH and gammaHL
  
  #adaptive tuning parameter
  # if(l < n.mcmc/2 & l %% 100 == 0){
  # 
  #   sigma = (2.38^2 / 2) * var(log(t(gamma.param[, 1:(l - 1)])))
  #   tau = sigma
  # }

  #proposing - beta.0, beta.1, alpha_0, alpha_1

    #proposal.log = rmvnorm(n = 1, mean = log(alpha.betas.param[c[1, 3], l - 1]), sigma = tau.exp)
    proposal.cov = rmvnorm(n=1, mean =c(log(alpha.betas.param[1, l-1]), 
                                        alpha.betas.param[2, l-1], 
                                        log(alpha.betas.param[3, l-1]), 
                                        alpha.betas.param[4, l-1]), 
                           sigma = tau)
  # proposal.pen = rmvnorm(n = 1, mean = log(penalty.param[, l - 1]), sigma = tau.pen)
  
   #unlog 
  proposal = c(exp(proposal.cov[1]), proposal.cov[2], #beta_0, beta_1
               #exp(proposal.cov[3]), proposal.cov[4] ) #alpha_0, alpha_1
                exp(proposal.cov[3]), 0) #ONLY USED FOR SIMULATION TESTING

  #alpha/betas -> gamma values   
  
  #Note order is not changed as above in commented out line
  gamma.star.LH.t = proposal[1] * exp(proposal[2] * exp(-covariate)) #LH #covariate[1]??? CHECK
  gamma.star.HL.t = proposal[3] * exp(proposal[4] * exp(-covariate)) #HL 
  
  
  
        #uncomment to hold gamma proposal fixed at true values - NEEDS UPDATING
        #
         # theta.star = c(.005, .005)
        # 
  
  
  
  ## Need to take the gamma values 
  ## we've proposed and caluate the 
  ## P matrix for probability of 'jumping' 
  ## between states (at each second)
  ## 
  P.12.star = gamma.star.LH.t * exp(-gamma.star.LH.t * seconds) /
    (gamma.star.LH.t / (-1 + exp(gamma.star.LH.t)) )
  
  
  P.11.star = 1 - P.12.star
  
  P.21.star = gamma.star.HL.t * exp(-gamma.star.HL.t * seconds) /
    (gamma.star.HL.t / (-1 + exp(gamma.star.HL.t)) )
  
  P.22.star = 1 - P.21.star
  

        #Uncommont to hold P matrix constant - NEEDS UPDATING
        #
        # P.star = P.param[,1]
        #
  
  ## NEED TO PULL OFF P param values from last run (for each second)
  ## 
  
  P.11.prev = P.11.param[, l-1]
  P.12.prev = P.12.param[, l-1]
  P.21.prev = P.21.param[, l-1]
  P.22.prev = P.22.param[, l-1]
  
  
  
 P.star = cbind(P.11.star, P.12.star, P.21.star, P.22.star)
 P.previous = cbind(P.11.prev, P.12.prev, P.21.prev, P.22.prev)
  
  #calculate probability
  MHprob = exp(log.fullcond(P.star, proposal, X.param, penalty) -
                 log.fullcond(P.previous, alpha.betas.param[, l-1], X.param, penalty))
  
  if(is.finite(MHprob) == FALSE){MHprob = 0}
  
  
  #accept/reject 
  
  if(runif(1) < MHprob){
    accept = accept + 1
    alpha.betas.param[, l] = proposal
    P.11.param[, l] = P.11.star[1:(hours*60*60)]
    P.12.param[, l] = P.12.star[1:(hours*60*60)]
    P.21.param[, l] = P.21.star[1:(hours*60*60)]
    P.22.param[, l] = P.22.star[1:(hours*60*60)]   
  }else{
    alpha.betas.param[, l] = alpha.betas.param[, l-1]
    P.11.param[, l] = P.11.param[, l-1]
    P.12.param[, l] = P.12.param[, l-1]
    P.21.param[, l] = P.21.param[, l-1]
    P.22.param[, l] = P.22.param[, l-1]  }
  
      
  
  #gibbs updates
  
  ## X Values over time
  
  ##X Parameters

  #Pull off just current probabilities (over time)
  #and then combine into matrix
  
  P.param = cbind(P.11.param[, l], P.12.param[, l], P.21.param[, l], P.22.param[, l])
  
  
  m = matrix(data = 0, nrow = 2, ncol = 2)
  rownames(m) <- c("low", "high")
  colnames(m) <- c("low", "high")
  # number states going from i to j, refreshes every run
  
  lambda.low = lambda.param[1, l - 1]
  lambda.high = lambda.param[2, l - 1] #+lambda.low
  
  ##X Parameters
  
  P.matrix = matrix(data = P.param[1, ], nrow = n, ncol = n, byrow = T)
  
  
  gam[1, 1] = lambda.low ^ data[1] * exp(-lambda.low) * 
    delta[1] * P.matrix[1, X.param[2, l - 1]]
  
  gam[1, 2] = lambda.high ^ data[1] * exp(-lambda.high) * 
    delta[1] * P.matrix[1, X.param[2, l - 1]]
  
  
  
  X.param[1, l] = sample(x = (1:n), size = 1, prob = gam[1,])
 
          ##################
        # X.param[1, l] = X.start[1] 
          ##################
   
  m[X.param[1, l], X.param[1, l]] = m[X.param[1, l], X.param[1, l]] + 1
  
  Y.1SA.param[1, l] = lambda.low * P.matrix[X.param[1, l], 1] + 
    lambda.high * P.matrix[X.param[1, l], 2]
  
  for(t in 2:(Time - 1)){
    P.matrix = matrix(data = P.param[t, ], nrow = n, ncol = n, byrow = T)
    
    
    gam[t, 1] = lambda.low ^ data[t] * exp(-lambda.low) * 
      P.matrix[X.param[t - 1, l - 1], 1] *
      P.matrix[1, X.param[t + 1, l - 1]]
    
    gam[t, 2] = lambda.high ^ data[t] * exp(-lambda.high) * 
      P.matrix[X.param[t - 1, l - 1], 2] *
      P.matrix[2, X.param[t + 1, l - 1]]
    
    
    
    X.param[t, l] = sample(x = (1:n), 1,  prob = gam[t, ]) 
              
                ##################
            # X.param[t, l] = X.start[t]
                ##################
          
    m[X.param[t - 1, l], X.param[t, l]] = m[X.param[t - 1, l], 
                                            X.param[t,l]] + 1
  
    Y.1SA.param[t, l] = lambda.low * P.matrix[X.param[t, l], 1] + 
      lambda.high * P.matrix[X.param[t, l], 2]
    
  }
  
  P.matrix = matrix(data = P.param[Time, ], nrow = n, ncol = n, byrow = T)
  
  
  gam[Time, 1] = lambda.low ^ data[Time] * exp(-lambda.low) * 
    P.matrix[X.param[Time - 1, l - 1], 1] 
  
  gam[Time, 2] = lambda.high ^ data[Time] * exp(-lambda.high) * 
    P.matrix[X.param[Time - 1, l - 1], 2] 
  
  
  X.param[Time, l] = sample(x = 1:n, 1,  prob = gam[Time, ])
  
        ##################
      # X.param[Time, l] = X.start[Time]
       ##################
      
  m[X.param[Time - 1, l], X.param[Time, l]] = m[X.param[Time - 1, l], 
                                                X.param[Time, l]] + 1
  Y.1SA.param[Time, l] = lambda.low * P.matrix[X.param[Time, l], 1] + 
    lambda.high * P.matrix[X.param[Time, l], 2]
  
  ## Lambda Parameters 
  
  lambda.param[1, l] = rgamma(n = 1, shape =
                                sum(data[which(X.param[, l] == 1)]) + a,
                              rate = sum(m[1, ]) + b )
  
  lambda.param[2, l] = rgamma(n = 1, shape =
                                sum(data[which(X.param[, l] == 2)]) + c,
                              rate = sum(m[2, ]) + d)
  
  #Uncomment to hold lambda values fixed at truth
  #
  # lambda.param[, l] = lambda.start
  #
}


#estimation
source("http://www.stat.psu.edu/~mharan/batchmeans.R")

lambda.high =  lambda.param[2, ] #+ lambda.param[1, ]

lambda.est = apply(rbind(lambda.param[1, ], lambda.high), 1, bm) 
lambda.var = apply(lambda.param, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 

# gamma.est = apply(gamma.param, 1, bm) 
# gamma.var = apply(gamma.param, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 


P.11.est = apply(P.11.param, 1, bm)
P.11.var = apply(P.11.param, 1 , quantile, probs = c(0.025, 0.975, na.rm = T))

P.12.est = apply(P.12.param, 1, bm)
P.12.var = apply(P.12.param, 1 , quantile, probs = c(0.025, 0.975, na.rm = T))

P.21.est = apply(P.21.param, 1, bm)
P.21.var = apply(P.21.param, 1 , quantile, probs = c(0.025, 0.975, na.rm = T))

P.22.est = apply(P.22.param, 1, bm)
P.22.var = apply(P.22.param, 1 , quantile, probs = c(0.025, 0.975, na.rm = T))

#Want to create combined estimates over time for visuals/estimation
X.est = matrix(NA, Time, 1, T)
P.11.est = matrix(NA, Time, 1)
P.12.est = matrix(NA, Time, 1)
P.21.est = matrix(NA, Time, 1)
P.22.est = matrix(NA, Time, 1)

for(t in 1:Time ){
  X.est[t, 1] = mean(X.param[t, ]) 
  P.11.est[t, 1] = mean(P.11.param[t, ])
  P.12.est[t, 1] = mean(P.12.param[t, ])
  P.21.est[t, 1] = mean(P.21.param[t, ])
  P.22.est[t, 1] = mean(P.22.param[t, ])
}

sum.it = 0 

for(i in 1:n.mcmc){
  for(t in 1:Time){
   sum.it = sum.it +  (Y.1SA.param[t, i] - data[t])^2
  }
}

MSPE.1SA = 1/n.mcmc * 1/Time * sum.it 

      


#visualization
#
if(fig.save == TRUE){
  pdf(file = paste("./output/", Sys.time(), ".pdf", sep = ""))
  
}

#plot the estimation runs.

col = c("#120d08", "#bc5356", "#538bbc", "#53bc84")

# if(fig.save == T){
#   pdf(file = paste("./output/", Sys.time(), ".pdf", sep = ""))
# }



#Rate Parameters
plot(0,0,xlab="MCMC Runs",
     ylab="Rates",
     ylim=c(0, 2*max(lambda.param)), 
     xlim=c(0,n.mcmc), 
     type="n",
     cex.lab = 1)
lines(1:n.mcmc, (lambda.param[1, ]), col = col[1])
lines(1:n.mcmc,  lambda.param[2, ], col = col[2])

par(mfrow = c(2,2),
    oma = c(0,0,2,0) + 1,
    mar = c(1,1,1,1) + 3)


choose = sample(1:Time, 1) #Choosing one second to explore through all the runs
P.11 = P.11.param[choose, ]
P.21 = P.21.param[choose, ]
P.12 = P.12.param[choose, ]
P.22 = P.22.param[choose, ]

P = cbind(P.11, P.21, P.12, P.22)

#Single P
plot(0,0,xlab="MCMC Runs",
     ylab = "Single Second Probability Matrix for State Switching",
     ylim = c(0, 1), xlim=c(0,n.mcmc),
     type="n", cex.lab = 1)
for(i in 1:(4)){
  lines(1:n.mcmc, P[, i], col = col[i])
}

#Probabilities over time
P.est = cbind(P.11.est, P.12.est, P.21.est, P.22.est)
plot(0,0,xlab="MCMC Runs",
     ylab = "Combined Probability Matrix for State Switching",
     ylim = c(0, 1), xlim=c(0,n.mcmc),
     type="n", cex.lab = 1)
for(i in 1:(4)){
  lines(1:Time, P.est[, i], col = col[i])
}


#Single X
X = X.param[choose, ]
plot(0, 0, xlab = "MCMC Runs", ylab = "Single X", ylim = c(0,max(X)), 
     xlim = c(0,n.mcmc), type = "n", cex.lab = 1)
lines(1:n.mcmc, X, col = col[4])

#States over time
#plot(X.est, type = "l")
plot(round(X.est), type = "l")


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
         col = cols[embedded.chain[j]], border = NA, density = NA)
    
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

if(fig.save == TRUE){
   dev.off()
}


list(X.est = X.est, lambda.est = lambda.est, 
     P.est = c(P.11.est, P.12.est, P.21.est, P.22.est), MSPE = MSPE.1SA, accept = accept)

}


