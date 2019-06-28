###############################################################################
## This script aims to initialize all simple model details
##
## Created: May 24, 2019
## Updated 1: fixed code so that dbern has right parameter value (#changed to get P[X_{t-1}, 2])
## Updated 2: fixed code so that the P denominator is sum{-i}(1 + e.beta) -- otherwise doesn't add to 1
###############################################################################
library(nimble)

source("./NIMBLE/vignettes/01.01_prepdata_colony2.R")

dat <- out$col2_low4_5$queen_starts_persec #created in 01_prepdata

seconds <- length(dat)

antsCode <- nimbleCode({
  #rates of interactions for low and high (low + diff) colonly-level states
  lambda_l ~ dgamma(a, b)
  lambda_diff ~ dgamma(c, d)

  lambda_h <- lambda_l + lambda_diff

  #P is a function of beta, a rate variable for each state
  e.beta[1:num.states] ~ dmnorm(mvnorm.mean[1:num.states],
                                cov = tau[1:num.states, 1:num.states])

  #probablility that the colony switches between states
  for(i in 1:num.states){

    P[i,i] <- 1 / (1 + sum(e.beta[1:num.states]) - e.beta[i])

    for(j in indx[[i]]){

      P[i, j] <- e.beta[j] / (1 + sum(e.beta[1:num.states]) - e.beta[i])
    }
  }


  state[1] ~ dbern(.5)

  for (t in 2:nSecs) {
    state[t] ~ dbern(prob = P[(state[t - 1] +1), 2])
    #changed to get P[X_{t-1}, 2]
  }

  for (t in 1:nSecs){
    y[t] ~ dpois(lambda_l + lambda_diff * (state[t]))
    # }

    ##need to calculate OSA MSPE within nimble - this way we don't need to monitor states
    # for (t in 1:nSecs){
    y_hat[t] <- (lambda_l * P[(state[t] + 1), 1] +
                   lambda_h * P[(state[t] + 1), 2])
    mspe_diff[t] <- ((y_hat[t]) - y[t])^2
  }

  mspe <- 1/nSecs * sum(mspe_diff[1:nSecs])



}) # end model

# constants: num.states, a, b, c, d, tau, n.Secs
# needs starting values: states, y, P?, lambdas,
## define constants, data, and initial values
num.states = 2
# cov.seconds = 1 #covariates NO
# cov.t = 1
# alpha = 1
# cov.t = cov_col2_low4$cov
# alpha = -0.0001
seconds = length(dat)

indx <- list()
for (i in 1:num.states){
  indx[[i]] <- which(1:num.states != i)
}


data <- list(y = dat)
x.init <- sample(x = c(0, 1), size = seconds, replace = T)
e.beta.init <- exp(c(-4.59, -4.59))

P <- nimArray(NA,
              dim = c(num.states, num.states),
              #cov.sec = 1 if no cov, else = seconds
              init = TRUE)
# tau <- nimMatrix(NA, num.states, num.states, init = TRUE)

# for( t in 1:cov.seconds){
for(i in 1:num.states){
  P[i,i] <- 1 /
    (1 + sum(e.beta.init[i != 1:num.states]))
  # tau[i, i] <- exp(10)

  for(j in which(1:num.states != i)){
    P[i, j] <- e.beta.init[j] /
      (1 + sum(e.beta.init[i != 1:num.states]))
    # tau[i, j] <- 0
  }
}
# }

p.init <- P

# constants <- list(a = 1, b = 1, c = 1, d = 1,
                  # mvnorm.mean = rep(0, num.states),
                  # # tau = tau,
                  # num.states = num.states,
                  # nSecs = seconds,
                  # indx = indx,
                  # tau = matrix(c(penalty, 0, 0, penalty), 2)
                  # )

inits <- list( lambda_l = 0.007, lambda_diff = 0.05, lambda_h = .007 + .05,
               e.beta = e.beta.init,  state = x.init, P = p.init, y_hat = dat, mspe = 0)
