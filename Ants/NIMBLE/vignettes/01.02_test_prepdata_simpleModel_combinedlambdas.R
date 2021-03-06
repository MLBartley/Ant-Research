
###############################################################################
## This script aims to initialize all simple model details
##
## Created: May 18, 2019
## Updated 1:
###############################################################################
library(nimble)

source("./NIMBLE/vignettes/01.01_prepdata_colony2.R")

dat <- out$col2_low4_5$queen_starts_persec #created in 01_prepdata

seconds <- length(dat)

# Define model
modelCode <- nimbleCode({
  # lambda_l ~ dgamma(a, b)
  # lambda_diff ~ dgamma(c, d)
  #
  # lambda_h <- lambda_l + lambda_diff

  for (i in 1:nStates){
    lambda[1] ~ dgamma(1, 1)
    P[i, ] ~ ddirch(alpha = theta[i, ] )
  }

  state[1] ~ dbern(.5)

  for (t in 2:nSecs) {
    state[t] ~ dbern(prob = P[(state[t - 1] +1), 2])
    # changed from [P_prev state, 1] to [P_ps , 2]....unsure if same mistake
    # was in previous simple model code
  }

  for (t in 1:nSecs){
    y[t] ~ dpois(lambda[state[t] + 1])


    ##need to calculate OSA MSPE within nimble - this way we don't need to monitor states
    # for (t in 1:nSecs){
    y_hat[t] <- (lambda[1] * P[(state[t] + 1), 1] +
                   lambda[2] * P[(state[t] + 1), 2])
    mspe_diff[t] <- ((y_hat[t]) - y[t])^2
  }

  mspe <- 1/nSecs * sum(mspe_diff[1:nSecs])

}
) # end model


## define constants, data, and initial values
nStates <- 2
theta.init <- matrix(c(100, 1, 1, 100), 2, 2)

# constants <- list(delta_t = 1, nSecs = seconds, nStates = nStates,
                  # a = 1, b = 1, c = 1, d = 1, theta = theta.init)
data <- list(y = dat)
x.init <- sample(x = c(0, 1), size = seconds, replace = T)
p.init <- matrix(c(gtools::rdirichlet(1, alpha = theta.init[1, ]),
                                     gtools::rdirichlet(1, alpha = theta.init[2, ] )),
                   2, 2)

# inits <- list( lambda_l = 0.007, lambda_diff = 0.05, lambda_h = .007 + .05,
inits <- list(lambda = c(.007, .057),
              state = x.init, P = p.init, y_hat = dat, mspe = 0)
