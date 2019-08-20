
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

  #priors for lambda values
  lambda_l ~ dgamma(shape = a, rate = b)
  lambda_diff ~ dgamma(shape = c, rate = d)

  lambda_h <- lambda_l + lambda_diff

  for (i in 1:nStates){
    P[i, ] ~ ddirch(alpha = theta[i, ] )
  }

  state[1] ~ dbern(.5)

  for (t in 2:nSecs) {
    state[t] ~ dbern(prob = P[(state[t - 1] +1), 2])
    # changed from [P_prev state, 1] to [P_ps , 2]....unsure if same mistake
    # was in previous simple model code
  }

  ## OBSERVED FEEDING INTERACTION DATA
  prob[1] <- lambda_l
  prob[2] <- lambda_diff

  for (t in 1:nSecs){

    # y_l[t] ~ dZIP(lambda_l, zeroProb = p)
    # y_l[t] ~ dpois(lambda_l)
    # y_diff[t] ~ dpois(lambda_diff)
    #
    #  y[t] <- y_l[t] + y_diff[t] * state[t]

    y[t] ~ dpois(lambda_l + (lambda_diff * (state[t]))) #same as y_l ~pois (lambda_l) and y_h ~ pois(lambda_h * I(state = H))

    split[t, 1:2] ~ dmultinom(size = y[t],
                              prob = prob[1:2])

    y_l[t] <- (state[t] == 0) * y[t] + (state[t] != 0) * split[t, 1] * (y[t] != 0)
    y_diff[t] <- (state[t] != 0) * split[t, 2] * (y[t] != 0)

    ##need to calculate OSA MSPE within nimble - this way we don't need to monitor states
    # for (t in 1:nSecs){
    y_hat[t] <- (lambda_l * P[(state[t] + 1), 1] +
                   lambda_h * P[(state[t] + 1), 2])
    mspe_diff[t] <- ((y_hat[t]) - y[t])^2
  }

  mspe <- 1/nSecs * sum(mspe_diff[1:nSecs])

}
) # end model


## define constants, data, and initial values
nStates <- 2
theta.init <- matrix(c(100, 1, 1, 100), 2, 2)

constants <- list(delta_t = 1,
                  nSecs = seconds,
                  nStates = nStates,
                  a = 1, b = 1, c = 1, d = 1, theta = theta.init)
data <- list(y = dat)
x.init <- sample(x = c(0, 1), size = seconds, replace = T)
p.init <- matrix(c(gtools::rdirichlet(1, alpha = theta.init[1, ]),
                                     gtools::rdirichlet(1, alpha = theta.init[2, ] )),
                   2, 2)

inits <- list( lambda_l = 0.007, lambda_diff = 0.05, lambda_h = .007 + .05,
               state = x.init, P = p.init, y_hat = dat, mspe = 0)
