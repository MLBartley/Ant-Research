
library(coda)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(igraph)
library(MCMCpack)
library(MCMCvis)
library(nimble)



dat <- col2_low4_5$queen_starts_persec#made in hsmm_test - need to save but meh
head(dat)
summary(dat)
plot(dat, type = "h")

seconds <- length(dat)


# Define model
modelCode <- nimbleCode({
  lambda_l ~ dgamma(a, b)
  lambda_diff ~ dgamma(c, d)

  lambda_h <- lambda_l + lambda_diff

  gamma_lh ~ T(dnorm(0, sd = tau), 0, 9999)
  gamma_hl ~ T(dnorm(0, sd = tau), 0, 9999)

  P[1, 2] <- gamma_lh * exp(-gamma_lh * delta_t)


  P[1, 1] <- 1 - P[1, 2]

  P[2, 1] <- gamma_hl * exp(-gamma_hl * delta_t)


  P[2, 2] <- 1 - P[2, 1]

  state[1] ~ dbern(.5)


  for (t in 2:nSecs) {
    state[t] ~ dbern(prob = P[(state[t - 1] +1), 1])
  }

  for (t in 1:nSecs){
        y[t] ~ dpois(lambda_l + lambda_diff * (state[t]))
  }

}
) # end model


## define constants, data, and initial values

constants <- list(delta_t = 1, nSecs = seconds,
                  a = 1, b = 1, c = 1, d = 1, tau = exp(-10))
data <- list(y = dat)
x.init <- sample(x = c(0, 1), size = seconds, replace = T)
p.init <- matrix(c(.99, 01, .01,.99), 2, 2)

inits <- list( lambda_l = 0.007, lambda_diff = 0.05, lambda_h = .007 + .05,
               gamma_lh = 0.005, gamma_hl =0.005,  state = x.init, P = p.init)

## create model object
Rmodel <- nimbleModel(code = modelCode, constants = constants, data = data,
                      inits = inits)


## Error in eval(code[[2]], constantsEnv) : object 'state' not found



## specify MCMC algorithm
spec <- configureMCMC(Rmodel, control = list(reflective = TRUE))
spec$printSamplers("lambda_l")
spec$printSamplers("lambda_diff")
spec$printSamplers("gamma_lh")

## build MCMC algorithm
Rmcmc <- buildMCMC(spec)

## compile model and MCMC
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

## run MCMC
set.seed(0)
Cmcmc$run(10000)

## extract samples
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
