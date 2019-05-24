
library(coda)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(igraph)
library(magrittr)
library(MCMCpack)
library(MCMCvis)
library(nimble)
library(tidyverse)

source("./R/data_prep_functions.R")
source("./R/summary_visual_functions.R")

source("./NIMBLE/vignettes/01_prepdata_colony2.R")


dat <- out$col2_low4_5$queen_starts_persec #created in 01_prepdata
head(dat)
summary(dat)
plot(dat, type = "h")

seconds <- length(dat)


# Define model
modelCode <- nimbleCode({
  lambda_l ~ dgamma(a, b)
  lambda_diff ~ dgamma(c, d)

  lambda_h <- lambda_l + lambda_diff

  for (i in 1:nStates){
    P[i, ] ~ ddirch(alpha = theta[i, ] )
  }

  state[1] ~ dbern(.5)

  for (t in 2:nSecs) {
    state[t] ~ dbern(prob = P[(state[t - 1] +1), 1])
  }

  for (t in 1:nSecs){
    y[t] ~ dpois(lambda_l + lambda_diff * (state[t]))


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
theta.init <- matrix(c(120000, 1, 1, 120000), 2, 2)

constants <- list(delta_t = 1, nSecs = seconds, nStates = nStates,
                  a = 1, b = 1, c = 1, d = 1, theta = theta.init)
data <- list(y = dat)
x.init <- sample(x = c(0, 1), size = seconds, replace = T)
p.init <- matrix(c(.99, 01, .01,.99), 2, 2)

inits <- list( lambda_l = 0.007, lambda_diff = 0.05, lambda_h = .007 + .05,
               state = x.init, P = p.init, y_hat = dat, mspe = 0)

#penalty range
range <- seq(17500, 22000, by =  500)

## create model object
set.seed(0)
n_mcmc <- 10000

mcmc.out <- for(i in 1:length(range)) {

  penalty <- range[i]

Rmodel <- nimbleModel(code = modelCode,
                      constants <- list(delta_t = 1,
                                        nSecs = seconds,
                                        nStates = nStates,
                                        a = 1, b = 1, c = 1, d = 1,
                                        theta = matrix(c(penalty, 1, 1, penalty), 2, 2)),
                      data = data,
                      inits = inits,
                      dimensions = list(theta = c(nStates, nStates)))

## specify MCMC algorithm
spec <- configureMCMC(Rmodel, control = list(reflective = TRUE))
# spec$printSamplers("lambda_l")
# spec$printSamplers("lambda_diff")
# spec$printSamplers("P")

spec$resetMonitors()
spec$addMonitors(c('lambda_l', 'lambda_diff', 'P', 'mspe')) #NOT monitoring X (states)


## build MCMC algorithm
Rmcmc <- buildMCMC(spec)

## compile model and MCMC
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

## run MCMC

Cmcmc$run(n_mcmc)

## extract samples
samples <- as.matrix(Cmcmc$mvSamples)
write.csv(samples, file =  paste("./NIMBLE/data-mcmc/", "simple_MCMC", "-",
                                 penalty, "-", n_mcmc, ".csv", sep = ""))

}



  # means <- apply(samples, 2, mean)
#
# ## plot samples
# df <- data.frame(samples)
# df_l <- df[, 1:6] %>% gather(key="parameter", value="value")
#
# ps <- df_l %>% ggplot(aes(x=seq_along(value), y = value)) + geom_line()
# ps + facet_wrap(~parameter, scales = "free")
#
# p <- ggplot(df_l,aes(value)) + geom_histogram(aes( y= ..density..),bins = 60)
# p + facet_wrap(~parameter, scales = "free")

# states.est <- round(means[-(1:6)], digits = 0)
#
# plot(means[-(1:6)], type = "l")
# plot(states.est, type = "l")
#
# sumvis_low <- sumvis_troph(data = col2_low4, entrance = F, hours = 4, density = "low")
# #Note: cumul.lowqueen is ggplot saved to environment
#
# cumul.lowqueen + geom_point(color = (states.est[cumul.lowqueen$data$start_time] + 1))

