
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

source("./NIMBLE/vignettes/01.01_prepdata_colony2.R")
source("./NIMBLE/vignettes/01.02_prepdata_simpleModel.R")


#penalty range
range <- seq(1000,15000, by =  2500) #next 16000 to 30000, then 33500 to 48500
range <- 300

## create model object
set.seed(0)
n_mcmc <- 80000

mcmc.out <- for(i in 1:length(range)) {

  penalty <- range[i]

Rmodel <- nimbleModel(code = modelCode,
                      constants <- list(delta_t = 1,
                                        nSecs = seconds,
                                        nStates = nStates,
                                        a = .005, b = .7, c = .05, d = .7,
                                        theta = matrix(c(penalty, 1, 1, penalty), 2, 2)),
                      data = data,
                      inits = inits,
                      dimensions = list(theta = c(nStates, nStates)))


## specify MCMC algorithm
spec <- configureMCMC(Rmodel, control = list(reflective = TRUE))
# spec$printSamplers("lambda_l")
# spec$printSamplers("lambda_diff")
# spec$printSamplers("P")

spec$monitors
spec$addMonitors2(c('y_l', 'y_diff', "mspe"))


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
