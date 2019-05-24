
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
range <- seq(35000, 120000, by =  5000)

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

