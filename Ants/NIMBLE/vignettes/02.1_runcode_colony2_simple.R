#
# library(coda)
# library(data.table)
# library(ggplot2)
# library(grid)
# library(gridExtra)
# library(igraph)
# library(magrittr)
# library(MCMCpack)
# library(MCMCvis)
# library(nimble)
# library(tidyverse)
library(doParallel)

source("./R/data_prep_functions.R")
source("./R/summary_visual_functions.R")

source("./NIMBLE/vignettes/01.01_prepdata_colony2.R")
source("./NIMBLE/vignettes/01.02_prepdata_simpleModel.R")


#penalty range
range <- seq(100, 200, by =  4)


## create model object
n_mcmc <- 10010

doParallel::registerDoParallel(cores = 5)


mcmc.out <- foreach(i = 1:length(range)) %dopar% {

# mcmc.out <- for(i in 1:length(range)) {

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

## check what needs to be initialized (if anything)
Rmodel$initializeInfo()

#check data is set
Rmodel$y[1:10]
Rmodel$isData('y')
Rmodel$setData(list(y = data$y))


## specify MCMC algorithm
## specify MCMC algorithm
spec <- configureMCMC(Rmodel)
spec$printSamplers("lambda_l")
spec$printSamplers("lambda_diff")
spec$printSamplers("P")

# spec$removeSamplers(c("y_l"))
# spec$addSampler(target = c("y_l"),
#                 type = "RW_block",
#                 control = list(adaptInterval = 100))
#
# spec$removeSamplers(c("y_diff"))
# spec$addSampler(target = c("y_diff"),
#                 type = "RW_block",
#                 control = list(adaptInterval = 100))

spec$monitors
 spec$resetMonitors()
#
 spec$addMonitors(c('lambda_l', "lambda_diff",
                   "P", "mspe"))
# spec$addMonitors2(c("y_l", "y_diff", 'mspe'))


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
