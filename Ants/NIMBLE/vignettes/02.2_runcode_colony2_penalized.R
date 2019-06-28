
# include any code here you don't want to show up in the document,
# e.g. package and dataset loading
library(doMC)
library(methods)  # otherwise new() not being found
library(dplyr)
   # devtools::install_github("nimble-dev/nimble", ref = "avoid-protect-stack-overflow", subdir = "packages/nimble")
library(nimble, lib.loc = "/usr/lib/R/site-library")
library(tidyr)
# library(doParallel)
 library(ggplot2)

source("./NIMBLE/vignettes/01.03_prepdata_penModel.R")

range = exp(seq(-4, 10, by =  1))

doParallel::registerDoParallel(cores = 5)

n_mcmc <- 60000


mcmc.out <- foreach(i = 1:length(range)) %dopar% {

  # mcmc.out <- for(i in 1:1) {
  # mcmc.out <- for(i in 1:length(range)) {

 penalty <- range[i]

                    model <-  nimbleModel(code = antsCode,
                                          constants <- list(a = 1, b = 1, c = 1, d = 1,
                                                            mvnorm.mean = rep(0, num.states),
                                                            num.states = num.states,
                                                            nSecs = seconds,
                                                            indx = indx,
                                                            tau = matrix(c(penalty, 0, 0, penalty), 2)),
                                  data = data,
                                  inits = inits,
                                  dimensions = list(P = c(num.states, num.states),
                                                    x.init = seconds,
                                                    e.beta = num.states,
                                                    tau = c(num.states, num.states)))

                      spec <- configureMCMC(model, control = list(reflective = TRUE))
                      spec$resetMonitors()
                      spec$addMonitors(c('lambda_l', 'lambda_diff', 'e.beta', 'mspe')) #NOT monitoring X (states)

                      ## build MCMC algorithm
                      Rmcmc <- buildMCMC(spec)
                      ## compile model and MCMC
                      Cmodel <- compileNimble(model, resetFunctions = T, showCompilerOutput = T)
                      Cmcmc<- compileNimble(Rmcmc, project = model)

                      Cmcmc$run(n_mcmc)

                      samples <- as.matrix(Cmcmc$mvSamples)
                      write.csv(samples, file =  paste("./NIMBLE/data-mcmc/", "pen_MCMC", "-",
                                                  log(penalty), "-", n_mcmc, ".csv", sep = ""))

                    }

