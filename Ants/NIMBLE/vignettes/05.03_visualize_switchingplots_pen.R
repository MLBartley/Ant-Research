###############################################################################
## This script aims to chose the best MSPE value, rerun MCMC to obtain latent
## state samples and to create the switching plots
##
## Created: May 27, 2019
## Updated 1:
###############################################################################
library(magrittr)
library(coda, lib.loc = "/usr/lib/R/site-library") #need for EH machine
# library(coda)
library(dplyr)

load("./NIMBLE/data-prepped/MSPE_penalized.Rdata")



MSPE_results_summary <- MSPE_results %>%
  group_by(penalty) %>%
  summarise(MSPE = mean(MSPE))


best <- which(MSPE_results_summary[,2] == min(MSPE_results_summary[,2]))
n_mcmc <- 60001

penalty <- exp(MSPE_results_summary[best,1])

nStates <- 2

source("./NIMBLE/vignettes/01.03_prepdata_penModel.R")

####RErun MCMC



Rmodel <- nimbleModel(code = antsCode,
                      constants <- constants,
                      data = dat,
                      inits = inits,
                      dimensions = list(P = c(num.states, num.states),
                                        x.init = seconds,
                                        e.beta = num.states,
                                        tau = c(num.states, num.states)))

## specify MCMC algorithm
spec <- configureMCMC(Rmodel, control = list(reflective = TRUE))

spec$resetMonitors()
spec$addMonitors(c('lambda_l',
                   'lambda_diff',
                   'P',
                   'state',
                   'mspe'))


## build MCMC algorithm
Rmcmc <- buildMCMC(spec)

## compile model and MCMC
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

## run MCMC

Cmcmc$run(n_mcmc)

## extract samples
samples <- as.matrix(Cmcmc$mvSamples[, -(1:1000)])
write.csv(samples, file =  paste("./NIMBLE/data-mcmc/", "pen_MCMC", "-",
                                 penalty, "-", n_mcmc, ".csv", sep = ""))

# samples <- read.csv(file =  paste("./NIMBLE/data-mcmc/", "simple_MCMC", "-",
# penalty, "-", n_mcmc, ".csv", sep = ""))



coda_samples <- mcmc(samples[, 1:20])
plot(coda_samples)

##load in ant data

ant_file = col2_low4_5$data
chamber = "queen"
hours <- 4

# Low Density - 4 Hours

location <- ant_file$Location
start <- ant_file$start_time
start = start[which(location == 1)] #queen chamber only
int.num <- length(start)
maxtime <- hours * 60 * 60


#state estimates
state_samples <- samples[, -c(1:7)]

states_est <-apply(state_samples, 2, mean)


plot(start, 1:int.num, xlab = "Seconds",
     ylab = "Cumulative Interaction Count",
     xlim = c(0, maxtime))
states <- states_est
rr <- rle(states)
rr$values <- round(rr$values, digits = 0)
embedded.chain <- rr$values
cs <- c(0, cumsum(rr$lengths)) - 1
cols <- c("#bc535644", "#538bbc44")
for (j in 1:length(embedded.chain)) {
  rect(cs[j], 0, cs[j + 1], int.num, col = cols[embedded.chain[j] + 1],
       density = 0)
}
points(start, 1:int.num, xlab = "Seconds", ylab = "Cumulative Interaction Count",
       xlim = c(0, maxtime))

##peeky peek

plot(states_est, type = "l")

