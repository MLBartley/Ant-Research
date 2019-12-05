###############################################################################
## This script aims to chose the best MSPE value, rerun MCMC to obtain latent
## state samples and to create the switching plots
##
## Created: April 15, 2019
## Updated 1:
###############################################################################
library(magrittr)
library(dplyr)
library(coda)

load("./NIMBLE/data-prepped/MSPE_simple.Rdata")



MSPE_results_summary <- MSPE_results %>%
  group_by(penalty) %>%
  summarise(MSPE = mean(MSPE))


best <- which(MSPE_results_summary[,2] == min(MSPE_results_summary[,2]))

n_mcmc <- 20003

penalty <- MSPE_results_summary[best,1]
 penalty <- 265 #slight raise in mspe after here (next try above line)
 penalty <- 1 #no penalty


source("./NIMBLE/vignettes/01.02_prepdata_simpleModel.R")

 ##fiddling around
 # penalty <- 1 #ie no penalty
 inits$lambda_diff <- .1

 inits$lambda_l <- .01

 ####RErun MCMC

Rmodel <- nimbleModel(code = modelCode,
                      constants <- list(delta_t = 1,
                                        nSecs = seconds,
                                        nStates = nStates,
                                        # a = .005, b = .7, c = .05, d = .7,
                                        a = 1, b = 1, c = 1, d = 1,
                                        theta = matrix(c(penalty, 1, 1, penalty), 2, 2)),
                      data = dat,
                      inits = inits,
                      dimensions = list(theta = c(nStates, nStates)))

Rmodel$setData(list(y = data$y))


## specify MCMC algorithm
spec <- configureMCMC(Rmodel)

spec$resetMonitors()
spec$addMonitors(c('lambda_l',
                   'lambda_diff',
                   'P',
                   'state',
                   'mspe'))

# spec$addMonitors(c('lambda',
#                    'P',
#                    'state',
#                    'mspe'))

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


samples <- read.csv(file =  paste("./NIMBLE/data-mcmc/", "simple_MCMC", "-",
 penalty, "-", n_mcmc, ".csv", sep = ""))



coda_samples <- mcmc(samples[, 1:10])
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
state_samples <- samples[, -c(1:8)]

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
plot(round(states_est), type = "l")

plot(samples[, 50], type = "l")
