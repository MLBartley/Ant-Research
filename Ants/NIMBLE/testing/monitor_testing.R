## does adding state monitor really chnage the MSPE/results??

library(magrittr)
# library(coda, lib.loc = "/usr/lib/R/site-library") #need for EH machine
library(coda)
library(dplyr)

load("./NIMBLE/data-prepped/MSPE_penalized.Rdata")



MSPE_results_summary <- MSPE_results %>%
  group_by(penalty) %>%
  summarise(MSPE = mean(MSPE))


best <- which(MSPE_results_summary[,2] == min(MSPE_results_summary[,2]))
n_mcmc <-50001

penalty <- exp(MSPE_results_summary[best,1])

nStates <- 2

source("./NIMBLE/vignettes/01.03_prepdata_penModel.R")

####RErun MCMC



Rmodel <- nimbleModel(code = antsCode,
                      constants <- list(a = 1, b = 1, c = 1, d = 1,
                                        mvnorm.mean = rep(0, num.states),
                                        num.states = num.states,
                                        nSecs = seconds,
                                        indx = indx,
                                        tau = matrix(c(penalty, 0, 0, penalty), 2)),
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
                   'e.beta',
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
samples_all <- as.matrix(Cmcmc$mvSamples)
write.csv(samples_all, file =  paste("./NIMBLE/data-mcmc/", "spec_test_all", "-",
                                 log(penalty), "-", n_mcmc, ".csv", sep = ""))



## reset spec

## specify MCMC algorithm
spec <- configureMCMC(Rmodel, control = list(reflective = TRUE))

spec$resetMonitors()
spec$addMonitors(c('lambda_l',
                   'lambda_diff',
                   'e.beta',
                   'P',
                   # 'state',
                   'mspe'))


## build MCMC algorithm
Rmcmc <- buildMCMC(spec)

## compile model and MCMC
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

## run MCMC

Cmcmc$run(n_mcmc)

## extract samples
samples_nostates <- as.matrix(Cmcmc$mvSamples)
write.csv(samples_nostates, file =  paste("./NIMBLE/data-mcmc/", "spec_test_nostates", "-",
                                 log(penalty), "-", n_mcmc, ".csv", sep = ""))


## specify MCMC algorithm
spec <- configureMCMC(Rmodel, control = list(reflective = TRUE))

# spec$resetMonitors()
# spec$addMonitors(c('lambda_l',
#                    'lambda_diff',
#                    'e.beta',
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
samples_noreset <- as.matrix(Cmcmc$mvSamples)
write.csv(samples_noreset, file =  paste("./NIMBLE/data-mcmc/", "spec_test_all_noreset", "-",
                                 log(penalty), "-", n_mcmc, ".csv", sep = ""))



## compare mean MSPE values

mspe_all <- mean(samples_all[, 9])
plot(samples_all[, 9], type= "l")
# mspe_norest <- mean(samples_noreset[, ])

mspe_nostates <- mean(samples_nostates[, 9])
plot(samples_nostates[, 9], type= "l")


coda_samples_all <- mcmc(samples_all[, 1:10])
plot(coda_samples_all)

coda_samples_noreset <- mcmc(samples_noreset[, 1:5])
plot(coda_samples_noreset)

coda_samples_nostates <- mcmc(samples_nostates[, 1:9])
plot(coda_samples_nostates)


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
state_samples <- samples[, -c(1:10)]

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

plot(state_samples[, 50], type = "l")

