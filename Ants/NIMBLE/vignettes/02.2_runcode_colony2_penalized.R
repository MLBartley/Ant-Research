
# include any code here you don't want to show up in the document,
# e.g. package and dataset loading
library(doMC)
library(methods)  # otherwise new() not being found
library(dplyr)
   # devtools::install_github("nimble-dev/nimble", ref = "avoid-protect-stack-overflow", subdir = "packages/nimble")
library(nimble)
library(tidyr)
library(doParallel)
 library(ggplot2)

# source/load the data
load(here::here("NIMBLE", "data-prepped", "col2preppeddata.Rdata"))
list2env(out,globalenv())

#paths for outputs


#data
dat <- col2_low4_5$starts_persec

antsCode <- nimbleCode({
  #rates of interactions for low and high (low + diff) colonly-level states
  lambda_l ~ dgamma(a, b)
  lambda_diff ~ dgamma(c, d)

  lambda_h <- lambda_l + lambda_diff

  #P is a function of beta, a rate variable for each state
  e.beta[1:num.states] ~ dmnorm(mvnorm.mean[1:num.states],
                                cov = tau[1:num.states, 1:num.states])

  #probablility that the colony switches between states
  for(i in 1:num.states){
    P[i,i] <- 1 / (1 + sum(e.beta[1:num.states]))

    for(j in indx[[i]]){
      P[i, j] <- e.beta[j] / (1 + sum(e.beta[1:num.states]))
    }
  }


  state[1] ~ dbern(.5)

  for (t in 2:nSecs) {
    state[t] ~ dbern(prob = P[(state[t - 1] +1), 1])
  }

  for (t in 1:nSecs){
    y[t] ~ dpois(lambda_l + lambda_diff * (state[t]))
  # }

  ##need to calculate OSA MSPE within nimble - this way we don't need to monitor states
# for (t in 1:nSecs){
  y_hat[t] <- (lambda_l * P[(state[t] + 1), 1] +
                 lambda_h * P[(state[t] + 1), 2])
  mspe_diff[t] <- ((y_hat[t]) - y[t])^2
}

  mspe <- 1/nSecs * sum(mspe_diff[1:nSecs])



}) # end model



# constants: num.states, a, b, c, d, tau, n.Secs
# needs starting values: states, y, P?, lambdas,
## define constants, data, and initial values
num.states = 2
# cov.seconds = 1 #covariates NO
# cov.t = 1
# alpha = 1
# cov.t = cov_col2_low4$cov
# alpha = -0.0001
seconds = length(dat)

indx <- list()
for (i in 1:num.states){
  indx[[i]] <- which(1:num.states != i)
}


data <- list(y = dat)
x.init <- sample(x = c(0, 1), size = seconds, replace = T)
e.beta.init <- exp(c(-4.59, -4.59))

P <- nimArray(NA,
              dim = c(num.states, num.states),
              #cov.sec = 1 if no cov, else = seconds
              init = TRUE)
# tau <- nimMatrix(NA, num.states, num.states, init = TRUE)

# for( t in 1:cov.seconds){
  for(i in 1:num.states){
    P[i,i] <- 1 /
      (1 + sum(e.beta.init[1:num.states]))
    # tau[i, i] <- exp(10)

    for(j in which(1:num.states != i)){
      P[i, j] <- e.beta.init[j] /
        (1 + sum(e.beta.init[1:num.states]))
      # tau[i, j] <- 0
    }
  }
# }

p.init <- P

constants <- list(a = 1, b = 1, c = 1, d = 1,
                  mvnorm.mean = rep(0, num.states),
                  # tau = tau,
                  num.states = num.states,
                  nSecs = seconds,
                  indx = indx)

inits <- list( lambda_l = 0.007, lambda_diff = 0.05, lambda_h = .007 + .05,
               e.beta = e.beta.init,  state = x.init, P = p.init, y_hat = dat, mspe = 0)

# constants.cov <- list(sigma.2 = 100, cov.t = cov.t)
# inits.cov <- list(alpha = alpha)



range = exp(seq(-20, 20, by =  10))
doMC::registerDoMC(cores = 5)


n_mcmc <- 500

mcmc.out <- foreach(i = range,
                    .errorhandling = "remove") %dopar% {

                      # antsCode <- nimbleCode({
                      #   #rates of interactions for low and high (low + diff) colonly-level states
                      #   lambda_l ~ dgamma(a, b)
                      #   lambda_diff ~ dgamma(c, d)
                      #
                      #   lambda_h <- lambda_l + lambda_diff
                      #
                      #   #P is a function of beta, a rate variable for each state
                      #   e.beta[1:num.states] ~ dmnorm(mvnorm.mean[1:num.states],
                      #                                 cov = tau[1:num.states, 1:num.states])
                      #
                      #   #probablility that the colony switches between states
                      #   for(i in 1:num.states){
                      #     P[i,i] <- 1 / (1 + sum(e.beta[1:num.states]))
                      #
                      #     for(j in indx[[i]]){
                      #       P[i, j] <- e.beta[j] / (1 + sum(e.beta[1:num.states]))
                      #     }
                      #   }
                      #
                      #
                      #   state[1] ~ dbern(.5)
                      #
                      #   for (t in 2:nSecs) {
                      #     state[t] ~ dbern(prob = P[(state[t - 1] +1), 1])
                      #   }
                      #
                      #   for (t in 1:nSecs){
                      #     y[t] ~ dpois(lambda_l + lambda_diff * (state[t]))
                      #     # }
                      #
                      #     ##need to calculate OSA MSPE within nimble - this way we don't need to monitor states
                      #     # for (t in 1:nSecs){
                      #     y_hat[t] <- (lambda_l * P[(state[t] + 1), 1] +
                      #                    lambda_h * P[(state[t] + 1), 2])
                      #   }
                      #
                      #
                      #
                      #   mspe[1:nSecs] <- 1/nSecs * sum(((y_hat[1:nSecs]) - y[1:nSecs])^2)
                      # })

                    temp <-  nimbleModel(code = antsCode,
                                  constants = list(a = 1, b = 1, c = 1, d = 1,
                                                   mvnorm.mean = rep(0, num.states),
                                                   num.states = num.states,
                                                   nSecs = seconds,
                                                   indx = indx,
                                                   tau = matrix(c(i, 0, 0, i), 2, 2)),
                                  data = data,
                                  inits = inits,
                                  dimensions = list(P = c(num.states, num.states),
                                                    x.init = seconds,
                                                    e.beta = num.states))

                      spec <- configureMCMC(temp, control = list(reflective = TRUE))
                      spec$resetMonitors()
                      spec$addMonitors(c('lambda_l', 'lambda_diff', 'e.beta')) #NOT monitoring X (states)
                      spec$addMonitors2("mspe")

                      ## build MCMC algorithm
                      Rmcmc <- buildMCMC(spec)
                      ## compile model and MCMC
                      Cmodel <- compileNimble(temp)
                      Cmcmc <- compileNimble(Rmcmc, project = temp, resetFunctions = T)

                      Cmcmc$run(n_mcmc)
                      Cmcmc

                      samples <- as.matrix(Cmcmc$mvSamples)
                      saveRDS(samples, file =  paste("./NIMBLE/data-mcmc/", "pen_MCMC", "-",
                                                  log(i), "-", n_mcmc, ".Rds", sep = ""))

                    }


mcmc.out$summary
mcmc.out$WAI


# df <- data.frame(tempout$samples)
# df_l <- df %>% dplyr::select(lambda_l, lambda_diff, e.beta.1., e.beta.2.) %>% gather(key="parameter", value="value")
#
# ps <- df_l %>% ggplot(aes(x=seq_along(value), y = value)) + geom_line()
#
# pdf(here::here("NIMBLE", "visuals", "tempvis.pdf"))
# ps + facet_wrap(~parameter, scales = "free")
#
#
# p <- ggplot(df_l,aes(value)) + geom_histogram(aes( y= ..density..),bins = 60)
# p + facet_wrap(~parameter, scales = "free")
#
# dev.off()

##testing alt approach

# temp <- nimbleModel(code = antsCode,
#                     constants = list(a = 1, b = 1, c = 1, d = 1,
#                                      mvnorm.mean = rep(0, num.states),
#                                      num.states = num.states,
#                                      nSecs = seconds,
#                                      indx = indx,
#                                      tau = matrix(c(i, 0, 0, i), 2, 2)),
#                     data = data,
#                     inits = inits,
#                     dimensions = list(P = c(num.states, num.states),
#                                       x.init = seconds,
#                                       e.beta = num.states))
#
# spec <- configureMCMC(temp, control = list(reflective = TRUE))
# ## build MCMC algorithm
# Rmcmc <- buildMCMC(spec)
#
#
# ## compile model and MCMC
# Cmodel <- compileNimble(temp)
# Cmcmc <- compileNimble(Rmcmc, project = temp)
#
# set.seed(0)
# Cmcmc$run(500)
Cmcmc$mvSamples

mcmcOut <- as.matrix(Cmcmc$mvSamples)

mcmcOut <- coda::as.mcmc(mcmcOut)

pdf(here::here("NIMBLE", "visuals", "tempvis.pdf"))
coda::traceplot(mcmcOut[,"e.beta[2]"],
          main = "Trace plot for Beta")
dev.off()

# tempout <-
#   nimbleMCMC(code = antsCode,
#              constants = list(a = 1, b = 1, c = 1, d = 1,
#                               mvnorm.mean = rep(0, num.states),
#                               num.states = num.states,
#                               nSecs = seconds,
#                               indx = indx,
#                               tau = matrix(c(i, 0, 0, i), 2, 2)),
#              data = data,
#              inits = inits,
#              # monitors=c("lambda_l", "lambda_diff",
#              #            "e.beta",  "state"),
#              nchains = 1, niter = 500,
#              summary = TRUE, WAIC = TRUE)
