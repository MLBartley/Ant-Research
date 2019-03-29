
# include any code here you don't want to show up in the document,
# e.g. package and dataset loading
library(doMC)
library(methods)  # otherwise new() not being found
# install_github("nimble-dev/nimble", ref = "avoid-protect-stack-overflow", subdir = "packages/nimble")
library(nimble)
library(tidyr)

# source/load the data
load("../data-prepped/col2preppeddata.Rdata")
list2env(out,globalenv())
#paths for outputs



antsCode <- nimbleCode({
  #rates of interactions for low and high (low + diff) colonly-level states
  lambda_l ~ dgamma(a, b)
  lambda_diff ~ dgamma(c, d)

  lambda_h <- lambda_l + lambda_diff

  #P is a function of beta, a rate variable for each state
  e.beta[1:num.states] ~ dmnorm(mvnorm.mean[1:num.states],
                                tau[1:num.states, 1:num.states])

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
  }

}
) # end model



# constants: num.states, a, b, c, d, tau, n.Secs
# needs starting values: states, y, P?, lambdas,
## define constants, data, and initial values
num.states = 2
cov.seconds = seconds #covariates YES
cov.t = cov_col2_low4$cov
alpha = -0.0001

indx <- list()
for (i in 1:num.states){
  indx[[i]] <- which(1:num.states != i)
}


data <- list(y = dat)
x.init <- sample(x = c(0, 1), size = seconds, replace = T)
e.beta.init <- exp(c(-4.59, -4.59))

P <- nimArray(NA,
              dim = c(num.states, num.states, cov.seconds),
              #cov.sec = 1 if no cov, else = seconds
              init = TRUE)
tau <- nimMatrix(NA, num.states, num.states, init = TRUE)

for( t in 1:cov.seconds){
  for(i in 1:num.states){
    P[i,i,t ] <- 1 /
      (1 + sum(exp(log(e.beta.init[1:num.states]) *cov.t[t] * alpha)))
    tau[i, i] <- exp(10)

    for(j in which(1:num.states != i)){
      P[i, j, t] <- exp(log(e.beta.init[j]) * cov.t[t] * alpha) /
        (1 + sum(exp(log(e.beta.init[1:num.states]) * cov.t[t] * alpha)))
      tau[i, j] <- 0
    }
  }
}

p.init <- P

constants <- list(a = 1, b = 1, c = 1, d = 1,
                  mvnorm.mean = rep(0, num.states),
                  # tau = tau,
                  num.states = num.states,
                  nSecs = seconds,
                  indx = indx)

inits <- list( lambda_l = 0.007, lambda_diff = 0.05, lambda_h = .007 + .05,
               e.beta = e.beta.init,  state = x.init, P = p.init)

constants.cov <- list(sigma.2 = 100, cov.t = cov.t)
inits.cov <- list(alpha = alpha)



range = exp(seq(-40, 20, by =  2))

mcmc.out <- foreach(i = range,
                    .errorhandling = "remove") %dopar%
  nimbleMCMC(code = antsCode,
             constants = c(constants, tau = i),
             data = data, inits = inits,
             monitors=c("lambda_l", "lambda_diff",
                        "e.beta",  "state"),
             nchains = 1, niter = 1000,
             summary = TRUE, WAIC = TRUE)

mcmc.out$summary
mcmc.out$WAI


df <- data.frame(mcmc.out$samples)
df_l <- df %>% select(lambda_l, lambda_diff, e.beta.1., e.beta.2.) %>% gather(key="parameter", value="value")

ps <- df_l %>% ggplot(aes(x=seq_along(value), y = value)) + geom_line()
ps + facet_wrap(~parameter, scales = "free")

p <- ggplot(df_l,aes(value)) + geom_histogram(aes( y= ..density..),bins = 60)
p + facet_wrap(~parameter, scales = "free")




