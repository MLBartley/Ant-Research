##
##
## Model Checking - simple
##
##

# outline

  ## Simulate data
  ## model with fixed latent states, X_t
  ## model with fixed event rates, lambda_{X_t}
  ##


library(nimble)

# Define model
modelCode <- nimbleCode({

  #priors for lambda values
  lambda_l ~ dgamma(shape = a, rate = b)
  lambda_diff ~ dgamma(shape = c, rate = d)

  lambda_h <- lambda_l + lambda_diff

  for (i in 1:nStates){
    P[i, ] ~ ddirch(alpha = theta[i, ] )
  }

  state[1] ~ dbern(.5)

  for (t in 2:nSecs) {
    state[t] ~ dbern(prob = P[(state[t - 1] +1), 2])
    # changed from [P_prev state, 1] to [P_ps , 2]....unsure if same mistake
    # was in previous simple model code
  }

  for (t in 1:nSecs){
    # y[t] ~ dpois(lambda_l + (lambda_diff * (state[t]))) #same as y_l ~pois (lambda_l) and y_h ~ pois(lambda_h * I(state = H))

    y[t] <- y_l[t] + y_diff[t] * state[t]

    y_l[t] ~ dpois(lambda_l)
    y_diff[t] ~ dpois(lambda_diff)

    ##need to calculate OSA MSPE within nimble - this way we don't need to monitor states
    # for (t in 1:nSecs){
    # y_hat[t] <- (lambda_l * P[(state[t] + 1), 1] + lambda_h * P[(state[t] + 1), 2])
    # mspe_diff[t] <- ((y_hat[t]) - y[t])^2
  }

  # mspe <- 1/nSecs * sum(mspe_diff[1:nSecs])

}
) # end model

penalty <- 300

## define constants, data, and initial values
nStates <- 2
seconds <- 3600

theta.init <- matrix(c(penalty, 1, 1, penalty), 2, 2)

data <- list(y = rep(0, seconds))
x.init <- sample(x = c(0, 1), size = seconds, replace = T)
p.init <- matrix(c(gtools::rdirichlet(1, alpha = theta.init[1, ]),
                   gtools::rdirichlet(1, alpha = theta.init[2, ] )),
                 2, 2)

inits <- list( lambda_l = 0.007,
               lambda_diff = 0.05,
               lambda_h = .007 + .05,
               state = x.init,
               P = p.init,
               y_hat = data,
               mspe = 0)

constants <- list(nSecs = seconds,
     nStates = nStates,
     a = .005, b = .7, c = .05, d = .7,
     theta = matrix(c(penalty, 1, 1, penalty), 2, 2))

## Nimble Model

model <- nimbleModel(modelCode,
                     constants = constants,
                     # data = data,
                     inits = inits,
                     dimensions = list(theta = c(nStates, nStates)))

## Simulate latent states

model$simulate("state")

plot(model$state, type = "l")

# model$simulate("lifted_lambda_diff_times_state_oBt_cB_L12")

## Simulate events

model$simulate(c("y_l", "y_diff"))
model$calculate("y")

simulatedData <- model$y

model$setData(list(y = simulatedData))

plot(1:seconds, model$y, type = "l")


## switching plots

plot(1:seconds, cumsum(model$y), xlab = "Seconds",
     ylab = "Cumulative Interaction Count",
     xlim = c(0, seconds))
states <- as.vector(model$state)
rr <- rle(states)
embedded.chain <- rr$values
cs <- c(0, cumsum(rr$lengths)) - 1
cols <- c("#bc535644", "#538bbc44")
for (j in 1:length(embedded.chain)) {
  rect(cs[j], 0, cs[j + 1], seconds,
       col = cols[embedded.chain[j] + 1]
      )
}


## specify MCMC algorithm
spec <- configureMCMC(model)
 spec$printSamplers("lambda_l")
 spec$printSamplers("lambda_diff")
 spec$printSamplers("P")
 # spec$printSamplers("state")
 spec$printSamplers("y_diff")[1:10]

 spec$removeSamplers(c("y_l"))
 spec$addSampler(target = c("y_l"),
                 type = "AF_slice",
                 control = list(adaptInterval = 100))

 spec$removeSamplers(c("y_diff"))
 spec$addSampler(target = c("y_diff"),
                 type = "AF_slice",
                 control = list(adaptInterval = 100))

 spec$monitors
 # spec$resetMonitors()
 # spec$addMonitors(c('lambda_l', 'lambda_diff'))

 spec$addMonitors2(c('y_l', 'y_diff'))

## build MCMC algorithm
Rmcmc <- buildMCMC(spec)

## compile model and MCMC
Cmodel <- compileNimble(model)
Cmcmc <- compileNimble(Rmcmc, project = model)

## run MCMC
n_mcmc <- 10000

Cmcmc$run(n_mcmc)

## extract samples
samples <- as.matrix(Cmcmc$mvSamples)
samples2 <- as.matrix(Cmcmc$mvSamples2)
write.csv(samples, file =  paste("./NIMBLE/data-mcmc/", "simple_testmodel", "-",
                                 penalty, "-", n_mcmc, ".csv", sep = ""))





coda_samples <- mcmc(samples)
coda_samples2 <- mcmc(samples2)

plot(coda_samples)
plot(coda_samples2[, 3600:3610])

effectiveSize(coda_samples[, 1:16])
summary(effectiveSize(coda_samples2))

autocorr(coda_samples[, 1:6])


acf(samples[, "lambda_l"])
acf(samples[, "lambda_diff"])

cor(samples[, c("lambda_l", "lambda_diff")])

y_low_est <- apply(samples2[,1:3600], 2, mean)
plot(y_low_est, type="l")

y_diff_est <- apply(samples2[,3601:7200], 2, mean)
plot(y_diff_est, type="l")

#state estimates
state_samples <- samples[, -c(1:6)]

states_est <-apply(state_samples, 2, mean)

plot(round(states_est), type = "l")
plot(states_est, type = "l")
plot(model$state, type = "l") #truth


