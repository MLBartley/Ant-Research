## user defined zero inflated poisson distribution

dZIP <- nimbleFunction(
  run = function(x = integer(), lambda = double(), zeroProb = double(), log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if(x != 0) {
      ## return the log probability if log = TRUE
      if(log) return(dpois(x, lambda, log = TRUE) + log(1-zeroProb))
      ## or the probability if log = FALSE
      else return((1-zeroProb) * dpois(x, lambda, log = FALSE))
    }
    ## From here down we know x is 0
    totalProbZero <- zeroProb + (1-zeroProb) * dpois(0, lambda, log = FALSE)
    if(log) return(log(totalProbZero))
    return(totalProbZero)
  })

rZIP <- nimbleFunction(
  run = function(n = integer(), lambda = double(), zeroProb = double()) {
    returnType(integer())
    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if(isStructuralZero) return(0)
    return(rpois(1, lambda))
  })

#register new function distributions

registerDistributions(list(
  dZIP = list(
    BUGSdist = "dZIP(lambda, zeroProb)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = integer()', 'lambda = double()', 'zeroProb = double()')
  )))

## keep latent states constant

library(nimble)
library(coda)

#load latent states (simulated with another NIMBLE model)

load("NIMBLE/model-checking/latentstates_simulated.Rdata")

plot(x.init, type = "l")


# Define model
modelCode <- nimbleCode({

  #priors for lambda values
  lambda_l ~ dgamma(shape = a, rate = b)
  lambda_diff ~ dgamma(shape = c, rate = d)

  lambda_h <- lambda_l + lambda_diff

  #penalty priors for P matrix
  for (i in 1:nStates){
    P[i, ] ~ ddirch(alpha = theta[i, ] )
  }

  ## LATENT STATES

  state[1] ~ dbern(.5)

  for (t in 2:nSecs) {
    state[t] ~ dbern(prob = P[(state[t - 1] +1), 2])
    # changed from [P_prev state, 1] to [P_ps , 2]....unsure if same mistake
    # was in previous simple model code
  }

  #zero inflated probability for dZIP
  p ~ dunif(0,1)

  ## OBSERVED FEEDING INTERACTION DATA

  for (t in 1:nSecs){
    # y[t] ~ dpois(lambda_l + (lambda_diff * (state[t]))) #same as y_l ~pois (lambda_l) and y_h ~ pois(lambda_h * I(state = H))

    y[t] <- y_l[t] + y_diff[t]

    y_l[t] ~ dZIP(lambda_l, zeroProb = p)
    y_diff[t] ~ dZIP(lambda_diff * state[t], zeroProb = p)

    ##need to calculate OSA MSPE within nimble - this way we don't need to monitor states
    # for (t in 1:nSecs){
    # y_hat[t] <- (lambda_l * P[(state[t] + 1), 1] +
    #                lambda_h * P[(state[t] + 1), 2])
    # # mspe_diff[t] <- ((y_hat[t]) - y[t])^2
  }

  # mspe <- 1/nSecs * sum(mspe_diff[1:nSecs])

}
) # end model

penalty <- 300

## define constants, data, and initial values
nStates <- 2
seconds <- 3600 #1 hour - same as lenght of simulated latent states

theta.init <- matrix(c(penalty, 1, 1, penalty), 2, 2)

# data <- list(y = rep(0, seconds))
p.init <- matrix(c(gtools::rdirichlet(1, alpha = theta.init[1, ]),
                   gtools::rdirichlet(1, alpha = theta.init[2, ] )),
                 2, 2)

inits <- list( lambda_l = 0.007,
               lambda_diff = 0.05,
               lambda_h = .007 + .05,
               # state = x.init, #loaded above
               P = p.init,#,
               # y_hat = dat,
               # mspe = 0
               p = .7
)

## Nimble Model

model <- nimbleModel(modelCode,
                     constants <- list(delta_t = 1,
                                       nSecs = seconds,
                                       nStates = nStates,
                                       a = .05, b = .7, c = .05, d = .7,
                                       state = x.init,
                                       theta = matrix(c(penalty, 1,
                                                        1, penalty),
                                                      2, 2)),
                     # data = data,
                     inits = inits,
                     dimensions = list(theta = c(nStates, nStates)))

## Simulate latent states and events

model$simulate(c("y_l", "y_diff"), includeData = TRUE)
model$calculate("y")
simulatedData <- model$y

plot(1:seconds, model$y, type = "l")

model$setData(list(y = simulatedData))

## switching plots

plot(1:seconds, cumsum(model$y), xlab = "Seconds",
     ylab = "Cumulative Interaction Count",
     xlim = c(0, seconds))
states <- as.vector(x.init)
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
spec <- configureMCMC(model, thin = 100)
spec$printSamplers("lambda_l")
spec$printSamplers("lambda_diff")
# spec$printSamplers("y_l")
# spec$pringSamplers("y_diff")
spec$printSamplers("P")

## tried alternative samplers - no mixing improvement

spec$removeSamplers(c("lambda_diff", "lambda_l"))
spec$addSampler(target = c("lambda_diff", "lambda_l"),
                type = "AF_slice")

spec$monitors
# spec$resetMonitors()
# spec$addMonitors(c('lambda_l', 'lambda_diff',
#                    'P', 'mspe'))

## build MCMC algorithm
Rmcmc <- buildMCMC(spec)

## compile model and MCMC
Cmodel <- compileNimble(model)
Cmcmc <- compileNimble(Rmcmc, project = model)

## run MCMC
n_mcmc <- 100000

# Cmcmc$run(n_mcmc, reset = F)
Cmcmc$run(n_mcmc)

## extract samples
samples <- as.matrix(Cmcmc$mvSamples)
write.csv(samples, file =  paste("simple_testmodel_holdlatentstates", "-",
                                 penalty, "-", n_mcmc, ".csv", sep = ""))


# samples <- read.csv(file =  paste("simple_testmodel_holdlatentstates", "-",
# penalty, "-", n_mcmc, ".csv", sep = ""))

coda_samples <- mcmc(samples)
plot(coda_samples)
plot(samples[,'lambda_l'], pch = '.', main = 'lambda low trace plot')
plot(samples[,'lambda_diff'], pch = '.', main = 'lambda trace plot')


effectiveSize(coda_samples[, 1:10])
autocorr(coda_samples)

acf(samples[, "lambda_l"])
acf(samples[, "lambda_diff"])

cor(samples[, c("lambda_l", "lambda_diff")])
