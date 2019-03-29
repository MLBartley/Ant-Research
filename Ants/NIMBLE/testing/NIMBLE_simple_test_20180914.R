
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


#data directory
dat.dir <- "../data-raw/"
#output directory
out.dir <- "../output/"

#column names - helps when only pulling in those columns, no extra

col_names <- c("Location", "Ant_ID", "Ant_ID_partner", "start_time", "end_time")

col2_low4 <- read.csv("./data-raw/Colony2_trophallaxis_low_density_4hr.csv")

#removes any extra columns, rows, and adds column names - depends on col_names being correct length
col2_low4 <- col2_low4[, 1:length(col_names)]
col2_low4 <- col2_low4 %>%
  tidyr::drop_na()
colnames(col2_low4) <- col_names


col2_low4_5 <- prep_troph_data(col2_low4, hours = 4, delta_t =  5)


dat <- col2_low4_5$queen_starts_persec
head(dat)
summary(dat)
plot(dat, type = "h")

seconds <- length(dat)


# Define model
modelCode <- nimbleCode({
  lambda_l ~ dgamma(a, b)
  lambda_diff ~ dgamma(c, d)

  lambda_h <- lambda_l + lambda_diff

  for (i in 1:nStates){
    P[i, ] ~ ddirch(alpha = theta[i, ] )
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


## define constants, data, and initial values
nStates <- 2
theta.init <- matrix(c(120000, 1, 1, 120000), 2, 2)

constants <- list(delta_t = 1, nSecs = seconds, nStates = nStates,
                  a = 1, b = 1, c = 1, d = 1, theta = theta.init)
data <- list(y = dat)
x.init <- sample(x = c(0, 1), size = seconds, replace = T)
p.init <- matrix(c(.99, 01, .01,.99), 2, 2)

inits <- list( lambda_l = 0.007, lambda_diff = 0.05, lambda_h = .007 + .05,
              state = x.init, P = p.init)

## create model object
Rmodel <- nimbleModel(code = modelCode, constants = constants, data = data,
                      inits = inits, dimensions = list(theta = c(nStates, nStates)))


## Error in eval(code[[2]], constantsEnv) : object 'state' not found



## specify MCMC algorithm
spec <- configureMCMC(Rmodel, control = list(reflective = TRUE))
spec$printSamplers("lambda_l")
spec$printSamplers("lambda_diff")
spec$printSamplers("P")
spec$printSamplers("y")

## build MCMC algorithm
Rmcmc <- buildMCMC(spec)

## compile model and MCMC
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

## run MCMC
set.seed(0)
Cmcmc$run(10000)

## extract samples
samples <- as.matrix(Cmcmc$mvSamples)
means <- apply(samples, 2, mean)

## plot samples
df <- data.frame(samples)
df_l <- df[, 1:6] %>% gather(key="parameter", value="value")

ps <- df_l %>% ggplot(aes(x=seq_along(value), y = value)) + geom_line()
ps + facet_wrap(~parameter, scales = "free")

p <- ggplot(df_l,aes(value)) + geom_histogram(aes( y= ..density..),bins = 60)
p + facet_wrap(~parameter, scales = "free")

states.est <- round(means[-(1:6)], digits = 0)

plot(means[-(1:6)], type = "l")
plot(states.est, type = "l")

sumvis_low <- sumvis_troph(data = col2_low4, entrance = F, hours = 4, density = "low")
#Note: cumul.lowqueen is ggplot saved to environment

cumul.lowqueen + geom_point(color = (states.est[cumul.lowqueen$data$start_time] + 1))

