###############################################################################
## This script aims to chose the best MSPE value, rerun MCMC to obtain latent
## state samples and to create the switching plots
##
## Created: April 15, 2019
## Updated 1:
###############################################################################
library(coda)

load("./NIMBLE/data-prepped/MSPE_simple.Rdata")


MSPE_results_summary <- MSPE_results %>%
  group_by(penalty) %>%
  summarise(MSPE = mean(MSPE))


best <- which(MSPE_results_summary[,2] == min(MSPE_results_summary[,2]))
n_mcmc <- 50001

penalty <- MSPE_results_summary[best,1]

source("./NIMBLE/vignettes/01.02_prepdata_simpleModel.R")

####RErun MCMC




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
  rect(cs[j], 0, cs[j + 1], int.num, col = cols[embedded.chain[j]],
       density = 0)
}
points(start, 1:int.num, xlab = "Seconds", ylab = "Cumulative Interaction Count",
       xlim = c(0, maxtime))

