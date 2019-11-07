###############################################################################
## This script aims to produce plots relevant to the "simple" model
##
## Created: September 19, 2019
## Updated 1:
###############################################################################

#outline:
## switching plot of unpenalized "simple" model (dirichlet priors)
## switching plot of penalized "simple" model

##load needed ant data for plots



## load in unpenalized "simpe" model
samples <- read.csv(file =  paste("./NIMBLE/data-mcmc/simple_MCMC", "-",
                                  1, "-", 20002, ".csv", sep = ""))


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
