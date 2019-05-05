###############################################################################
## This script aims to use the model samples to check for best (smallest)
## MSPE value
##
## Created: April 15, 2019
## Updated 1:
###############################################################################

#outline:
## load all simple_MCMC*.csv files
## pull mspe estimate from samples
## plot graph

#load data files
path = "NIMBLE/data-mcmc/"
temp = list.files(path = path, pattern="simple*")
mysamples = lapply(paste(path, temp, sep = ""), read.csv)

MSPE_results <-  data.frame(penalty = numeric(0),
                            MSPE = numeric(0),
                            num_iter = numeric(0))


for(i in 1:length(mysamples)) {
  temp.mspe <- mean(mysamples[[i]][, 8])


  penalty <- abs(parse_number(temp[i]))
  iterations <- abs(parse_number(substring(temp[i], 18)))


  MSPE_results <- rbind(MSPE_results, c(penalty, temp.mspe, iterations))
}

colnames(MSPE_results) <- c("penalty", "MSPE", "num_iter")

save(MSPE_results, file = "./NIMBLE/data-prepped/MSPE_simple.Rdata")

