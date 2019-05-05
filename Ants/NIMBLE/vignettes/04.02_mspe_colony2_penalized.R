###############################################################################
## This script aims to use the model samples to check for best (smallest)
## MSPE value
##
## Created: May 5, 2019
## Updated 1:
###############################################################################

#outline:
## load all simple_MCMC*.csv files
## pull mspe estimate from samples
## plot graph

#load data files
path = "NIMBLE/data-mcmc/"
temp = list.files(path = path, pattern="pen*")
mysamples = lapply(paste(path, temp, sep = ""), read.csv)

MSPE_results <-  data.frame(penalty = numeric(0),
                            MSPE = numeric(0),
                            num_iter = numeric(0))


for(i in 1:length(mysamples)) {
  temp.mspe <- mean(mysamples[[i]][, 6])


  penalty <- readr::parse_number(substring(temp[i], 10))
  iterations <- abs(readr::parse_number(stringr::str_sub(temp[i], start= -9)))


  MSPE_results <- rbind(MSPE_results, c(penalty, temp.mspe, iterations))
}

colnames(MSPE_results) <- c("penalty", "MSPE", "num_iter")

save(MSPE_results, file = "./NIMBLE/data-prepped/MSPE_penalized.Rdata")

