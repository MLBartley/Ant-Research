###############################################################################
## This script aims to use the model samples to check for convergance
##
## Created: April 5, 2019
## Updated 1:
###############################################################################

#outline:
  ## load all pen_MCMC*.csv files
  ## coda package check trace + histograms
  ##

library(coda)

#load data files
path = "NIMBLE/data-mcmc/"
temp = list.files(path = path, pattern="*.csv")
mysamples = lapply(paste(path, temp, sep = ""), read.csv)


##checking
coda_samples <- list()

for(i in 1:length(mysamples)){
  coda_samples[[i]] <- mcmc(mysamples[[i]])
  pdf(file = here::here("NIMBLE", "visuals", paste("codaplot_", temp[[i]], ".pdf", sep = '')))
  plot(coda_samples[[i]])
  dev.off()
}

