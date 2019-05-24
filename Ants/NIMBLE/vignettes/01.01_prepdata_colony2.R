
library(dplyr)
library(magrittr)

source(here::here("R", "data_prep_functions.R"))




#data directory
dat.dir <- "../data-raw/"


#column names - helps when only pulling in those columns, no extra

col_names <- c("Location", "Ant_ID", "Ant_ID_partner", "start_time", "end_time")


col2_high4 <- read.csv(here::here("data-raw", "Colony2_trophallaxis_high_density_4hr.csv"))
col2_low4 <- read.csv(here::here("data-raw", "Colony2_trophallaxis_low_density_4hr.csv"))

#removes any extra columns, rows, and adds column names - depends on col_names being correct length
col2_high4 <- col2_high4[, 1:length(col_names)]
col2_high4 <- col2_high4 %>%
  tidyr::drop_na()
colnames(col2_high4) <- col_names

col2_low4 <- col2_low4[, 1:length(col_names)]
col2_low4 <- col2_low4 %>%
  tidyr::drop_na()
colnames(col2_low4) <- col_names

#check for correct class for data (numberic, etc)


col2_high4_5 <- prep_troph_data(col2_high4, hours = 4, delta_t = 5)
col2_low4_5 <- prep_troph_data(col2_low4, hours = 4, delta_t =  5)


load(here::here("data-raw", "Colony2_covariate_low_density_4hr.Rda"))
cov_col2_low4 <- prep_inout_data(covariate, delta_t = 1, hours = 4)

out <- list(col2_high4_5 = col2_high4_5,
            col2_low4_5 = col2_low4_5,
            cov_col2_low4 = cov_col2_low4)
save(out, file = "./NIMBLE/data-prepped/col2preppeddata.Rdata")


