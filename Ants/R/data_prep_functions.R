#' Data Preparation for Ant Trophallaxis Data File
#'
#' The purpose of this function is to take in .csv files and
#' apply necessary changes to the data format. Changes include
#' removing rows/columns of NA, removing redundant entries,
#' sorting by entrance times, and binning data into smaller
#' chunks.
#'
#' @param data Ant trophallaxis data file in .csv format, already loaded into R.
#' @param hours Total time ants are observed.
#' @param delta_t Time increments for (starting time) data to be binned into.
#'
#' @return This function will return one or three (depending on high or low density,
#'   respectively) file lists with amount of interaction starting within
#'   binned time chunk (delta_t).
#' @export
#' @examples
#'  prep_troph_data(high4, 60)
#'
#'

prep_troph_data <- function(data, hours, delta_t) {

  seconds <- hours * 60 * 60

  # remove columns with NA values
  data <- data[, colSums(is.na(data)) < nrow(data)]

  # remove rows with NA values
  data <- data[rowSums(is.na(data)) == 0,]

  # remove duplicate entries
  data <- data[seq(1, nrow(data), by = 2), ]

  # order data frame by start time so plot works better
  data %<>% arrange(start_time)

  # combining the (start time) data into delta_t increments
  if (length(unique(data$Location)) == 1) {

    # High Density Data
      start <- data$start_time
      start_persec <- rep(0, seconds)
      bin_starts <- rep(0, floor(seconds/delta_t))
      mint <- 0

      for (t in 1:length(start_persec)) {
        start_persec[t] <- length(which(start > mint & start <= mint + 1))
        #how many interactions start in second increment
        mint <- mint + 1
      }

      for (t in 1:length(bin_starts)) {
        bin_starts[t] <- length(which(start > mint & start <= mint + delta_t))
        #how many interactions start in delta_t increment
        mint <- mint + delta_t
      }



     list <- list(data = data, starts_persec = start_persec, bin_starts = bin_starts)
  } else {
    ## Separate Low Density by Location
    low_1 <- data[which(data$Location == 1), ]

    low_4 <- data[which(data$Location == 4), ]
    low_4 %<>% arrange(start_time)

    # Low Density Data: Both Locations
      start <- data$start_time
      start_persec <- rep(0, seconds)
      bin_starts <- rep(0, floor(seconds/delta_t))  #max time same ~7200
      mint <- 0

      for (t in 1:length(start_persec)) {
        start_persec[t] <- length(which(start > mint & start <= mint + 1))
        #how many interactions start in second increment
        mint <- mint + 1
      }
      for (t in 1:length(bin_starts)) {
        bin_starts[t] <- length(which(start > mint & start <= mint + delta_t))
        mint <- mint + delta_t
      }


    # Low Density Data: Location 1
      start <- low_1$start_time
      start_persec_l1 <- rep(0, seconds)
      bin_starts_l1 <- rep(0, floor(seconds/delta_t))
      mint <- 0

      for (t in 1:length(start_persec)) {
        start_persec_l1[t] <- length(which(start > mint & start <= mint + 1))
        #how many interactions start in second increment
        mint <- mint + 1
      }
      for (t in 1:length(bin_starts)) {
        bin_starts_l1[t] <- length(which(start > mint & start <= mint + delta_t))
        mint <- mint + delta_t
      }


    # Low Density Data: Location 4
      start <- low_4$start_time
      start_persec_l4 <- rep(0, seconds)
      bin_starts_l4 <- rep(0, floor(seconds/delta_t))
      mint <- 0

      for (t in 1:length(start_persec)) {
        start_persec_l4[t] <- length(which(start > mint & start <= mint + 1))
        #how many interactions start in second increment
        mint <- mint + 1
      }
      for (t in 1:length(bin_starts)) {
        bin_starts_l4[t] <- length(which(start > mint & start <= mint + delta_t))
        mint <- mint + delta_t
      }



   list =  list(data = data, starts_persec = start_persec, queen_starts_persec = start_persec_l1, entrance_start_persec = start_persec_l4, bin_starts = bin_starts, queen_bin_starts = bin_starts_l1, entrance_bin_starts = bin_starts_l4)
  }
  return(list)
}


#' Data Preparation for Ant In/Out Data File
#'
#' The purpose of this function is to take a .csv file with
#' entrance/exit data and format it for use in mcmc model.
#'
#' @param data
#' @param delta_t
#' @param hours
#' @return (1) List file for binned covariate (forager arrival) data.
#' @export
#' @examples
#'
#' prep.inout.data(inout.high.4, 60, 4)


prep_inout_data = function(data, delta_t, hours) {

  data = data[which(data$Action == "Exit"), ]
  forager_arrivals = sort(data$time)
  forager_arrivals = forager_arrivals[which(duplicated(forager_arrivals) ==
      FALSE)]

  forager_arrivals = c(forager_arrivals, 999999)

  ## Time since Arrivals

  time = hours * 60 * 60
  covariate = rep(NA, time)
  covariate[1] = 3000  #10 minutes since arrival

  arrival = 1

  for (i in 2:time) {

    if (forager_arrivals[arrival] == i) {
      covariate[i] = 0
      arrival = arrival + 1
    } else {
      covariate[i] = covariate[i - 1] + 1
    }

  }

  # bin covariates into similar chunks as the data
  max.time = time

  if (delta_t != 1) {

    c = rep(0, max.time/delta_t)
    mint = 1
    for (t in 1:length(c)) {
      c[t] = min(covariate[mint:(mint + delta_t - 1)])
      mint = mint + delta_t
    }

  } else {
    c = covariate
  }


  list(cov = c)
}
