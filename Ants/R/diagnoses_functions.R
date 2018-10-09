
#' Diagnoses of Penalized HMM Model
#'
#' @param mcmc_matrix
#' @param states
#' @param ant_file
#' @param chamber
#' @param Time
#' @param fig_path
#' @param fig_name
#' @param penalty
#' @param covariate
#'
#'
#'
#' @return Diagnoses plots of parameters in model + state estimates plot. Works with all models.
#' @export
#'


penalty_diagnosis <- function(mcmc_matrix, states, ant_file, chamber,
                              Time, fig_path, fig_name, penalty, covariate = NULL){


n <- states


  source("http://www.stat.psu.edu/~mharan/batchmeans.R")

  #batch means (mean, s.error)
  bmeans <- bmmat(t(mcmc_matrix))

pdf( file = paste(fig_path, fig_name, round(penalty, 11), ".diagnostics", ".pdf", sep = ""))

  #estimates vs sample size

    #lamdas - trop.rate(L/H)

    par(mfrow = c(1, 2))

    estvssamp(samp = mcmc_matrix[1, ], plotname = "Trophallaxis Rate: Low")

    if (n == 2) {
      estvssamp(samp = (mcmc_matrix[1, ] + mcmc_matrix[2, ]), plotname = "Trophallaxis Rate: High")
    }else{
      estvssamp(samp = (mcmc_matrix[1, ] + mcmc_matrix[2, ]),
        plotname = "Trophallaxis Rate: Medium")
      estvssamp(samp = (mcmc_matrix[1, ] + mcmc_matrix[2, ] + mcmc_matrix[3, ]),
        plotname = "Trophallaxis Rate: High")
    }

  if (is.null(covariate) == TRUE ){
    #gammas - st_switch_rate(LH/HL)
    if (n == 2) {
      estvssamp(samp = mcmc_matrix[4, ], plotname = "State Switching Rate: Low to High")
    estvssamp(samp = mcmc_matrix[5, ], plotname = "State Switching Rate: High to Low")
    }else{
      estvssamp(samp = mcmc_matrix[5, ], plotname = "State Switching Rate: Low to High")
      estvssamp(samp = mcmc_matrix[6, ], plotname = "State Switching Rate: High to Low")
      estvssamp(samp = mcmc_matrix[7, ], plotname = "State Switching Rate: Medium to Low")
      estvssamp(samp = mcmc_matrix[8, ], plotname = "State Switching Rate: Medium to High")
      estvssamp(samp = mcmc_matrix[9, ], plotname = "State Switching Rate: Low to Medium")
      estvssamp(samp = mcmc_matrix[10, ], plotname = "State Switching Rate: High to Medium")
    }

    #Probabilities

    par(mfrow = c(2, 2))
if (n == 2) {
    estvssamp(samp = mcmc_matrix[6, ], plotname = "Probability Transitions: Low to Low")
    estvssamp(samp = mcmc_matrix[7, ], plotname = "Probability Transitions: Low to High")
    estvssamp(samp = mcmc_matrix[8, ], plotname = "Probability Transistions: High to Low")
    estvssamp(samp = mcmc_matrix[9, ], plotname = "Probability Transistions: High to High")
}else{
  estvssamp(samp = mcmc_matrix[11, ], plotname = "Probability Transitions: Low to Low")
  estvssamp(samp = mcmc_matrix[12, ], plotname = "Probability Transitions: Low to Medium")
  estvssamp(samp = mcmc_matrix[13, ], plotname = "Probability Transitions: Low to High")

  estvssamp(samp = mcmc_matrix[14, ], plotname = "Probability Transitions: Medium to Low")
  estvssamp(samp = mcmc_matrix[15, ], plotname = "Probability Transitions: Medium to Medium")
  estvssamp(samp = mcmc_matrix[16, ], plotname = "Probability Transitions: Medium to High")

  estvssamp(samp = mcmc_matrix[17, ], plotname = "Probability Transitions: High to Low")
  estvssamp(samp = mcmc_matrix[18, ], plotname = "Probability Transitions: High to Medium")
  estvssamp(samp = mcmc_matrix[19, ], plotname = "Probability Transitions: High to High")

}
  }else{
    estvssamp(samp = mcmc_matrix[4, ], plotname = "State Switching Rates: e^beta_0LH")
    estvssamp(samp = mcmc_matrix[5, ], plotname = "State Switching Rates: e^beta_0HL")
    estvssamp(samp = mcmc_matrix[6, ], plotname = "State Switching Rates: beta_1LH")
    estvssamp(samp = mcmc_matrix[7, ], plotname = "State Switching Rates: beta_1HL")
    estvssamp(samp = mcmc_matrix[8, ], plotname = "State Swithc Rates: cov_exponent")
  }




    #Random time for state chain -

    par(mfrow = c(1, 1))

    sample_x <- sample(x = 10:nrow(mcmc_matrix), 1)
    if (n != 2) { sample_x <- sample(20:nrow(mcmc_matrix), 1)}
    estvssamp(samp = mcmc_matrix[sample_x, ], plotname = "Behavior State: Random Time Chain")


    # estvssamp(samp = mcmc_matrix[sample_x, ], plotname = "Behavior State: Random Time Chain")
    #

  # visualization

  #lambdas - troph.rate
    if (n == 2) {
      plot(coda::mcmc(t(mcmc_matrix[1:2, ])))

  if (is.null(covariate) == TRUE){
    #gammas - st_switch_rate
  plot(coda::mcmc(t(mcmc_matrix[4:5, ])))

  #ptm
  plot(coda::mcmc(t(mcmc_matrix[6:9, ])))
  }else{
    plot(coda::mcmc(t(mcmc_matrix[4:8, ]))) #e^betas and betas
  }


  #rounded estimates over time for sto process

   plot(round(bmeans[-(1:9), 1]), type = "l")
    }else{
      plot(coda::mcmc(t(mcmc_matrix[1:3, ])))

      #gammas - st_switch_rate
      plot(coda::mcmc(t(mcmc_matrix[5:10, ])))

      #ptm
      plot(coda::mcmc(t(mcmc_matrix[11:19, ])))

      #rounded estimates over time for sto process

      plot(round(bmeans[-(1:19), 1]), type = "l")
    }


   #fancy plot


   # plot(round(bmeans[-(1:9), 1 ]), type = "l")


   #fancy plot

   # needed for final graphic
   location <- ant_file$Location
   start <- ant_file$start_time
   # chamber <- chamber

   if (length(unique(location)) != 1){

     if (chamber == "queen") {
       start = start[which(location == 1)]
     }else{
       start = start[which(location == 4)]
     }
   }

   start <- sort(start)
   int.num <- length(start)
   maxtime <- Time


   par(mfrow = c(1, 1))


   ## High Density - 4 Hours
   plot(start, 1:int.num, xlab = "Seconds", ylab = "Cumulative Interaction Count",
     xlim = c(0, Time))
   states <- bmeans[-(1:9), 1]  #from code above
   if (n != 2) { states <- bmeans[-(1:19), 1]}
   rr <- rle(states)
   rr$values <- round(rr$values, digits = 0)
   embedded.chain <- rr$values
   cs <- c(0, cumsum(rr$lengths)) * delta_t - delta_t
<<<<<<< HEAD
   cols <- c("#bc535644", "#538bbc44")
   if (n!= 2) {cols <- c("#bc535644", "#bcb85344", "#5356bc44")}
=======
   cols <- c("#bc535644", "#538bbc44", "#56bc5344")
>>>>>>> 8cc3341e4c9336c9a157eb5736a8c32b9c22651a
   for (j in 1:length(embedded.chain)) {
     rect(cs[j], 0, cs[j + 1], int.num, col = cols[embedded.chain[j]],
       density = NA)

   }
   points(start, 1:int.num,
     xlim = c(0, maxtime))

  if (is.null(covariate) == F) {
<<<<<<< HEAD
      points(which(covariate == 0), rep(0, length(which(covariate == 0))), pch = 8, col = "black")
=======
    points(which(covariate == 0), rep(0, length(which(covariate == 0))), pch = 8, col = "#53bc84")
>>>>>>> 8cc3341e4c9336c9a157eb5736a8c32b9c22651a
  }
    dev.off()

  return(bmeans)

}
