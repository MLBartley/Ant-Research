
#' Diagnoses of Penalized HMM Model 
#'
#' @param mcmc_matrix 
#' @param Time
#' @param fig_path
#' @param fig_name
#' @param penalty
#'
#' 
#'
#' @return
#' @export
#'


penalty_diagnosis <- function(mcmc_matrix, Time, fig_path, fig_name, penalty){
  
  source("http://www.stat.psu.edu/~mharan/batchmeans.R")
  
  #batch means (mean, s.error)
  bmeans <- bmmat(t(mcmc_matrix))
  
pdf( file = paste(fig_path, fig_name, round(penalty, 11), ".diagnostics", ".pdf", sep = ""))  
  #estimates vs sample size
  
    #lamdas - trop.rate(L/H)
    
    par(mfrow = c(1, 2))  
    
    estvssamp(samp = mcmc_matrix[1, ], plotname = "Trophallaxis Rate: Low")
    estvssamp(samp = (mcmc_matrix[1, ] + mcmc_matrix[2, ]), plotname = "Trophallaxis Rate: High")
    
    #gammas - st_switch_rate(LH/HL)
    
    estvssamp(samp = mcmc_matrix[4, ], plotname = "State Switching Rate: Low to High")
    estvssamp(samp = mcmc_matrix[5, ], plotname = "State Switching Rate: High to Low")
    
    #Probabilities
    
    par(mfrow = c(2, 2))
    
    estvssamp(samp = mcmc_matrix[6, ], plotname = "Probability Transitions: Low to Low")
    estvssamp(samp = mcmc_matrix[7, ], plotname = "State Switching Rate: Low to High")
    estvssamp(samp = mcmc_matrix[8, ], plotname = "State Switching Rate: High to Low")
    estvssamp(samp = mcmc_matrix[9, ], plotname = "State Switching Rate: High to High")
    
    
    #Random time for state chain - 

    par(mfrow = c(1, 1))
    
    sample_x <- sample(x = 10:nrow(mcmc_matrix), 1)
    estvssamp(samp = mcmc_matrix[sample_x, ], plotname = "Behavior State: Random Time Chain")
    
  
  # visualization
    
  #lambdas - troph.rate  
  plot(coda::mcmc(t(mcmc_matrix[1:2, ])))
  
  #gammas - st_switch_rate
  plot(coda::mcmc(t(mcmc_matrix[4:5, ])))
  
  #ptm 
  plot(coda::mcmc(t(mcmc_matrix[6:9, ])))
  
  #random time for state chain
  # plot(coda::mcmc(t(mcmc_matrix[sample_x, ])))
    

  #   # jpeg( file = paste(fig_path, fig_name, round(penalty, 11), ".diagnostics", ".jpg", sep = ""))
  # 
  # 
  # # plot the estimation runs.
  # 
  # col <- c("#120d08", "#bc5356", "#538bbc", "#53bc84")
  # 
  # 
  # 
  # par(mfrow = c(2, 2), oma = c(0, 0, 2, 0) + 1, mar = c(1, 1, 1, 1) + 
  #     3)
  # 
  # 
  # # Rate Parameters
  # plot(0, 0, xlab = "MCMC Runs", ylab = "Rates (per second)", ylim = c(0, 
  #   max(switch_rate_param, st_rates_param[1:2, ])/delta_t), xlim = c(0, n_mcmc), 
  #   type = "n", cex.lab = 1)
  # lines(1:n_mcmc, (st_rates_param[1, ] + st_rates_param[2, ])/delta_t, 
  #   col = col[1])
  # lines(1:n_mcmc, st_rates_param[1, ], col = col[2])
  # 
  # lines(1:n_mcmc, switch_rate_param[1, ], col = col[3])
  # lines(1:n_mcmc, switch_rate_param[2, ], col = col[4])
  # 
  # # X params
  # 
  # # P
  # plot(0, 0, xlab = "MCMC Runs", ylab = "Probability Matrix for State Switching", 
  #   ylim = c(0, max(st_ptm_param)), xlim = c(0, n_mcmc), type = "n", cex.lab = 1)
  # for (i in 1:(4)) {
  #   lines(1:n_mcmc, st_ptm_param[i, ], col = col[i])
  # }
  # 
  # # Single X
  # X <- states_param[sample(1:Time, 1), ]
  # 
  # plot(0, 0, xlab = "MCMC Runs", ylab = "Single X", ylim = c(0, max(X)), 
  #   xlim = c(0, n_mcmc), type = "n", cex.lab = 1)
  # lines(1:n_mcmc, X, col = col[4])
  # 
  # # States over time plot(X.est, type = 'l')
  # plot(round(X.est), type = "l")
  # 
  # 
  # title(main = "Diagnostic Plots", outer = T)
  
    dev.off()
  
  return(bmeans)

}