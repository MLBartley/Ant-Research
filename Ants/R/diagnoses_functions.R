
#' Diagnoses of Penalized HMM Model 
#'
#' @param mcmc_matrix 
#' @param Time
#'
#' @return
#' @export
#'

penalty_diagnosis <- function(mcmc_matrix, Time){
  
  #separate back out the results matrix
  
  st_rates_param  <- mcmc_matrix[1:2, ]
  switch_rate_param <- mcmc_matrix[4:5, ]
  st_ptm_param <- mcmc_matrix[6:9, ]
  states_param <- mcmc_matrix[-(1:9), ]
  
  source("http://www.stat.psu.edu/~mharan/batchmeans.R")
  
  st_rate_high <- st_rates_param[1, ] + st_rates_param[2, ]
  
  lambda.est <- apply(rbind(st_rates_param[1, ], st_rate_high), 1, bm)
  lambda.var <- apply(st_rates_param, 1, quantile, probs = c(0.025, 0.975), 
    na.rm = TRUE)
  
  gamma.est <- apply(switch_rate_param, 1, bm)
  gamma.var <- apply(switch_rate_param, 1, quantile, probs = c(0.025, 0.975), 
    na.rm = TRUE)
  
  
  P.est <- apply(st_ptm_param, 1, bm)
  P.var <- apply(st_ptm_param, 1, quantile, probs = c(0.025, 0.975, na.rm = T))
  
  
  X.est <- matrix(NA, Time, 1, T)
  
  for (t in 1:Time) {
    X.est[t, 1] <- mean(states_param[t, ])
  }
  
  # sum.it <- 0
  # 
  # for (i in 1:n_mcmc) {
  #   for (t in 1:Time) {
  #     sum.it <- sum.it + (osa_param[t, i] - data[t])^2
  #   }
  # }
  # 
  # MSPE.1SA <- 1/n_mcmc * 1/Time * sum.it
  
  # visualization
  if (fig_save == TRUE) {
    jpeg( file = paste(fig_path, fig_name, round(penalty, 11), ".diagnostics", ".jpg", sep = ""))
  }
  
  # plot the estimation runs.
  
  col <- c("#120d08", "#bc5356", "#538bbc", "#53bc84")
  
  
  
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0) + 1, mar = c(1, 1, 1, 1) + 
      3)
  
  
  # Rate Parameters
  plot(0, 0, xlab = "MCMC Runs", ylab = "Rates (per second)", ylim = c(0, 
    max(switch_rate_param, st_rates_param[1:2, ])/delta_t), xlim = c(0, n_mcmc), 
    type = "n", cex.lab = 1)
  lines(1:n_mcmc, (st_rates_param[1, ] + st_rates_param[2, ])/delta_t, 
    col = col[1])
  lines(1:n_mcmc, st_rates_param[1, ], col = col[2])
  
  lines(1:n_mcmc, switch_rate_param[1, ], col = col[3])
  lines(1:n_mcmc, switch_rate_param[2, ], col = col[4])
  
  # X params
  
  # P
  plot(0, 0, xlab = "MCMC Runs", ylab = "Probability Matrix for State Switching", 
    ylim = c(0, max(st_ptm_param)), xlim = c(0, n_mcmc), type = "n", cex.lab = 1)
  for (i in 1:(4)) {
    lines(1:n_mcmc, st_ptm_param[i, ], col = col[i])
  }
  
  # Single X
  X <- states_param[sample(1:Time, 1), ]
  
  plot(0, 0, xlab = "MCMC Runs", ylab = "Single X", ylim = c(0, max(X)), 
    xlim = c(0, n_mcmc), type = "n", cex.lab = 1)
  lines(1:n_mcmc, X, col = col[4])
  
  # States over time plot(X.est, type = 'l')
  plot(round(X.est), type = "l")
  
  
  title(main = "Diagnostic Plots", outer = T)
  
  if (fig_save == TRUE) {
    dev.off()
  }

}