
DT_mcmc_troph <- function(starts_data, ant_file, chamber, title, a, b, c, d, theta,
                          states = 2, n_mcmc, delta_t, hours, param_start,
                          data_out, modelsrun_out,
                          fig_save = TRUE, fig_path, fig_name, plot_title) {

  data <- starts_data
  Time <- length(data)
  n <- states
  delta <- rep(1/n, n)

  states_start <- param_start$X
  st_rates_start <- param_start$lambda
  state_ptm_start <- param_start$P

  # needed for final graphic
  location <- ant_file$Location
  start <- ant_file$start_time
  # chamber <- chamber

  if (length(unique(location)) != 1){

    if(chamber == "queen") {
      start = start[which(location == 1)]
    }else {
      start = start[which(location == 4)]
    }
  }

  start <- sort(start)
  int.num <- length(start)
  maxtime <- hours * 60 * 60


  # homes Build Homes for X(1:T), lambda(1:n), and P(nXn) and gam
  # vectors

  states_param <- matrix(data = rep(NA, Time * 1001), nrow = Time, ncol = 1001,
                         byrow = T)
  row.names(states_param) <- rep("State at Random Time", nrow(states_param))

  st_rates_param <- st_rates_param <- matrix(data = NA, nrow = n + 1, ncol = 1001,
                                             byrow = T)
  row.names(st_rates_param) <- c("trop.rate.low", "trop.rate.change",
                                 "trop.rate.high")

  st_ptm_param <- matrix(data = NA, nrow = n * n, ncol = 1001,
                         byrow = T)
  row.names(st_ptm_param) <- c("LL", "LH", "HL", "HH")

  gam <- matrix(NA, nrow = Time, ncol = n, byrow = T)

  starts_low <- matrix(NA, nrow = Time, ncol = 1001)
  starts_high <- matrix(NA, nrow = Time, ncol = 1001)

  osa_param <- matrix(NA, Time, 1001, T)
  osa_final <- 0

  ## Initialize parameters - MOVE START VALUES TO OUTSIDE FUNCTION?

  states_param[, 1] <- states_start

  osa_param[, 1] <- data

  st_rates_param[1, 1] <- st_rates_start[1]  #lambda low
  st_rates_param[2, 1] <- st_rates_start[2] - st_rates_start[1]  #change in lambda

  st_rates_param[3, 1] <- st_rates_param[1, 1] + st_rates_param[2, 1]  #lambda high, not needed, just a reminder

  ptm_matrix <- state_ptm_start

  st_ptm_param[, 1] <- as.vector(t(ptm_matrix))
  # holds all st_ptm_parameter values over runs


  # HOW BEST TO CODE SPLITTING OF DATA?
  for (t in 1:Time) {
    if (states_param[t, 1] == 1) {
      starts_low[t, 1] <- data[t]
      starts_high[t, 1] <- 0
    } else {
      split <- rmultinom(1, size = data[t], prob = c(st_rates_param[1,
                                                                    1], st_rates_param[2, 1]))
      starts_low[t, 1] <- split[1]
      starts_high[t, 1] <- split[2]
    }

  }

  l.1000 <- 2 #one is prefilled with 'previous' run

  #start each save file
  results <- rbind(st_rates_param,
                   st_ptm_param,
                   states_param)

  osa_results <- osa_param

  write.csv(t(results[, 1]), file = data_out)

  ## Gibbs Updates



  for (l in 2:n_mcmc) {

    # print out every 10 iterations completed
    if (l %% 1000 == 0)
      cat(paste("iteration", l, "complete\n"))


    m <- matrix(data = rep(0, n * n), nrow = n, ncol = n)
    rownames(m) <- c("Low", "High")
    colnames(m) <- c("Low", "High")
    # number states going from i to j, refreshes every run

    ptm_matrix <- matrix(data = c(st_ptm_param[, l.1000 - 1]), nrow = n, ncol = n,
                         byrow = T)

    st_rates_low <- st_rates_param[1, l.1000 - 1]
    st_rates_high <- st_rates_low + st_rates_param[2, l.1000 - 1]



    # X Parameters, split into X_1, X_{2:Time-1}, X_Time

    gam[1, 1] <- st_rates_low^data[1] * exp(-st_rates_low) * delta[1] *
      ptm_matrix[1, states_param[2, l.1000 - 1]]

    gam[1, 2] <- st_rates_high^data[1] * exp(-st_rates_high) * delta[1] *
      ptm_matrix[1, states_param[2, l.1000 - 1]]


    states_param[1, l.1000] <- sample(x = (1:n), size = 1, prob = gam[1,
                                                                      ])

    m[states_param[1, l.1000], states_param[1, l.1000]] <- m[states_param[1, l.1000], states_param[1,
                                                                                                   l.1000]] + 1

    osa_param[1, l.1000] <- st_rates_low * ptm_matrix[states_param[1, l.1000], 1] +
      st_rates_high * ptm_matrix[states_param[1, l.1000], 2]

    for (t in 2:(Time - 1)) {


      gam[t, 1] <- st_rates_low^data[t] * exp(-st_rates_low) * ptm_matrix[states_param[t -
                                                                                         1, l.1000 - 1], 1] * ptm_matrix[1, states_param[t + 1, l.1000 - 1]]

      gam[t, 2] <- st_rates_high^data[t] * exp(-st_rates_high) *
        ptm_matrix[states_param[t - 1, l.1000 - 1], 2] * ptm_matrix[2, states_param[t +
                                                                                      1, l.1000 - 1]]


      states_param[t, l.1000] <- sample(x = (1:n), 1, prob = gam[t, ])

      m[states_param[t - 1, l.1000], states_param[t, l.1000]] <- m[states_param[t - 1,
                                                                                l.1000], states_param[t, l.1000]] + 1

      osa_param[t, l.1000] <- st_rates_low * ptm_matrix[states_param[t, l.1000],
                                                        1] + st_rates_high * ptm_matrix[states_param[t, l.1000], 2]

    }

    gam[Time, 1] <- st_rates_low^data[Time] * exp(-st_rates_low) *
      ptm_matrix[states_param[Time - 1, l.1000 - 1], 1]

    gam[Time, 2] <- st_rates_high^data[Time] * exp(-st_rates_high) *
      ptm_matrix[states_param[Time - 1, l.1000 - 1], 2]


    states_param[Time, l.1000] <- sample(x = 1:n, 1, prob = gam[Time, ])

    m[states_param[Time - 1, l.1000], states_param[Time, l.1000]] <- m[states_param[Time -
                                                                                      1, l.1000], states_param[Time, l.1000]] + 1

    osa_param[Time, l.1000] <- st_rates_low * ptm_matrix[states_param[Time, l.1000],
                                                         1] + st_rates_high * ptm_matrix[states_param[Time, l.1000], 2]


    # Split data (N_t) into N_Ht, N_Lt
    for (t in 1:Time) {
      if (states_param[t, l.1000] == 1) {
        starts_low[t, l.1000] <- data[t]
        starts_high[t, l.1000] <- 0
      } else {
        split <- rmultinom(1, size = data[t], prob = c(st_rates_param[1,
                                                                      l.1000 - 1], st_rates_param[2, l.1000 - 1]))
        starts_low[t, l.1000] <- split[1]
        starts_high[t, l.1000] <- split[2]
      }

    }

    # Lambda and P parameters
    for (h in 1:n) {

      ptm_matrix[h, ] <- (rdirichlet(n = 1, alpha = theta[h, ] +
                                       m[h, ]))

    }

    st_rates_param[1, l.1000] <- rgamma(n = 1, shape = sum(starts_low[, l.1000]) + a,
                                        rate = Time + b)

    st_rates_param[2, l.1000] <- rgamma(n = 1, shape =
                                          sum(starts_high[which(states_param[, l.1000] == 2)]) +
                                          c, rate = sum(m[2, ]) + d)

    st_rates_param[3,l.1000] <- st_rates_param[1, l.1000] + st_rates_param[2, l.1000]


    st_ptm_param[, l.1000] <- as.vector(t(ptm_matrix))

    # move iteration forward
    l.1000 <- l.1000 + 1
    ###################################### Every 1000 iterations

    if (l %% 1000 == 0) {
      # print out every 1000 iterations completed
      cat(paste("iteration", l, "complete\n"))

      #save off chunk of 1000 chains

      results <- rbind(st_rates_param,
                       st_ptm_param,
                       states_param)

      data_mat <- matrix(rep(data, ncol(osa_param[, c(-1, -1001)])),
                         byrow = F, nrow = nrow(osa_param[, c(-1, -1001)]))
      sum.it <- sum((osa_param[, c(-1, -1001)] - data_mat)^2)
      osa_final <- osa_final + sum.it

      # pdf( file = paste(fig_path, fig_name, log(penalty), ".tracesnapshot", ".pdf", sep = ""))
      # for (i in 2:9) {
      #   plot((results[i, ]), type = "l") #snapshot trace plots for
      #
      # }
      # dev.off()

      write.table(t(results[, 2:1000 ]), file = data_out,
                  append = T, col.names = F, sep = ',')

      # write.table(t(osa_results[, 2:1000]), file = osa_out,
      #   append = T, col.names = F, sep = ',')


      #reset chain home
      #save most recent

      st_rates_param[, 1] <- st_rates_param[, 1000]
      st_ptm_param[ , 1] <- st_ptm_param[ , 1000]
      states_param[ , 1] <- states_param[ , 1000]
      osa_param[, 1] <- osa_param[ , 1000]

      starts_low[ , 1] <- starts_low[ , 1000]
      starts_high[ , 1] <- starts_high[ , 1000]


      #reset chain home index
      l.1000 <- 2

    }
  }


  ## Rescale Lambda parameters into per minute segments (instead of
  ## delta_t time segments)
  #
  #   st_rates_scale <- st_rates_param/delta_t
  #
  #
  #
  #   ## Compile the Estimates
  #
  #   ## X1:XT, Lambda, Pmatrix
  #
  #   # homes
  #   states_est <- matrix(data = rep(NA, Time), nrow = Time, ncol = 1)
  #   st_rates_est <- matrix(data = NA, nrow = n + 1, ncol = 1)
  #   st_ptm_est <- matrix(data = rep(NA, n * n), nrow = n * n, ncol = 1)
  #
  #   # estimation
  #   source("http://www.stat.psu.edu/~mharan/batchmeans.R")
  #
  #
  #   for (t in 1:Time) {
  #     states_est[t, 1] <- mean(states_param[t, ])
  #   }
  #
  #
  #   st_rates_high <- st_rates_scale[2, ] + st_rates_scale[1, ]
  #
  #   st_rates_est <- apply(rbind(st_rates_scale[1, ], st_rates_high), 1, bm)
  #   st_rates_var <- apply(st_rates_scale, 1, quantile, probs = c(0.025,
  #     0.975), na.rm = TRUE)
  #
  #   st_ptm_est <- apply(st_ptm_param, 1, bm)
  #   st_ptm_var <- apply(st_ptm_param, 1, quantile, probs = c(0.025, 0.975), na.rm = T)
  #
  #   sum.it <- sum((osa_param - data)^2)
  #
  #   # for (i in 1:n_mcmc) {
  #   #   for (t in 1:Time) {
  #   #     sum.it <- sum.it +  (osa_param[t, i] - data[t])^2
  #   #   }
  #   # }
  #
  #   MSPE.1SA <- 1/n_mcmc * 1/Time * sum.it
  #
  #   # plot the estimation runs.
  #   col <- c("#120d08", "#bc5356", "#538bbc", "#53bc84")
  #
  #   if (fig_save == TRUE) {
  #     jpeg(file = paste(fig_path, fig_name, ".diagnostics" , theta[1], ".jpg", sep = ""))
  #
  #   }
  #
  #   # lambda
  #   par(mfrow = c(2, 2), oma = c(0, 0, 2, 0) + 1, mar = c(1, 1, 1,
  #     1) + 3)
  #
  #
  #   plot(0, 0, xlab = "MCMC Runs", ylab = "Lambda (scaled per second)",
  #     ylim = c(0, max(st_rates_high)), xlim = c(0, n_mcmc), type = "n",
  #     cex.lab = 1)
  #   lines(1:n_mcmc, st_rates_scale[1, ], col = col[1])
  #   lines(1:n_mcmc, st_rates_high, col = col[2])
  #
  #
  #   # for(i in 1:n){ lines(1:n_mcmc, st_rates_scale[i, ], col = col[i])
  #   # }
  #
  #
  #   # P
  #   plot(0, 0, xlab = "MCMC Runs", ylab = "P", ylim = c(0, max(st_ptm_param)),
  #     xlim = c(0, n_mcmc), type = "n", cex.lab = 1)
  #   for (i in 1:(n * n)) {
  #     lines(1:n_mcmc, st_ptm_param[i, ], col = col[i])
  #   }
  #
  #   # Single X
  #   X <- states_param[sample(1:Time, 1), ]
  #   plot(0, 0, xlab = "MCMC Runs", ylab = "Single X", ylim = c(0,
  #     max(X)), xlim = c(0, n_mcmc), type = "n", cex.lab = 1)
  #   lines(1:n_mcmc, X, col = col[2])
  #
  #   # States over time
  #   plot(round(states_est), type = "l", lwd = 3, cex.lab = 1, col = col[1])
  #
  #   title(main = title, outer = T)
  #
  #   if (fig_save == TRUE) {
  #     dev.off()
  #   }
  #   ######################################################### Fancy Plots with Background Colors
  #   if (fig_save == TRUE) {
  #     jpeg(file = paste(fig_path, fig_name, ".states" , theta[1], ".jpg", sep = ""))
  #
  #   }
  #
  #   par(mfrow = c(1, 1))
  #
  #
  #   if (length(unique(location)) == 1) {
  #
  #     ## High Density - 4 Hours
  #     plot(start, 1:int.num, main = plot_title, xlab = "Seconds", ylab = "Cumulative
  #       Interaction Count",
  #       xlim = c(0, maxtime))
  #     states <- states_est  #from code above
  #     rr <- rle(states[, 1])
  #     rr$values <- round(rr$values, digits = 0)
  #     embedded.chain <- rr$values
  #     cs <- c(0, cumsum(rr$lengths)) * delta_t - delta_t
  #     cols <- c("#bc535644", "#538bbc44")
  #
  #     for (j in 1:length(embedded.chain)) {
  #       rect(cs[j], 0, cs[j + 1], int.num, col = cols[embedded.chain[j]],
  #         density = NA)
  #     }
  #
  #     points(start, 1:int.num, main = plot_title, xlab = "Seconds", ylab = "Cumulative
  #       Interaction Count",
  #       xlim = c(0, maxtime))
  #   } else {
  #     # Low Density - 4 Hours
  #
  #     # if (chamber == "queen") {
  #     #   start <- start_1
  #     #   int.num <- length(start)
  #     #
  #     # } else {
  #     #   start <- start_4
  #     #   int.num <- length(start)
  #     #
  #     # }
  #
  #     plot(start, 1:int.num, main = plot_title, xlab = "Seconds", ylab = "Cumulative
  #       Interaction Count",
  #       xlim = c(0, maxtime))
  #     states <- states_est
  #     rr <- rle(states[, 1])
  #     rr$values <- round(rr$values, digits = 0)
  #     embedded.chain <- rr$values
  #     cs <- c(0, cumsum(rr$lengths)) * delta_t - delta_t
  #     cols <- c("#bc535644", "#538bbc44")
  #     for (j in 1:length(embedded.chain)) {
  #       rect(cs[j], 0, cs[j + 1], int.num, col = cols[embedded.chain[j]],
  #         density = NA)
  #     }
  #     points(start, 1:int.num, main = plot_title, xlab = "Seconds", ylab = "Cumulative
  #       Interaction Count",
  #       xlim = c(0, maxtime))
  #
  #   }
  #
  #
  #   if(fig_save == TRUE){
  #     dev.off()
  #   }
  #
  #   list(states_est = states_est, st_rates_est = st_rates_est,
  #     st_ptm_est = st_ptm_est, P.run = st_ptm_param, MSPE = MSPE.1SA)

  MSPE = 1/n_mcmc * 1/Time * osa_final

  out <- c(theta[1, 1], MSPE, n_mcmc)


  write.table(t(out), file = modelsrun_out,
              append = T, col.names = F, sep = ',')


}
