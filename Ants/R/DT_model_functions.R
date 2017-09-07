 #' Discrete Time MCMC Estimation for Trophallaxis data
#' 
#' The purpose of this function is to find MCMC generated 
#' estimates of (1) - the state (X_t = high/low troph rates) of 
#' the colony at time t (2) - the specific rates of interaction 
#' start times (lambda) of each state (3) - the probability of 
#' moving from one state to another (P matrix)
#' 
#' @param starts_data Data file of number of trophallaxis 
#'   interactions starting per delta_t delta_t.
#' @param ant_file Data file of ant trophallaxis interaction 
#'   information, including Location and start times.
#' @param title Title of diagnostic plots used to check model 
#'   results.
#' @param a Hyperprior for rates of starting interactions during 
#'   "low" state. (lambda_L)
#' @param b Hyperprior for rates of starting interactions during 
#'   "low" state. (lambda_L)
#' @param c Hyperprior for rates of starting interactions during 
#'   "high" state. (lambda_H)
#' @param d Hyperprior for rates of starting interactions during 
#'   "high" state. (lambda_H)
#' @param theta Hyperprior for state probability transition 
#'   matrix  (P)
#' @param states  Number of states, generally 2 ("low", "high").
#' @param n_mcmc  Number of mcmc iterations
#' @param delta_t  Time segment the start data is binned into.
#' @param hours Hours of ant observations reflected in data file.
#' @param param_start Starting values for the chains of estimated
#'   parameters.
#' @param fig_save If "TRUE", plots will save in path provided as
#'   .jpeg files.
#' @param fig_path Path needed to send plot figures to folder.
#' @param fig_name Base name of plot files to be saved.
#'   
#' @return  (1) - estimates of X, lambda, P (2) - 2x2 visual of 
#'   estimates over time (runs)
#' @export
#' @examples
#' theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 
#' out.high = DT.mcmc.troph(data = high.y, title = 'High Density',
#'                          a = 5, b = 2, 
#' theta = theta, states = 2, n_mcmc = 3000)


DT_mcmc_troph <- function(starts_data, ant_file, chamber, title, a, b, c, d, theta, 
                          states = 2, n_mcmc, delta_t, hours, param_start,
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
  
  states_param <- matrix(data = NA, nrow = Time, ncol = n_mcmc, byrow = T)
  
  st_rates_param <- matrix(data = NA, nrow = n + 1, ncol = n_mcmc, byrow = T)
  row.names(st_rates_param) <- c("trop.rate.low", "trop.rate.change", 
    "trop.rate.high")
  
  st_ptm_param <- matrix(data = rep(NA, n * n * n_mcmc), nrow = n * n, 
    ncol = n_mcmc, byrow = T)
  row.names(st_ptm_param) <- c("LL", "LH", "HL", "HH")
  
  gam <- matrix(NA, nrow = Time, ncol = n, byrow = T)
  
  starts_low <- matrix(NA, nrow = Time, ncol = n_mcmc)
  starts_high <- matrix(NA, nrow = Time, ncol = n_mcmc)
  
  osa_param <- matrix(NA, Time, n_mcmc, T)
  
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
  
  
  ## Gibbs Updates
  
  for (l in 2:n_mcmc) {
    
    # print out every 10 iterations completed
    if (l %% 100 == 0) 
      cat(paste("iteration", l, "complete\n"))
    
    
    m <- matrix(data = rep(0, n * n), nrow = n, ncol = n)
    rownames(m) <- c("Low", "High")
    colnames(m) <- c("Low", "High")
    # number states going from i to j, refreshes every run
    
    ptm_matrix <- matrix(data = c(st_ptm_param[, l - 1]), nrow = n, ncol = n, 
      byrow = T)
    
    st_rates_low <- st_rates_param[1, l - 1]
    st_rates_high <- st_rates_low + st_rates_param[2, l - 1]
    
    
    
    # X Parameters, split into X_1, X_{2:Time-1}, X_Time
    
    gam[1, 1] <- st_rates_low^data[1] * exp(-st_rates_low) * delta[1] * 
      ptm_matrix[1, states_param[2, l - 1]]
    
    gam[1, 2] <- st_rates_high^data[1] * exp(-st_rates_high) * delta[1] * 
      ptm_matrix[1, states_param[2, l - 1]]
  
    
    states_param[1, l] <- sample(x = (1:n), size = 1, prob = gam[1, 
      ])
    
    m[states_param[1, l], states_param[1, l]] <- m[states_param[1, l], states_param[1, 
      l]] + 1
    
    osa_param[1, l] <- st_rates_low * ptm_matrix[states_param[1, l], 1] + 
      st_rates_high * ptm_matrix[states_param[1, l], 2]
    
    for (t in 2:(Time - 1)) {
      
      
      gam[t, 1] <- st_rates_low^data[t] * exp(-st_rates_low) * ptm_matrix[states_param[t - 
          1, l - 1], 1] * ptm_matrix[1, states_param[t + 1, l - 1]]
      
      gam[t, 2] <- st_rates_high^data[t] * exp(-st_rates_high) * 
        ptm_matrix[states_param[t - 1, l - 1], 2] * ptm_matrix[2, states_param[t + 
            1, l - 1]]
      
      
      
      states_param[t, l] <- sample(x = (1:n), 1, prob = gam[t, ])
      
      m[states_param[t - 1, l], states_param[t, l]] <- m[states_param[t - 1, 
        l], states_param[t, l]] + 1
      
      osa_param[t, l] <- st_rates_low * ptm_matrix[states_param[t, l], 
        1] + st_rates_high * ptm_matrix[states_param[t, l], 2]
      
    }
    
    gam[Time, 1] <- st_rates_low^data[Time] * exp(-st_rates_low) * 
      ptm_matrix[states_param[Time - 1, l - 1], 1]
    
    gam[Time, 2] <- st_rates_high^data[Time] * exp(-st_rates_high) * 
      ptm_matrix[states_param[Time - 1, l - 1], 2]
    
    
    states_param[Time, l] <- sample(x = 1:n, 1, prob = gam[Time, ])
    
    m[states_param[Time - 1, l], states_param[Time, l]] <- m[states_param[Time - 
        1, l], states_param[Time, l]] + 1
    
    osa_param[Time, l] <- st_rates_low * ptm_matrix[states_param[Time, l], 
      1] + st_rates_high * ptm_matrix[states_param[Time, l], 2]
    
    
    # Split data (N_t) into N_Ht, N_Lt
    for (t in 1:Time) {
      if (states_param[t, l] == 1) {
        starts_low[t, l] <- data[t]
        starts_high[t, l] <- 0
      } else {
        split <- rmultinom(1, size = data[t], prob = c(st_rates_param[1, 
          l - 1], st_rates_param[2, l - 1]))
        starts_low[t, l] <- split[1]
        starts_high[t, l] <- split[2]
      }
      
    }
    
    # Lambda and P parameters
    for (h in 1:n) {
      
    ptm_matrix[h, ] <- (rdirichlet(n = 1, alpha = theta[h, ] + 
          m[h, ]))
      
    }
    
    st_rates_param[1, l] <- rgamma(n = 1, shape = sum(starts_low[, l]) + a,
                                  rate = Time + b)
    
    st_rates_param[2, l] <- rgamma(n = 1, shape = 
                                  sum(starts_high[which(states_param[, l] == 2)]) +
                                  c, rate = sum(m[2, ]) + d)
    
    st_rates_param[3,l] <- st_rates_param[1, l] + st_rates_param[2, l]
    
    
    st_ptm_param[, l] <- as.vector(t(ptm_matrix))
  }
  
  
  ## Rescale Lambda parameters into per minute segments (instead of
  ## delta_t time segments)
  
  st_rates_scale <- st_rates_param/delta_t
  
  
  
  ## Compile the Estimates
  
  ## X1:XT, Lambda, Pmatrix
  
  # homes
  states_est <- matrix(data = rep(NA, Time), nrow = Time, ncol = 1)
  st_rates_est <- matrix(data = NA, nrow = n + 1, ncol = 1)
  st_ptm_est <- matrix(data = rep(NA, n * n), nrow = n * n, ncol = 1)
  
  # estimation
  source("http://www.stat.psu.edu/~mharan/batchmeans.R")
  
  
  for (t in 1:Time) {
    states_est[t, 1] <- mean(states_param[t, ])
  }
  
  
  st_rates_high <- st_rates_scale[2, ] + st_rates_scale[1, ]
  
  st_rates_est <- apply(rbind(st_rates_scale[1, ], st_rates_high), 1, bm)
  st_rates_var <- apply(st_rates_scale, 1, quantile, probs = c(0.025, 
    0.975), na.rm = TRUE)
  
  st_ptm_est <- apply(st_ptm_param, 1, bm)
  st_ptm_var <- apply(st_ptm_param, 1, quantile, probs = c(0.025, 0.975), na.rm = T)
  
  sum.it <- 0 
  
  for (i in 1:n_mcmc) {
    for (t in 1:Time) {
      sum.it <- sum.it +  (osa_param[t, i] - data[t])^2
    }
  }
  
  MSPE.1SA <- 1/n_mcmc * 1/Time * sum.it 
  
  # plot the estimation runs.
  col <- c("#120d08", "#bc5356", "#538bbc", "#53bc84")
  
  if (fig_save == TRUE) {
    jpeg(file = paste(fig_path, fig_name, ".diagnostics" , theta[1], ".jpg", sep = ""))
    
  }
  
  # lambda
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0) + 1, mar = c(1, 1, 1, 
    1) + 3)
  
  
  plot(0, 0, xlab = "MCMC Runs", ylab = "Lambda (scaled per second)", 
    ylim = c(0, max(st_rates_high)), xlim = c(0, n_mcmc), type = "n", 
    cex.lab = 1)
  lines(1:n_mcmc, st_rates_scale[1, ], col = col[1])
  lines(1:n_mcmc, st_rates_high, col = col[2])
  
  
  # for(i in 1:n){ lines(1:n_mcmc, st_rates_scale[i, ], col = col[i])
  # }
  
  
  # P
  plot(0, 0, xlab = "MCMC Runs", ylab = "P", ylim = c(0, max(st_ptm_param)), 
    xlim = c(0, n_mcmc), type = "n", cex.lab = 1)
  for (i in 1:(n * n)) {
    lines(1:n_mcmc, st_ptm_param[i, ], col = col[i])
  }
  
  # Single X
  X <- states_param[sample(1:Time, 1), ]
  plot(0, 0, xlab = "MCMC Runs", ylab = "Single X", ylim = c(0, 
    max(X)), xlim = c(0, n_mcmc), type = "n", cex.lab = 1)
  lines(1:n_mcmc, X, col = col[2])
  
  # States over time
  plot(round(states_est), type = "l", lwd = 3, cex.lab = 1, col = col[1])
  
  title(main = title, outer = T)
  
  if (fig_save == TRUE) {
    dev.off()
  }
  ######################################################### Fancy Plots with Background Colors
  if (fig_save == TRUE) {
    jpeg(file = paste(fig_path, fig_name, ".states" , theta[1], ".jpg", sep = ""))
    
  }
  
  par(mfrow = c(1, 1))
  
  
  if (length(unique(location)) == 1) {
    
    ## High Density - 4 Hours

    plot(start, 1:int.num, main = plot_title, xlab = "Seconds", ylab = "Cumulative
      Interaction Count", 
      xlim = c(0, maxtime))
    states <- states_est  #from code above
    rr <- rle(states[, 1])
    rr$values <- round(rr$values, digits = 0)
    embedded.chain <- rr$values
    cs <- c(0, cumsum(rr$lengths)) * delta_t - delta_t
    cols <- c("#bc535644", "#538bbc44")
    
    for (j in 1:length(embedded.chain)) {
      rect(cs[j], 0, cs[j + 1], int.num, col = cols[embedded.chain[j]], 
        density = NA)
    }
    

    points(start, 1:int.num, main = plot_title, xlab = "Seconds", ylab = "Cumulative 
      Interaction Count", 
      xlim = c(0, maxtime))
  } else {
    # Low Density - 4 Hours
    
    # if (chamber == "queen") {
    #   start <- start_1
    #   int.num <- length(start)
    #   
    # } else {
    #   start <- start_4
    #   int.num <- length(start)
    #   
    # }
    

    plot(start, 1:int.num, main = plot_title, xlab = "Seconds", ylab = "Cumulative 
         Interaction Count", 
      xlim = c(0, maxtime))
    states <- states_est
    rr <- rle(states[, 1])
    rr$values <- round(rr$values, digits = 0)
    embedded.chain <- rr$values
    cs <- c(0, cumsum(rr$lengths)) * delta_t - delta_t
    cols <- c("#bc535644", "#538bbc44")
    for (j in 1:length(embedded.chain)) {
      rect(cs[j], 0, cs[j + 1], int.num, col = cols[embedded.chain[j]], 
        density = NA)
    }

    points(start, 1:int.num, main = plot_title, xlab = "Seconds", ylab = "Cumulative 
      Interaction Count", 
      xlim = c(0, maxtime))
    
  }
  
  
  if(fig_save == TRUE){
    dev.off()
  }
  
  list(states_est = states_est, st_rates_est = st_rates_est, 
       st_ptm_est = st_ptm_est, P.run = st_ptm_param, MSPE = MSPE.1SA)
  
}


#' Penalized Discrete Time MCMC Estimation for Trophallaxis data
#' 
#' The purpose of this function is to find MCMC generated 
#' estimates of (1) - the state (X_t = high/low troph rates) of 
#' the colony at time t (2) - the specific rates of interaction 
#' (lambda) of each state (3) - the specific rates of state 
#' switching (gamma) for the ants (3a) - the probability of 
#' moving from one state to another (P matrix) generated from 
#' gamma values.
#' 
#' @param starts_data Data file of number of trophallaxis 
#'   interactions starting per delta_t time interval.
#' @param ant_file Data file of ant trophallaxis information, 
#'   including Location and start times.
#' @param hours Hours of ant observations
#' @param title Title of diagnositic plots used to check model 
#'   results
#' @param a Hyperprior for rates of starting interactions during 
#'   "low" state. (lambda_L)
#' @param b Hyperprior for rates of starting interactions during 
#'   "low" state. (lambda_L)
#' @param c Hyperprior for rates of starting interactions during 
#'   "high" state. (lambda_H)
#' @param d Hyperprior for rates of starting interactions during 
#'   "high" state. (lambda_H)
#' @param tau Tuning parameter for the proposal function.
#' @param penalty Penalty parameter for the switching rates 
#'   (gamma in model write- up).
#' @param states Number of unobserved states, generally 2. 
#'   ("low", "high")
#' @param n_mcmc Number of MCMC iteractions.
#' @param delta_t Time segemnts the start data is binned into.
#' @param start Starting values for the chains.
#' @param fig_save If "TRUE", plots will saved in path proved as
#'   .jpeg files
#' @param fig_path Path needed to sent plot figures to folder. 
#' @param fig_name Base name of plot files to be saved. 
#'   
#' @return  (1) - estimates of X, lambda, gamma, P (2) - 2x2 
#'   visual of estimates over time (runs) (3) - color block state
#'   switching graph
#' @export
#' @examples
#' theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 
#' out.high = DT.pen_mcmc.troph()
#' 
#' 


DT_pen_mcmc <- function(penalty, starts_data, states, ant_file, chamber, hours, 
  a, b, c, d, tau, tau.pen, n_mcmc, delta_t, start, fig_save, fig_path, fig_name) {
  
  # starting values - mostly to keep this all in one place to easily
  # check
  X.start <- start$X  #same number of values as number of delta_t
  lambda.start <- start$lambda  #two values
  gamma.start <- start$gamma  #two values
  
  # any changes to data
  
  data <- starts_data
  Time <- length(data)
  n <- states
  delta <- rep(1/n, n)
  
  # needed for final graphic
  location <- ant_file$Location
  start <- ant_file$start_time
  
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
  
  # build homes for chains
  
  ## Build Homes for X(1:T), lambda(1:n) vector , and P(nXn), 
  ## gamma (nxn) matrices
  
  states_param <- matrix(data = rep(NA, Time * n_mcmc), nrow = Time, ncol = n_mcmc, 
    byrow = T)

  osa_param <- matrix(NA, Time, n_mcmc, T)
  
  st_rates_param <- matrix(data = NA, nrow = n + 1, ncol = n_mcmc, 
    byrow = T)
  row.names(st_rates_param) <- c("trop.rate.low", "trop.rate.change", 
    "trop.rate.high")
  
  st_ptm_param <- matrix(data = NA, nrow = n * n, ncol = n_mcmc, 
    byrow = T)
  row.names(st_ptm_param) <- c("LL", "LH", "HL", "HH")
  
  
  # note that we're just collecting the off diagonal (non zero) values, needed to calculate P
  switch_rate_param <- matrix(NA, nrow = n, ncol = n_mcmc, byrow = T)
  
  # probability home for generating states_param,
  gam <- matrix(NA, nrow = Time, ncol = n, byrow = T)
  
  starts_low <- matrix(NA, nrow = Time, ncol = n_mcmc)
  starts_high <- matrix(NA, nrow = Time, ncol = n_mcmc)
  
  # penalty.param = matrix(NA, 1, n_mcmc, T)
  
  ## Initialize parameters
  
  states_param[, 1] <- X.start
  
  osa_param[, 1] <- data
  
  st_rates_param[1, 1] <- lambda.start[1]  #lambda low
  st_rates_param[2, 1] <- lambda.start[2] - lambda.start[1]  #change in lambda
  
  st_rates_param[3, 1] <- st_rates_param[1, 1] + st_rates_param[2, 1]  #lambda high, not needed, just a reminder  
  
  switch_rate_param[1, 1] <- gamma.start[1]  #LH
  switch_rate_param[2, 1] <- gamma.start[2]  #HL 
  
  ptm_matrix <- matrix(NA, nrow = n, ncol = n, byrow = T)
  
  
  ptm_matrix[1, 2] <- switch_rate_param[1, 1] * exp(-switch_rate_param[1, 1] * delta_t) #not normalized because within this context will not be greater than one!
  
  ptm_matrix[1, 1] <- 1 - ptm_matrix[1, 2]
  
  ptm_matrix[2, 1] <- switch_rate_param[2, 1] * exp(-switch_rate_param[2, 1] * delta_t) 
  
  ptm_matrix[2, 2] <- 1 - ptm_matrix[1, 2]
  
  st_ptm_param[, 1] <- as.vector(t(ptm_matrix))
  
  
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

  # penalty.param[1, 1] = penalty

    # log likelihood
  
  log.fullcond <- function(st_ptm_param, params, states_param, penalty) {
    
    ptm_matrix <- matrix(data = st_ptm_param, nrow = n, ncol = n, byrow = T)
    
    sumX <- 0
    
    for (t in 2:length(data)) {
      
      sumX <- sumX + log(ptm_matrix[states_param[t - 1, l - 1], states_param[t, 
        l - 1]])
    }
    
    loglike <- sumX - # log(penalty) -
      (1/(2 * penalty) * (params[1]^2 + params[2]^2))
    
    return(loglike)
  }
  
  accept <- 0
  sigma <- tau
  
  for (l in 2:n_mcmc) {
    
    # print out every 100 iterations completed
    if (l %% 100 == 0) 
      cat(paste("iteration", l, "complete\n"))
    
    
    # MH updates - want to propose/accept/reject gammaLH and gammaHL
    
    # adaptive tuning parameter 
    # if (l < n_mcmc/2 & l %% 100 == 0) { 
    #   sigma <- c(sigma,  (2.38^2 / 2) * var(log(t(switch_rate_param[, 1:(l - 1)]))))
    #   tau <- matrix(tail(sigma, n=4), 2, 2) }
    
    #propose switch rates for LH and HL (gamma_LH, gamma_HL)
    proposal <- rmvnorm(n = 1, mean = log(switch_rate_param[, l - 1]), sigma = tau)
    
    # proposal.pen = rmvnorm(n = 1, mean = log(penalty.param[, l - 1]),
    # sigma = tau.pen)
    
    # unlog
    theta.star <- exp(proposal)
    # theta.star.pen = exp(proposal.pen)
    
    
    
    ## Need to take the gamma values we've proposed and caluate the P
    ## matrix for probability of 'jumping' between states
    
    ptm_matrix[1, 2] <- theta.star[1] * exp(-theta.star[1] * delta_t)
    
    ptm_matrix[1, 1] <- 1 - ptm_matrix[1, 2]
    
    ptm_matrix[2, 1] <- theta.star[2] * exp(-theta.star[2] * delta_t)
    
    ptm_matrix[2, 2] <- 1 - ptm_matrix[2, 1]
    
    P.star <- as.vector(t(ptm_matrix))
    

    # calculate probability
    MHprob <- exp(log.fullcond(P.star, theta.star, states_param, penalty) - 
        log.fullcond(st_ptm_param[, l - 1], switch_rate_param[, l - 1], states_param, 
          penalty))
    
    if (is.finite(MHprob) == FALSE) {
      MHprob <- 0
    }
    
    
    # accept/reject
    
    if (runif(1) < MHprob) {
      accept <- accept + 1
      switch_rate_param[, l] <- theta.star
      st_ptm_param[, l] <- as.vector(t(P.star))
      
    } else {
      switch_rate_param[, l] <- switch_rate_param[, l - 1]
      st_ptm_param[, l] <- st_ptm_param[, l - 1]
    }
    
    
    
    # #calculate probability MHprob.pen = exp(log.fullcond(st_ptm_param[,
    # l-1], switch_rate_param[, l-1], states_param, theta.star.pen) -
    # log.fullcond(st_ptm_param[, l - 1], switch_rate_param[, l-1], states_param,
    # penalty.param[ , l - 1])) if(is.finite(MHprob.pen) ==
    # FALSE){MHprob.pen = 0} #accept/reject if(runif(1) < MHprob.pen){ #
    # accept = accept + 1 penalty.param[, l] = theta.star.pen }else{
    # penalty.param[, l] = penalty.param[, l - 1] }
    
    
    
    # gibbs updates
    
    ## X Values over time
    
    ## X Parameters
    
    ptm_matrix <- matrix(data = c(st_ptm_param[, l]), nrow = n, ncol = n, 
      byrow = T)
    
    m <- matrix(data = 0, nrow = 2, ncol = 2)
    rownames(m) <- c("low", "high")
    colnames(m) <- c("low", "high")
    # number states going from i to j, refreshes every run
    
    st_rate_low <- st_rates_param[1, l - 1]
    st_rate_high <- st_rate_low + st_rates_param[2, l - 1]
    
    ## X Parameters
    
    gam[1, 1] <- st_rate_low^data[1] * exp(-st_rate_low) * delta[1] * 
      ptm_matrix[1, states_param[2, l - 1]]
    
    gam[1, 2] <- st_rate_high^data[1] * exp(-st_rate_high) * delta[1] * 
      ptm_matrix[1, states_param[2, l - 1]]
    
    
    
    states_param[1, l] <- sample(x = (1:n), size = 1, prob = gam[1, ])
    
    ################## states_param[1, l] = X.start[1]
    
    m[states_param[1, l], states_param[1, l]] <- m[states_param[1, l], states_param[1, 
      l]] + 1
    
    osa_param[1, l] <- st_rate_low * ptm_matrix[states_param[1, l], 1] + 
      st_rate_high * ptm_matrix[states_param[1, l], 2]
    
    for (t in 2:(Time - 1)) {
      
      gam[t, 1] <- st_rate_low^data[t] * exp(-st_rate_low) * ptm_matrix[states_param[t - 
          1, l - 1], 1] * ptm_matrix[1, states_param[t + 1, l - 1]]
      
      gam[t, 2] <- st_rate_high^data[t] * exp(-st_rate_high) * ptm_matrix[states_param[t - 
          1, l - 1], 2] * ptm_matrix[2, states_param[t + 1, l - 1]]
      
      
      
      states_param[t, l] <- sample(x = (1:n), 1, prob = gam[t, ])
      
      ################## states_param[t, l] = X.start[t]
      
      m[states_param[t - 1, l], states_param[t, l]] <- m[states_param[t - 1, l], 
        states_param[t, l]] + 1
      
      osa_param[t, l] <- st_rate_low * ptm_matrix[states_param[t, l], 
        1] + st_rate_high * ptm_matrix[states_param[t, l], 2]
      
    }
    
    gam[Time, 1] <- st_rate_low^data[Time] * exp(-st_rate_low) *
      ptm_matrix[states_param[Time - 1, l - 1], 1]
    
    gam[Time, 2] <- st_rate_high^data[Time] * exp(-st_rate_high) * 
      ptm_matrix[states_param[Time - 1, l - 1], 2]
    
    
    states_param[Time, l] <- sample(x = 1:n, 1, prob = gam[Time, ])
    
    ################## states_param[Time, l] = X.start[Time]
    
    m[states_param[Time - 1, l], states_param[Time, l]] <- m[states_param[Time - 
        1, l], states_param[Time, l]] + 1
   
    osa_param[Time, l] <- st_rate_low * ptm_matrix[states_param[Time, l], 
      1] + st_rate_high * ptm_matrix[states_param[Time, l], 2]
    
   
   # Data augmentation step - split data into "low" and "high" (N_t <- N_Lt + N_Ht I{X_t = H})
    
    for (t in 1:Time) {
      if (states_param[t, l] == 1) {
        starts_low[t, l] <- data[t]
        starts_high[t, l] <- 0
      } else {
        split <- rmultinom(1, size = data[t], prob = c(st_rates_param[1, 
          l - 1], st_rates_param[2, l - 1]))
        starts_low[t, l] <- split[1]
        starts_high[t, l] <- split[2]
      }
      
    }
    
   ## Lambda Parameters
    
    st_rates_param[1, l] <- rgamma(n = 1, shape = sum(starts_low[, 
      l]) + a, rate = Time + b)
    
    st_rates_param[2, l] <- rgamma(n = 1, shape = sum(starts_high[which(states_param[, 
      l] == 2)]) + c, rate = sum(m[2, ]) + d)
  
  }
  
  
  # estimation
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
  
  sum.it <- 0
  
  for (i in 1:n_mcmc) {
    for (t in 1:Time) {
      sum.it <- sum.it + (osa_param[t, i] - data[t])^2
    }
  }
  
  MSPE.1SA <- 1/n_mcmc * 1/Time * sum.it
  
  
  
  
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
  
  ######################################################### Fancy Plots with Background Colors
  if (fig_save == TRUE) {
    jpeg( file = paste(fig_path, fig_name, round(penalty, 11), ".states", ".jpg", sep = ""))
  }
  
  par(mfrow = c(1, 1))
  
  
  if (length(unique(location)) == 1) {
    
    ## High Density - 4 Hours
    plot(start, 1:int.num, xlab = "Seconds", ylab = "Cumulative Interaction Count", 
      xlim = c(0, maxtime))
    states <- X.est  #from code above
    rr <- rle(states[, 1])
    rr$values <- round(rr$values, digits = 0)
    embedded.chain <- rr$values
    cs <- c(0, cumsum(rr$lengths)) * delta_t - delta_t
    cols <- c("#bc535644", "#538bbc44")
    for (j in 1:length(embedded.chain)) {
      rect(cs[j], 0, cs[j + 1], int.num, col = cols[embedded.chain[j]], 
        density = NA)
      
    }
    points(start, 1:int.num, xlab = "Seconds", ylab = "Cumulative Interaction Count", 
      xlim = c(0, maxtime))
  } else {
    # Low Density - 4 Hours
    
    plot(start, 1:int.num, xlab = "Seconds", ylab = "Cumulative Interaction Count", 
      xlim = c(0, maxtime))
    states <- X.est
    rr <- rle(states[, 1])
    rr$values <- round(rr$values, digits = 0)
    embedded.chain <- rr$values
    cs <- c(0, cumsum(rr$lengths)) * delta_t - delta_t
    cols <- c("#bc535644", "#538bbc44")
    for (j in 1:length(embedded.chain)) {
      rect(cs[j], 0, cs[j + 1], int.num, col = cols[embedded.chain[j]], 
        density = NA)
    }
    points(start, 1:int.num, xlab = "Seconds", ylab = "Cumulative Interaction Count", 
      xlim = c(0, maxtime))
  }
  
  if (fig_save == TRUE) {
    dev.off()
  }
  
  
  list(X.est = X.est, st_rates_est = lambda.est, sw_rates_est = gamma.est, 
    st_ptm_est = P.est, P.run = st_ptm_param, MSPE = MSPE.1SA, accept = accept, sigma = sigma)
  
}


#' Penalized Discrete Time MCMC Estimation function for 
#' Trophallaxis data with Covariates
#' 
#' The purpose of this function is to find MCMC generated 
#' estimates of (1) - the state (X_t = high/low troph rates) of 
#' the colony at time t (2) - the specific rates of interaction 
#' (lambda) of each state (3) - the specific rates of state 
#' switching (gamma) for the ants (3a) - the probability of 
#' moving from one state to another (P matrix) generated from 
#' gamma values.
#' 
#' @param starts_data Data file of number of trophallais 
#'   interactions starting per delta_t time increments.
#' @param ant_file Data file of ant trophallaxis interaction 
#'   information, including Location and start/end times.
#' @param title Title of diagnostic plots used to check model 
#'   results.
#' @param start  Starting values for the chains of parameters to 
#'   be estimated.
#' @param a Hyperprior for rates of starting interactions during 
#'   "low" state. (lambda_L)
#' @param b Hyperprior for rates of starting interactions during 
#'   "low" state. (lambda_L)
#' @param c Hyperprior for rates of starting interactions during 
#'   "high" state. (lambda_H)
#' @param d Hyperprior for rates of starting interactions during 
#'   "high" state. (lambda_H)
#' @param hours Hours of ant obseravtions reflected in data 
#'   files.
#' @param title Title of diagnostic plots.
#' @param tau Tuning parameter for proposal function.
#' @param penalty Penalty parameter for the switching rates 
#'   (alpha/beta in model write- up).
#' @param covariate Vector of second by second covariate data 
#'   used to estimate parameters.
#' @param states Number of unobserved states, generally 2 ("low",
#'   "high").
#' @param n_mcmc Number of MCMC iterations.
#' @param delta_t Time segments that the starting times are 
#'   binned into.
#' @param fig_save If "TRUE", plots will saved in path proved as 
#'   .jpeg files
#' @param fig_path Path needed to sent plot figures to folder.
#' @param fig_name Base name of plot files to be saved.
#'   
#' @return  (1) - estimates of betas, alpha, X, lambda, gamma, P 
#'   (2) - 2x2 visual of estimates over time (runs) (3) - color 
#'   block state switching graph
#' @export
#' 
#' 
#' 


DT_pencov_mcmc <- function(penalty, covariate, starts_data, states, ant_file, chamber,
                          hours, title, a, b, c, d, tau, #tau.pen,
                          n_mcmc, delta_t, fig_save, fig_path, fig_name, start){
  
  #starting values - mostly to keep this all in one place to easily check
 
  X.start <- start$X #same number of values as number of seconds (hours * 60 * 60)
  lambda.start <- start$lambda #two values
  betas.start <- start$alpha.beta #four values (for now)
  
  #any changes to data
  
  data <- starts_data
  Time <- length(data)
  n <- states
  delta <- rep(1/n, n)
  
  #needed for final graphic 
  location <- ant_file$Location 
  start <- ant_file$start_time
  
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
  
  #build homes for chains
  
  ## Build Homes for X(1:T), lambda(1:n) vector , and P(nXn),  gamma (nxn) matrices
  
  states_param <- matrix(data = rep(NA, Time * n_mcmc), nrow = Time, 
    ncol = n_mcmc, byrow = T)
  
  osa_param <- matrix(NA, Time, n_mcmc, T)
  
  st_rates_param <- matrix(NA, nrow = n + 1, ncol = n_mcmc, byrow = T)
  row.names(st_rates_param) <- c("trop.rate.low", "trop.rate.change", 
    "trop.rate.high")
  
  ## P matrix now varies over time, need new homes
  
  ptm_11_param <- matrix(NA, nrow = Time, ncol = n_mcmc)
  ptm_12_param <- matrix(NA, nrow = Time, ncol = n_mcmc)
  ptm_21_param <- matrix(NA, nrow = Time, ncol = n_mcmc)
  ptm_22_param <- matrix(NA, nrow = Time, ncol = n_mcmc)
  
  switch_rates_param <- matrix(data = NA, 
    nrow = 4, #NOT GENERALIZED 
    ncol = n_mcmc)
  #if one covariate this is two by two
  
  #note that we're just collecting the off diagonal (non zero) values 
  #needed to calculate P
  #Since gammas now vary second-to-second based on covariate(s) 
  #we've split this into two separate variables
  
  switch_rate_calc_LH <- matrix(NA, nrow = Time, ncol = n_mcmc, 
    byrow = T)
  
  switch_rate_calc_HL <- matrix(NA, Time, n_mcmc)
  
  
  
  #probability home for generating states_param, 
  gam <- matrix(NA, nrow = Time, ncol = n, 
    byrow = T)
  
  starts_low <- matrix(NA, nrow = Time, ncol = n_mcmc)
  starts_high <- matrix(NA, nrow = Time, ncol = n_mcmc)

  # penalty.param = matrix(NA, 1, n_mcmc, T) #only used when penalty is not fixed
  
  ## Initialize parameters
  
  states_param[,1] <- X.start
  
  osa_param[, 1] <- data
  
  st_rates_param[1, 1] <- lambda.start[1] #lambda low
  st_rates_param[2, 1] <- lambda.start[2]# - lambda.start[1] #change in lambda
  
  #st_rate_high <- st_rates_param[1,1] + st_rates_param[2, 1] #lambda high, not needed, just a reminder  
  
  switch_rates_param[1, 1] <- betas.start[1] #exp{beta_0}
  switch_rates_param[2, 1] <- betas.start[2] #beta_1
  
  switch_rates_param[3, 1] <- betas.start[3] #exp{alpha_0}
  switch_rates_param[4, 1] <- betas.start[4] #alpha_1
  
  ## IN FUTURE NEED TO GENERALIZE ABOVE CODE FOR ANY NUMBER OF COVARIATES!
  ## Current Meridith apologizes to Future Meridith, but this is on Her.
  
  #gammas now function of betas (and alpha?)
  switch_rate_calc_LH[1, 1] <- switch_rates_param[1,1] * exp(switch_rates_param[2,1] * covariate[1]) #LH #covariate[1]??? CHECK
  switch_rate_calc_HL[1, 1] <- switch_rates_param[3,1] * exp(switch_rates_param[4,1] * covariate[1]) #HL 
  
  #Don't forget that now Pmatrix and Gamma values change over time
  P.matrix <- matrix(NA, nrow = n, ncol = n, byrow = T) 
  
  
  P.matrix[1, 2] <- switch_rate_calc_LH[1, 1] * exp(-switch_rate_calc_LH[1, 1] * delta_t) 
  P.matrix[1, 1] <- 1 - P.matrix[1,2]
  P.matrix[2, 1] <- switch_rate_calc_HL[1, 1] * exp(-switch_rate_calc_HL[1, 1] * delta_t) 
  P.matrix[2, 2] <- 1 - P.matrix[2,1]
  
  #need to send each start value to full P variables (same for all runs)
  
  ptm_11_param[, 1] <- P.matrix[1, 1]
  ptm_12_param[, 1] <- P.matrix[1, 2]
  ptm_21_param[, 1] <- P.matrix[2, 1]
  ptm_22_param[, 1] <- P.matrix[2, 2]
  
  
  
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
  

  # penalty.param[1, 1] = penalty
  
  #log likelihood
  
  log.fullcond <- function(P.param, params, states_param, penalty){
    
    #LogFullCond: [X|P(gamma)] + [e^beta_0] + [betas] + [e^alpha_b] + [alphas]
    #Note that the P.matrix varies over time
    # and is represented as a matrix now 
    # P.param <- [P11, P12, P.21, P.22] are column heads
    
    sumX <- 0
    
    for (t in 2:length(data)) {
      
      #pull off P matrix for each time
      P.matrix <- matrix(data = P.param[t - 1, ], nrow = n, ncol = n, byrow = T)
      
      sumX <- sumX + log(P.matrix[states_param[t - 1, l - 1], states_param[t, l - 1]])
    }
    
    loglike <- sumX -
      1/penalty * (params[1]^2)  + log(params[1]) -  #e^beta_0 
      1/penalty * (params[3]^2) + log(params[3]) -  #e^alpha_0
      1/penalty * (params[2]^2) - #betas
      1/penalty * (params[4]^2)  #alphas
    return(loglike)
  }
  
  accept <- 0
  
  
  for (l in 2:n_mcmc) {
    
    # print out every 100 iterations completed
    if ( l %% 100 == 0 ) cat(paste("iteration", l, "complete\n")) 
    
    
    #MH updates - want to propose/accept/reject gammaLH and gammaHL
    
    #adaptive tuning parameter
    # if(l < n_mcmc/2 & l %% 100 == 0){
    # 
    #   sigma <- (2.38^2 / 2) * var(log(t(gamma.param[, 1:(l - 1)])))
    #   tau <- sigma
    # }
    
    #proposing - beta.0, beta.1, alpha_0, alpha_1
    
    #proposal.log <- rmvnorm(n = 1, mean = log(switch_rates_param[c[1, 3], l - 1]), sigma = tau.exp)
    proposal_cov <- rmvnorm(n = 1, mean = c(log(switch_rates_param[1, l - 1]), 
                          switch_rates_param[2, l - 1], 
                          log(switch_rates_param[3, l - 1]), 
                          switch_rates_param[4, l - 1]), 
                          sigma = tau)
   
     # proposal.pen <- rmvnorm(n = 1, mean = log(penalty.param[, l - 1]), sigma = tau.pen)
    
    #unlog 
    proposal <- c(exp(proposal_cov[1]), proposal_cov[2], #beta_0, beta_1
      exp(proposal_cov[3]), proposal_cov[4])
    
    #alpha/betas <- gamma values   
    
    #Note order is not changed as above in commented out line
    sw_rate_star_LH_t <- proposal[1] * exp(proposal[2] * exp(-covariate)) #LH #covariate[1]??? CHECK
    sw_rate_star_HL_t <- proposal[3] * exp(proposal[4] * exp(-covariate)) #HL 
    
    
    
    #uncomment to hold gamma proposal fixed at true values - NEEDS UPDATING
    #
    # theta.star <- c(.005, .005)
    # 
    
    
    
    ## Need to take the gamma values 
    ## we've proposed and caluate the 
    ## P matrix for probability of 'jumping' 
    ## between states (at each second)
    ## 
    P.12.star <- sw_rate_star_LH_t * exp(-sw_rate_star_LH_t * delta_t) 
     
    
    
    P.11.star <- 1 - P.12.star
    
    P.21.star <- sw_rate_star_HL_t * exp(-sw_rate_star_HL_t * delta_t) 
     
    
    P.22.star <- 1 - P.21.star
    
    
    #Uncommont to hold P matrix constant - NEEDS UPDATING
    #
    # P.star <- P.param[,1]
    #
    
    ## NEED TO PULL OFF P param values from last run (for each second)
    ## 
    
    P.11.prev <- ptm_11_param[, l - 1]
    P.12.prev <- ptm_12_param[, l - 1]
    P.21.prev <- ptm_21_param[, l - 1]
    P.22.prev <- ptm_22_param[, l - 1]
    
    
    
    P.star <- cbind(P.11.star, P.12.star, P.21.star, P.22.star)
    P.previous <- cbind(P.11.prev, P.12.prev, P.21.prev, P.22.prev)
    
    #calculate probability
    MHprob <- exp(log.fullcond(P.star, proposal, states_param, penalty) -
        log.fullcond(P.previous, switch_rates_param[, l - 1], states_param, penalty))
    
    if (is.finite(MHprob) == FALSE) {MHprob <- 0}
    
    
    #accept/reject 
    
    if (runif(1) < MHprob) {
      accept <- accept + 1
      switch_rates_param[, l] <- proposal
      ptm_11_param[, l] <- P.11.star[1:(hours*60*60)]
      ptm_12_param[, l] <- P.12.star[1:(hours*60*60)]
      ptm_21_param[, l] <- P.21.star[1:(hours*60*60)]
      ptm_22_param[, l] <- P.22.star[1:(hours*60*60)]   
    }else{
      switch_rates_param[, l] <- switch_rates_param[, l - 1]
      ptm_11_param[, l] <- ptm_11_param[, l - 1]
      ptm_12_param[, l] <- ptm_12_param[, l - 1]
      ptm_21_param[, l] <- ptm_21_param[, l - 1]
      ptm_22_param[, l] <- ptm_22_param[, l - 1]  }
    
    
    
    #gibbs updates
    
    ## X Values over time
    
    ##X Parameters
    
    #Pull off just current probabilities (over time)
    #and then combine into matrix
    
    P.param <- cbind(ptm_11_param[, l], ptm_12_param[, l], ptm_21_param[, l], ptm_22_param[, l])
    
    
    m <- matrix(data = 0, nrow = 2, ncol = 2)
    rownames(m) <- c("low", "high")
    colnames(m) <- c("low", "high")
    # number states going from i to j, refreshes every run
    
    st_rate_low <- st_rates_param[1, l - 1]
    st_rate_high <- st_rates_param[2, l - 1] #+st_rate_low
    
    ##X Parameters
    
    P.matrix <- matrix(data = P.param[1, ], nrow = n, ncol = n, byrow = T)
    
    
    gam[1, 1] <- st_rate_low ^ data[1] * exp(-st_rate_low) * 
      delta[1] * P.matrix[1, states_param[2, l - 1]]
    
    gam[1, 2] <- st_rate_high ^ data[1] * exp(-st_rate_high) * 
      delta[1] * P.matrix[1, states_param[2, l - 1]]
    
    
    
    states_param[1, l] <- sample(x = (1:n), size = 1, prob = gam[1,])
    
    ##################
    # states_param[1, l] <- X.start[1] 
    ##################
    
    m[states_param[1, l], states_param[1, l]] <- m[states_param[1, l], states_param[1, l]] + 1
    
    osa_param[1, l] <- st_rate_low * P.matrix[states_param[1, l], 1] + 
      st_rate_high * P.matrix[states_param[1, l], 2]
    
    for(t in 2:(Time - 1)){
      P.matrix <- matrix(data = P.param[t, ], nrow = n, ncol = n, byrow = T)
      
      
      gam[t, 1] <- st_rate_low ^ data[t] * exp(-st_rate_low) * 
        P.matrix[states_param[t - 1, l - 1], 1] *
        P.matrix[1, states_param[t + 1, l - 1]]
      
      gam[t, 2] <- st_rate_high ^ data[t] * exp(-st_rate_high) * 
        P.matrix[states_param[t - 1, l - 1], 2] *
        P.matrix[2, states_param[t + 1, l - 1]]
      
      
      
      states_param[t, l] <- sample(x = (1:n), 1,  prob = gam[t, ]) 
      
      ##################
      # states_param[t, l] = X.start[t]
      ##################
      
      m[states_param[t - 1, l], states_param[t, l]] <- m[states_param[t - 1, l], 
        states_param[t,l]] + 1
      
      osa_param[t, l] <- st_rate_low * P.matrix[states_param[t, l], 1] + 
        st_rate_high * P.matrix[states_param[t, l], 2]
      
    }
    
    P.matrix <- matrix(data = P.param[Time, ], nrow = n, ncol = n, byrow = T)
    
    
    gam[Time, 1] <- st_rate_low ^ data[Time] * exp(-st_rate_low) * 
      P.matrix[states_param[Time - 1, l - 1], 1] 
    
    gam[Time, 2] <- st_rate_high ^ data[Time] * exp(-st_rate_high) * 
      P.matrix[states_param[Time - 1, l - 1], 2] 
    
    
    states_param[Time, l] <- sample(x = 1:n, 1,  prob = gam[Time, ])
    
    ##################
    # states_param[Time, l] <- X.start[Time]
    ##################
    
    m[states_param[Time - 1, l], states_param[Time, l]] <- m[states_param[Time - 1, l], 
      states_param[Time, l]] + 1
    
    osa_param[Time, l] <- st_rate_low * P.matrix[states_param[Time, l], 1] + 
      st_rate_high * P.matrix[states_param[Time, l], 2]
    
    # Data augmentation step - split data into "low" and "high" (N_t <- N_Lt + N_Ht I{X_t = H})
    
    for (t in 1:Time) {
      if (states_param[t, l] == 1) {
        starts_low[t, l] <- data[t]
        starts_high[t, l] <- 0
      } else {
        split <- rmultinom(1, size = data[t], prob = c(st_rates_param[1, 
          l - 1], st_rates_param[2, l - 1]))
        starts_low[t, l] <- split[1]
        starts_high[t, l] <- split[2]
      }
      
    }
    
    ## Lambda Parameters
    
    st_rates_param[1, l] <- rgamma(n = 1, shape = sum(starts_low[, 
      l]) + a, rate = Time + b)
    
    st_rates_param[2, l] <- rgamma(n = 1, shape = sum(starts_high[which(states_param[, 
      l] == 2)]) + c, rate = sum(m[2, ]) + d)
  
  }
  
  
  #estimation
  source("http://www.stat.psu.edu/~mharan/batchmeans.R")
  
  st_rate_high <-  st_rates_param[2, ] #+ st_rates_param[1, ]
  
  lambda.est <- apply(rbind(st_rates_param[1, ], st_rate_high), 1, bm) 
  lambda.var <- apply(st_rates_param, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 
  
  # gamma.est <- apply(gamma.param, 1, bm) 
  # gamma.var <- apply(gamma.param, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 
  
  
  P.11.est <- apply(ptm_11_param, 1, bm)
  P.11.var <- apply(ptm_11_param, 1 , quantile, probs = c(0.025, 0.975, na.rm = T))
  
  P.12.est <- apply(ptm_12_param, 1, bm)
  P.12.var <- apply(ptm_12_param, 1 , quantile, probs = c(0.025, 0.975, na.rm = T))
  
  P.21.est <- apply(ptm_21_param, 1, bm)
  P.21.var <- apply(ptm_21_param, 1 , quantile, probs = c(0.025, 0.975, na.rm = T))
  
  P.22.est <- apply(ptm_22_param, 1, bm)
  P.22.var <- apply(ptm_22_param, 1 , quantile, probs = c(0.025, 0.975, na.rm = T))
  
  #Want to create combined estimates over time for visuals/estimation
  X.est <- matrix(NA, Time, 1, T)
  P.11.est <- matrix(NA, Time, 1)
  P.12.est <- matrix(NA, Time, 1)
  P.21.est <- matrix(NA, Time, 1)
  P.22.est <- matrix(NA, Time, 1)
  
  for (t in 1:Time) {
    X.est[t, 1] <- mean(states_param[t, ]) 
    P.11.est[t, 1] <- mean(ptm_11_param[t, ])
    P.12.est[t, 1] <- mean(ptm_12_param[t, ])
    P.21.est[t, 1] <- mean(ptm_21_param[t, ])
    P.22.est[t, 1] <- mean(ptm_22_param[t, ])
  }
  
  sum.it <- 0 
  
  for (i in 1:n_mcmc) {
    for (t in 1:Time) {
      sum.it <- sum.it +  (osa_param[t, i] - data[t])^2
    }
  }
  
  MSPE.1SA <- 1/n_mcmc * 1/Time * sum.it 
  
  
  #visualization
  #
  if (fig_save == TRUE) {
    jpeg(file = paste(fig_path, fig_name, round(penalty, digits = 11), ".diagnostics", ".jpg", sep = ""))
    
  }
  
  #plot the estimation runs.
  
  col <- c("#120d08", "#bc5356", "#538bbc", "#53bc84")
  
  # if(fig_save == T){
  #   pdf(file = paste("./output/", Sys.time(), ".pdf", sep = ""))
  # }
  
  
  
  #Rate Parameters
  plot(0,0,xlab = "MCMC Runs",
    ylab = "Rates",
    ylim = c(0, 2*max(st_rates_param[1:2, ])), 
    xlim = c(0,n_mcmc), 
    type = "n",
    cex.lab = 1)
  lines(1:n_mcmc, (st_rates_param[1, ]), col = col[1])
  lines(1:n_mcmc,  st_rates_param[2, ], col = col[2])
  
  par(mfrow = c(2,2),
    oma = c(0,0,2,0) + 1,
    mar = c(1,1,1,1) + 3)
  
  
  choose <- sample(1:Time, 1) #Choosing one second to explore through all the runs
  P.11 <- ptm_11_param[choose, ]
  P.21 <- ptm_21_param[choose, ]
  P.12 <- ptm_12_param[choose, ]
  P.22 <- ptm_22_param[choose, ]
  
  P <- cbind(P.11, P.21, P.12, P.22)
  
  #Single P
  plot(0,0,xlab = "MCMC Runs",
    ylab = "Single Second Probability Matrix for State Switching",
    ylim = c(0, 1), xlim = c(0,n_mcmc),
    type = "n", cex.lab = 1)
  for (i in 1:(4)) {
    lines(1:n_mcmc, P[, i], col = col[i])
  }
  
  #Probabilities over time
  P.est <- cbind(P.11.est, P.12.est, P.21.est, P.22.est)
  plot(0,0,xlab = "MCMC Runs",
    ylab = "Combined Probability Matrix for State Switching",
    ylim = c(0, 1), xlim = c(0,n_mcmc),
    type = "n", cex.lab = 1)
  for (i in 1:(4)) {
    lines(1:Time, P.est[, i], col = col[i])
  }
  
  
  #Single X
  X <- states_param[choose, ]
  plot(0, 0, xlab = "MCMC Runs", ylab = "Single X", ylim = c(0,max(X)), 
    xlim = c(0,n_mcmc), type = "n", cex.lab = 1)
  lines(1:n_mcmc, X, col = col[4])
  
  #States over time
  #plot(X.est, type = "l")
  plot(round(X.est), type = "l")
  
  
  title(main = title, outer = T)
  
  if (fig_save == TRUE) {
    dev.off()
  }
  
  #########################################################
  ##
  ## Fancy Plots with Background Colors
  ##
  #########################################################
 
  if (fig_save == TRUE) {
    jpeg(file = paste(fig_path, fig_name, round(penalty, digits = 11), ".states", ".jpg", sep = ""))
  }
  
  par(mfrow = c(1, 1))
  
  
  if (length(unique(location)) == 1) {
    
    ##High Density - 4 Hours
    plot(start, 1:int.num, main = "High", 
      xlab = "delta_t", ylab = "Cumulative Interaction Count", 
      xlim = c(0, maxtime))
    states <- X.est #from code above
    rr <- rle(states[,1])
    rr$values <- round(rr$values, digits = 0)
    embedded.chain <- rr$values
    cs <- c(0,cumsum(rr$lengths))*delta_t - delta_t
    cols <- c('#bc535644','#538bbc44')
    
    for (j in 1:length(embedded.chain)) {
      rect(cs[j],0,cs[j + 1],int.num, 
        col = cols[embedded.chain[j]], border = NA, density = NA)
      
    }
    points(start, 1:int.num, main = "Low", xlab = "delta_t",
      ylab = "Cumulative Interaction Count", 
      xlim = c(0,maxtime))
  }else{
    #Low Density - 4 Hours
    
    plot(start, 1:int.num, main = "Low", xlab = "delta_t",
      ylab = "Cumulative Interaction Count", 
      xlim = c(0,maxtime))
    states <- X.est
    rr <- rle(states[,1])
    rr$values <- round(rr$values, digits = 0)
    embedded.chain <- rr$values
    cs <- c(0,cumsum(rr$lengths))*delta_t - delta_t
    cols <- c('#bc535644','#538bbc44')
    for (j in 1:length(embedded.chain)) {
      rect(cs[j],0,cs[j + 1],int.num, col = cols[embedded.chain[j]] , density = NA)
    }
    points(start, 1:int.num, main = "Low", xlab = "delta_t",
      ylab = "Cumulative Interaction Count", 
      xlim = c(0,maxtime))
  }
  
  if (fig_save == TRUE) {
    dev.off()
  }
  
  
  list(X.est = X.est, st_rates_est = lambda.est, 
    P.est = c(P.11.est, P.12.est, P.21.est, P.22.est), MSPE = MSPE.1SA, accept = accept)
  
}



