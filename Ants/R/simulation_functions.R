#' 
#' Simulation of Ant Trophallaxis Start Times - Simple Model
#' 
#' \code{sim_DT_troph} simulates number of ant interactions
#' starting per second and per delta_t seconds.
#' 
#' This function allows you to simulate data from the following
#' model: X_t ~ Markov Chain ( X_t-1 , P ) [two states] N_t ~
#' Pois (lambda_L + lambda_Change) lambda_L ~ Gam(a, b)
#' lambda_Change ~ Gam(c, d) lambda_H = lambda_l + lambda_Change
#' 
#' @param time_max Number of seconds for which to simulate data.
#' @param delta_t Length of interval to bin interaction and state 
#'   data into.
#' @param start_state Starting trophallaxis rate state (X_1 in model
#'   write - up). Defaults to state 1 (low).
#' @param switch_rate Rate at which the rates at which ants trophallax 
#'   switch (between high and low). Denoted by gamma in model write-up. 
#'   Used to calculate probability transition matrix (P).
#' @param state_tpm Probably transition matrix for switching between 
#'   trophallaxis states (denoted by P in write-up). 
#' @param int_rate Rate of starting interactions. Denoted by lambda in 
#'   model write-up).
#' @param num_locations Number of chambers/locations in which ants are 
#'   found trophallaxis.
#'   
#' @return This function will return the following:
#'   \enumerate{
#'     \item inter_persec: (0, ...) 'observed' number of
#'       interactions per 1 second
#'     \item state: (1 = low, 2 = high) unobserved two state process
#'     \item cumu_inter: (0, ...) cumulative count of interactions
#'     \item bin_inter: Trophallaxis interactions binned into smaller 
#'       intervals determined by sum over delta_t
#'     \item bin_state: Trophallaxis rate state binned into smaller 
#'       intervals, determined by average over delta_t
#'     \item bin_sec: (0, ...) time (in seconds) binned by delta_t
#'     \item start_time: vector of interaction start times
#'     \item location: location of each ant interaction
#'     \item Visuals: 1x1 visual of cumlative counts over time, colored 
#'       by state, separated by location if applicable.
#'   }
#'  
#'   
#'   
#' @keywords simulation, HMM, Poisson Process
#' @export
#' @examples 
#' P <- matrix(c(.99, .01, .01, .99), nrow = 2, byrow = T)
#' lambda <- k <- c(1, 4)
#' delta_t <- 5 
#' sim <- sim.DT.troph(time_max = 3600, delta_t = delta_t, 
#'                     start_state = 1, switch_rate = 0, state_tpm = P, 
#'                     int_rate = lambda, num_locations = 1)
#' 

sim_DT_troph <- function(time_max, delta_t, start_state = 1, 
                        switch_rate = c(0, 0), state_tpm, int_rate,
                        num_locations = 1) {
    
  #if given switch rates (gamma in model) use them to calculate transition 
  #probability matrices 
  
   if (switch_rate[1] == 0) {
        state_tpm <- state_tpm
    } else {
        state_tpm[1, 2] <- switch_rate[1] * exp(-switch_rate[1] * 1)/(switch_rate[1]/(-1 + exp(switch_rate[1])))
        state_tpm[1, 1] <- 1 - state_tpm[1, 2]
        state_tpm[2, 1] <- switch_rate[2] * exp(-switch_rate[2] * 1)/(switch_rate[2]/(-1 + exp(switch_rate[2])))
        state_tpm[2, 2] <- 1 - state_tpm[1, 2]
    }
    
    
    #create homes for simulated states and interactions 
  
    state <- rep(NA, time_max)
    interactions <- rep(NA, time_max)
    state[1] <- start_state
    interactions[1] <- 0
    
    for (t in 2:time_max) {
        
        
        ## sample latent state (low, high)
        state[t] <- sample(1:2, 1, prob = state_tpm[state[t - 1], ])
        
        ## sample observed events (number of interactions started at time t)
        
        interactions[t] <- rpois(1, lambda = int_rate[state[t]] * 1)
        # want to simulate data every second, but then bin by delta_t
        # aftwards
        
    }
    
    t <- (0:(time_max - 1))  #vector of seconds through time_max
     
    
    # bin number of interactions started within delta_t time intervals

    bin_inters <- rep(NA, length(interactions)/delta_t) 
    
    tintv <- 1
    
    
    for (i in 1:length(bin_inters)) {
        
        bin_inters[i] <- sum(interactions[tintv:(tintv + delta_t - 1)])
        tintv <- tintv + delta_t
    }
    
    #bin trophallaxis states within delta_t time intervals
    
    bin_state <- rep(NA, length(bin_inters))
    
    tintv <- 1
    
  
    for (i in 1:length(bin_state)) {
        bin_state[i] <- round(mean(state[tintv:(tintv + delta_t - 1)]))
        tintv <- tintv + delta_t
    }
    
    
    
    start_time <- t[which(interactions >= 1)]
    
    #Visualize simulated ant interaction data - POSSIBLY PULL INTO VISUAL FUNCTION LATER
    
    par(mfrow = c(1, 1))
    
    plot(1:(length(cumsum(interactions))), cumsum(interactions), 
        type = "p", pch = ".", 
        cex = 2, col = state, xlab = "Seconds", 
        ylab = "Cumulative Interactions", 
        main = "Simulated Data")
    
    plot(1:time_max, interactions, type = "l", cex = 2, col = state,
        xlab = "Seconds", ylab = "Number of Interactions", 
        main = "Simulated Data")
    
    #simulate location (1 = queen's chamber, 4 = entrance)
    
    if (num_locations == 1) {
        location <- rep(1, length(start_time))
    } else {
        location <- sample(c(1, 4), size = length(start_time), replace = T, 
            prob <- c(0.5, 0.5))
    }
    
    
    
    list(inter_persec = interactions, state = state, 
        cumu_inter = cumsum(interactions),  
        bin_inter = bin_inters, 
        bin_state = bin_state, 
        bin_sec = (0:(T - 1)) * delta_t, 
        start_time = start_time, Location = location)
    
}
