# ` Simulation of Poisson Proccess HMM with dynamic P matrix
#'
#' This function allows you to simulate data from the following moddel: 
#' X_t ~ Markov Chain ( X_t-1 , P_t ) [two states]
#' Y_t ~ Pois (lambda_{X_t})
#' 
#' 
#'  P_1,2_t = (exp(beta.0 + beta.1 * covariate[t])) / 
#'            (1 + exp(beta.0 + beta.1 * covariate[i])) 
#' P_1,1_t = 1 - P_1,2_t
#' P_2,1_t = (exp(alpha)) / (1 + exp(alpha))
#' P_2,2_t = 1 - P_2,1_t
#' 
#' @param tmax, start.state, P11, P12, P21, P22, lambda
#' @return #Output: (1) - x: (1, 2) unobserved two state process 
#'        (2) - y: (0, ...) 'observed' number of interactions at time t
#'        (3) - N: (0, ...) cumulative count of interactions
#'        (4) - delta.t
#'        (5) - t: (0, ...) time (in seconds)
#'        (6) - 1x1 visual of cumlative counts over time, colored by state
#' @keywords simulation, HMM, Poisson Process, covariate, 
#' @export 
#' @examples 
#' 
#'
#' lambda = k = c(1, 4)
#' sim = sim.mcmc.dynamP(7200, start.state = 1, P11, P12, P21, P22, lambda)
#' 


sim.mcmc.dynamP <- function(tmax, start.state = 1, P11, P12, P21, 
    P22, lambda) {
    T = tmax
    x = rep(NA, T)
    y = rep(NA, T)
    x[1] = start.state
    y[1] = 0
    for (t in 2:T) {
        P = matrix(NA, 2, 2)
        
        P[1, 2] = P12[t, 1]
        P[1, 1] = P11[t, 1]
        P[2, 2] = P22[t, 1]
        P[2, 1] = P21[t, 1]
        
        ## sample latent state
        x[t] = sample(1:2, 1, prob = P[x[t - 1], ])
        ## sample observed events
        y[t] = rpois(1, lambda = lambda[x[t]])
    }
    par(mfrow = c(1, 1))
    plot((0:(T - 1)), cumsum(y), type = "p", pch = ".", cex = 2, col = x, 
        xlab = "Time", ylab = "Interactions")
    
    list(y = y, x = x, N = cumsum(y), t = (0:(T - 1)))
}
