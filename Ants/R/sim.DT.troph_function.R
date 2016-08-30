#` Simulation of Poisson Proccess HMM
#'
#' This function allows you to simulate data from the following moddel: 
#' X_t ~ Markov Chain ( X_t-1 , P ) [two states]
#' Y_t ~ Pois (lambda_{X_t})
#' @param tmax, delta.t, start.state, P, lambda, num.location
#' @return #Output: (1) - state: (1, 2) unobserved two state process 
#'        (2) - inter.persec: (0, ...) "observed" number of interactions per 1 second
#'        (3) - cumu.inter: (0, ...) cumulative count of interactions
#'        (4) - delta.t
#'        (5) - bin.sec: (0, ...) time (in seconds) binned by delta.t
#'        (6) - 1x1 visual of cumlative counts over time, colored by state, 
#'              separated by location if applicable
#' @keywords simulation, HMM, Poisson Process
#' @export
#' @examples 
#' P = matrix(c(.99, .01, .01, .99), nrow = 2, byrow = T)
#' lambda = k = c(1, 4)
#' delta.t = 1 #needs to be 1, else observations dependent on time
#' sim = sim.DT.troph(7200, delta.t, start.state = 1, P, lambda, num.location = 1)
#' 



sim.DT.troph <- function(tmax, delta.t, start.state = 1 , P, lambda, num.locations = 1){
   T = floor(tmax/delta.t) + 1
  # x=rep(NA,T)
  # y=rep(NA,T)
  # 
  x = rep(NA, tmax)
  y = rep(NA, tmax)
  x[1] = start.state
  y[1] = 0
  
  #for(t in 2:T){
  for(t in 2:tmax){
    ## sample latent state
    x[t]=sample(1:2,1,prob = P[x[t - 1],])
    
    ## sample observed events
    #y[t]=rpois(1,lambda=lambda[x[t]]*delta.t)
    y[t] = rpois(1,lambda = lambda[x[t]]*1)
    #want to simulate data every second, but then bin by delta.t aftwards
    #number of interactions per second cannot exceed 1
  }
  
  t = (0:(tmax - 1)) #vector of seconds through tmax
  
  bin.y =  rep(NA, length(y)/delta.t)
  
        tint = 1 
        
        for(i in 1:length(bin.y)){
          
          bin.y[i] = sum(y[tint:(tint + delta.t - 1)])
          tint = tint + delta.t
        }
        
  bin.x = rep(NA, length(bin.y))
  
      tint = 1     
  
      for(i in 1:length(bin.x)){
        bin.x[i] = round(mean(x[tint:(tint + delta.t - 1)]))
        tint = tint + delta.t
      }
        
  start.time = t[which(y >= 1)] #only works with per second data
 
   par(mfrow = c(1, 1))
  plot(1:(length(cumsum(bin.y))), cumsum(bin.y), 
       type = "p", pch = ".", cex = 2, col = bin.x,
       xlab = "Time", ylab = "Interactions")
  
  if(num.locations == 1){
     Location = rep(1, length(start.time))
  }

  else{
    Location = sample(1:2, size = length(start.time), replace = T, prob = c(.5, .5))
    }
  
  
  
  list(inter.persec = y,state = x, cumu.inter = cumsum(y), delta.t = delta.t,
       bin.inter = bin.y, bin.state = bin.x, 
       bin.sec = (0:(T - 1))*delta.t, start.time = start.time, Location = Location)
}
