#` Simulation of Poisson Proccess HMM - with forager times covariate
#'
#' This function allows you to simulate data from the following model: 
#' P ~ function of covariates, changes over time
#' X_t ~ Markov Chain ( X_t-1 , P ) [two states]
#' Y_t ~ Pois (lambda_{X_t})
#' @param tmax, delta.t, start.state, alpha.beta, lambda, num.location
#' @return #Output: (1) - state: (1, 2) unobserved two state process 
#'        (2) - inter.persec: (0, ...) "observed" number of interactions per 1 second
#'        (3) - cumu.inter: (0, ...) cumulative count of interactions
#'        (4) - delta.t
#'        (5) - bin.sec: (0, ...) time (in seconds) binned by delta.t
#'        (6) - 1x1 visual of cumlative counts over time, colored by state, 
#'              separated by location if applicable
#' @keywords simulation, HMM, Poisson Process, covariate
#' @export
#' @examples 
#' alpha.beta = c(.005, .0001, .005, .0001)
#' lambda = k = c(1, 4)
#' delta.t = 1 #needs to be 1, else observations dependent on time
#' sim = sim.DT.troph(7200, delta.t, start.state = 1, alpha.beta, lambda, num.location = 1)
#' 



sim.DT.cov.troph <- function(tmax, delta.t, start.state = 1,
                          alpha.beta, lambda, num.locations = 1){
  

  covariate = rep(NA, tmax)

  #based on data on entrance times for 4 hours L/H density, 
  #assuming about 20 ants per hour enter next

  
  time = 0
  num.arrive = 0
  inter.arrive = 0
  arrive.time = 0
  
  while(time < tmax){
    
    wait = rexp(1, rate = .005)
    num.arrive = num.arrive + 1
    inter.arrive = append(inter.arrive, wait)
    arrive.time = append(arrive.time, time + wait)
  
    time = time + wait
  }
  
  
  
  covariate = rep(NA, tmax)
  covariate[1] = 300 #five minutes since arrival
  
  arrival = 1
  
    
    for(i in 2:tmax){
      
      if(round(arrive.time[arrival+1]) == i){ #first 'time' is just zero placeholder
        covariate[i] = 0
        arrival = arrival + 1
      }else{
        covariate[i] = covariate[i - 1] + 1
      }
      
    }
  
 
  
  
  ## IN FUTURE NEED TO GENERALIZE ABOVE CODE FOR ANY NUMBER OF COVARIATES!
  ## Current Meridith apologizes to Future Meridith, but this is on Her.
  
  #gammas now function of betas (and alpha?)
  gamma.star.LH.t = alpha.beta[1] * exp(alpha.beta[2] * covariate) #LH #covariate[1]??? CHECK
  gamma.star.HL.t = alpha.beta[3] * exp(alpha.beta[4] * covariate) #HL 
  
  #Don't forget that now Pmatrix and Gamma values change over time
  #Want to first simulate second by second data then optionally
  #bin into delta.t segments
  
  P.12.star = gamma.star.LH.t * exp(-gamma.star.LH.t * 1) /
    (gamma.star.LH.t / (-1 + exp(gamma.star.LH.t)) )
  
  
  P.11.star = 1 - P.12.star
  
  P.21.star = gamma.star.HL.t * exp(-gamma.star.HL.t * 1) /
    (gamma.star.HL.t / (-1 + exp(gamma.star.HL.t)) )
  
  P.22.star = 1 - P.21.star
  
  P.star = cbind(P.11.star, P.12.star, P.21.star, P.22.star)
  

  
  
 
  x = rep(NA, tmax) #states
  y = rep(NA, tmax) # number of interactions 
  x[1] = start.state
  y[1] = 0
  
  
  for(t in 2:tmax){
    
    P = matrix(P.star[t, ], 2, 2, byrow = T)
    
    ## sample latent state
    x[t]=sample(1:2, 1, prob = P[x[t - 1],])
    
    ## sample observed events
    #y[t]=rpois(1,lambda=lambda[x[t]]*delta.t)
    y[t] = rpois(1,lambda = lambda[x[t]]*1)
    #want to simulate data every second, but then bin by delta.t aftwards
  }
  
  ## AM I SIMULATING ONLY START TIMES OR START AND WAIT TIMES?!!
  
   t = (0:(tmax - 1)) #vector of seconds through tmax
   
if(delta.t > 1){
  
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
}

  
  start.time = t[which(y >= 1)] #only works with per second data
  
   if(num.locations == 1){
    Location = rep(1, length(start.time))
  }
  
  else{
    Location = sample(1:2, size = length(start.time), replace = T, prob = c(.5, .5))
  }
  
  
  par(mfrow = c(1, 1))
  
  plot(1:(length(cumsum(y))), cumsum(y), 
       type = "p", pch = ".", cex = 2, col = x,
       xlab = "Time", ylab = "Interactions", 
       main = "Full Timeline, Cumulative Interactions")
  points(arrive.time, rep(0, length(arrive.time)), pch=8, col=	"#53bc84")
  
  plot(1:tmax, y, 
       type = "l", cex = 2, col = x,
       xlab = "Time", ylab = "Number of Interactions", 
       main = "Full Timeline")
  
  
  
 if(delta.t > 1){
    plot(1:(length(cumsum(bin.y))), cumsum(bin.y), 
       type = "p", pch = ".", cex = 2, col = bin.x,
       xlab = "Time", ylab = "Interactions", main = "Binned Intervals")
 
   list(inter.persec = y,state = x, cumu.inter = cumsum(y), delta.t = delta.t,
       bin.inter = bin.y, bin.state = bin.x, 
       bin.sec = (0:(T - 1))*delta.t, 
       start_time = start.time, Location = Location, covariate = covariate)
  
 }
else{
  list(inter.persec = y,state = x, cumu.inter = cumsum(y),
       delta.t = delta.t, start_time = start.time, Location = Location, covariate = covariate)
}
  
  
 
  
  
  
  
}
