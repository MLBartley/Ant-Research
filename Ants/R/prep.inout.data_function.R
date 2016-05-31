#' Data Preparation for Ant In/Out Data File
#'
#' The purpose of this function is to take a .csv file with 
#' entrance/exit data and format it for use in mcmc model. 
#' 
#' @param data, delta.t, hours
#' @return (1) List file for binned covariate (forager arrival) data.
#' @export
#' @examples
#' 
#' prep.inout.data(inout.high.4, 60, 4)


prep.inout.data = function(data, delta.t, hours){
  forager.arrivals = sort(data$time)
  forager.arrivals = forager.arrivals[which(duplicated(forager.arrivals) == FALSE)]
  forager.arrivals = c(forager.arrivals, 999999)
  
  ## Time since Arrivals
  
  time = hours * 60 * 60
  covariate = rep(NA, time)
  covariate[1] = 300 #five minutes since arrival
  
  arrival = 1
  
  for(i in 2:time){
    
    if(forager.arrivals[arrival] == i){
      covariate[i] = 0
      arrival = arrival + 1
    }else{
      covariate[i] = covariate[i - 1] + 1
    }
    
  }
  
  #bin covariates into similar chunks as the data
  max.time = time
  
  if(delta.t != 1){
    for(i in 1:length(covariate)){
    c = rep(0, max.time / delta.t)
    mint = 1
    for(t in 1:length(c)){
      c[t] = min(covariate[mint:(mint + delta.t - 1)])
      mint = mint + delta.t
    }
    }
  }
  
  
  list(cov = c)
}