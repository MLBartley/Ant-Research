#' Alternative Data Preparation for Ant Trophallaxis Data File
#'
#' The purpose of this function is to take in .csv files 
#' and apply necissary changes to the data format. Changes include 
#' 
#' 
#' @param data,
#' @return (1) .
#' @export
#' @examples
#'  prep.troph.pairs(high4)
#' 
#' 


prep.troph.pairs = function(data){
  
  hours = ceiling(max(data$end_time) / 60 / 60)
  
  
  num.high = unique(data$Ant_ID)
  ant = length(num.high)  #number of ants
  
  #Embedded Chain and Rate Transitions
  Master = matrix(0, ant, hours * 60 * 60)
  n = 0
  
  for(i in sort(num.high)){
    n = n+1
    high.i = data[which(data$Ant_ID == i), ]
    
    for(l in 1:nrow(high.i)){
      start.i = high.i$start_time[l]
      end.i = high.i$end_time[l]
      Master[n, start.i:end.i] = 1
    }
    
  }
  
  Master.sum = colSums(Master)
  
  cont.time = rle(Master.sum)
  
  EC = cont.time$values # embedded chain
  RT = cont.time$lengths #time in embedded chain state
  
   N = c()
  
  for(i in 1:length(EC)){
    N = c(N, rep(EC[i], RT[i])) # pairs in troph at time t for every second
  }
   
  N.2 = c()
  
  for(i in 1:length(N)){
    if(N[i] %% 2 != 0){
      N.2[i] = N[i] - 1
    }else(N.2[i] = N[i])
  }
   
  ##plot of interactions
  t = cumsum(RT)  
  X = EC
  
  plot(0,0,xlab="t",ylab="State",ylim=c(1,max(X)),xlim=c(0,hours * 60 * 60),type="n")
  lines(x = c(0, t[1]), y = c(X[1], X[1]))
  
  for(i in 2:length(t)){
    lines(x = c(t[i-1], t[i]), y = c(X[i], X[i]))
  }
  
  
 list(EC = EC, RT = RT, sum = cumsum(RT), pairs3 = N, pairs2 = N.2, hours = hours)
  
}