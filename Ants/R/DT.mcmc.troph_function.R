#' Discrete Time MCMC Estimation function for Trophallaxis data
#'
#' The purpose of this function is to find MCMC generated estimates of 
#' (1) - the state (X_t = high/low troph rates) of the colony at time t
#' (2) - the specific rates of interaction (lambda) of each state 
#' (3) - the probability of moving from one state to another (P matrix)
#'
#' @param y.data, ant.file, title, a, b, theta, states, n.mcmc, delta.t
#' @return  (1) - estimates of X, lambda, P
#'          (2) - 2x2 visual of estimates over time (runs)
#' @export
#' @examples
#' theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 
#' out.high = DT.mcmc.troph(data = high.y, title = "High Density", a = 5, b = 2, 
#' theta = theta, states = 2, n.mcmc = 3000)


DT.mcmc.troph = function(y.data, ant.file, title, a, b, c, d,
                      theta, states = 2, n.mcmc, delta.t, hours){
  data = y.data
  Time = length(data)
  n = states
  delta = rep(1/n, n)
  
  #needed for final graphic 
  location = ant.file$Location 
  start = ant.file$start_time
  start = sort(start)
  int.num = length(start)
  maxtime = hours * 60 * 60
  

  #homes
  ## Build Homes for X(1:T), lambda(1:n), and P(nXn) and gam vectors
  
  X.param = matrix(data = rep(NA, Time * n.mcmc), nrow = Time, 
                   ncol = n.mcmc, byrow = T)
  
  lambda.param = matrix(data = rep(NA, n * n.mcmc), nrow = n, 
                        ncol = n.mcmc, byrow = T)
  
  P.param = matrix(data = rep(NA, n * n * n.mcmc), nrow = n * n, 
                   ncol = n.mcmc, byrow = T)
  
  gam = matrix(NA, nrow = Time, ncol = n, 
               byrow = T)
  
  
  ## Initialize parameters
  
  X.param[,1] = sample(1:2, replace = T, Time)
  
  #lambda.param[, 1] = rgamma(n = n, shape = a, rate = b)
  
  lambda.param[1, 1] = rgamma(n = 1, shape = a, rate = b) #lambda low
  lambda.param[2, 1] = rgamma(n = 1, shape = c, rate = d) #change in lambda
  
  lambda.high = lambda.param[1,1] + lambda.param[2, 1] #lambda high, not needed, just a reminder  

  P.matrix = matrix(data = theta/(sum(theta[1,])), nrow = n, ncol = n, byrow = T) 
  
  P.param[,1] = as.vector(t(P.matrix))
  #holds all P.parameter values over runs
  
  ## Gibbs Updates
  
  for(l in 2:n.mcmc) {
    
    # print out every 10 iterations completed
    if( l %% 100 == 0 ) cat(paste("iteration", l, "complete\n")) 
    
    
    m = matrix(data = rep(0, n * n), nrow = n, ncol = n) 
    # number states going from i to j, refreshes every run
    
    P.matrix = matrix(data = c(P.param[, l - 1]), nrow = n, ncol = n, byrow = T)
    
    lambda.low = lambda.param[1, l - 1]
    lambda.high = lambda.low + lambda.param[2, l - 1]
    
    ##X Parameters
  
      gam[1, 1] = lambda.low ^ data[1] * exp(-lambda.low) * 
        delta[1] * P.matrix[1, X.param[2, l - 1]]
      
      gam[1, 2] = lambda.high ^ data[1] * exp(-lambda.high) * 
        delta[1] * P.matrix[1, X.param[2, l - 1]]
      
    
    
    X.param[1, l] = sample(x = (1:n), size = 1, prob = gam[1,])
    
    m[X.param[1, l], X.param[1, l]] = m[X.param[1, l], X.param[1, l]] + 1
    
    for(t in 2:(Time - 1)){
      
      # for(k in 1:n){
      #   gam[t, k] = (lambda.param[k, l - 1] ^ data[t]) * 
      #     exp(-lambda.param[k, l - 1]) * P.matrix[X.param[t - 1, l - 1], k] *
      #     P.matrix[k, X.param[t + 1, l - 1]]
      # }
      # 
      
      gam[t, 1] = lambda.low ^ data[t] * exp(-lambda.low) * 
         P.matrix[X.param[t - 1, l - 1], 1] *
            P.matrix[1, X.param[t + 1, l - 1]]
      
      gam[t, 2] = lambda.high ^ data[t] * exp(-lambda.high) * 
        P.matrix[X.param[t - 1, l - 1], 2] *
        P.matrix[2, X.param[t + 1, l - 1]]
      
      
      
      X.param[t, l] = sample(x = (1:n), 1,  prob = gam[t, ]) 
      
      m[X.param[t - 1, l], X.param[t, l]] = m[X.param[t - 1, l], 
                                              X.param[t,l]] + 1
    }
    
    # for(k in 1:n){
    #   gam[Time, k] = lambda.param[k, l - 1] ^ data[Time] * 
    #     exp(-lambda.param[k, l - 1]) * P.matrix[X.param[(Time - 1), l - 1], k]
    # }
    
    gam[Time, 1] = lambda.low ^ data[Time] * exp(-lambda.low) * 
      P.matrix[X.param[Time - 1, l - 1], 1] 
    
    gam[Time, 2] = lambda.high ^ data[Time] * exp(-lambda.high) * 
      P.matrix[X.param[Time - 1, l - 1], 2] 
    
    
    X.param[Time, l] = sample(x = 1:n, 1,  prob = gam[Time, ])
    
    m[X.param[Time - 1, l], X.param[Time, l]] = m[X.param[Time - 1, l], 
                                                  X.param[Time, l]] + 1
    
    #Lambda and P parameters
    for(h in 1:n) {

      # lambda.param[h, l] = rgamma(n = 1, shape =
      #                               sum(data[which(X.param[, l] == h)]) + a,
      #                             rate = sum(m[h, ]) + b )





      P.matrix[h, ] = (rdirichlet(n = 1 , alpha = theta[h, ] + m[h, ]))
      
    }
    
    lambda.param[1, l] = rgamma(n = 1, shape =
                        sum(data[which(X.param[, l] == 1)]) + a,
                        rate = sum(m[1, ]) + b )
    
    lambda.param[2, l] = rgamma(n = 1, shape =
                                  sum(data[which(X.param[, l] == 2)]) + c,
                                rate = sum(m[2, ]) + d)
    
    
    P.param[, l] = as.vector(t(P.matrix))
  }
  
  
  ## Rescale Lambda parameters into per minute 
  ## segments (instead of delta.t time segments)
  
  lambda.scale = lambda.param / delta.t * 60
  
  
  
  ## Compile the Estimates
  
  ## X1:XT, Lambda, Pmatrix
  
  #homes
  X.est = matrix(data = rep(NA, Time), nrow = Time, ncol = 1)
  lambda.est = matrix(data = rep(NA, n), nrow = n, ncol = 1)
  P.est = matrix(data = rep(NA, n * n), nrow = n * n, ncol = 1)
  
  
  for(t in 1:Time ){
    X.est[t, 1] = mean(X.param[t, ])  
  }
  
  # for(i in 1:n ){
  #   lambda.est[i, 1] = mean(lambda.scale[i, ])
  # }  
  # 
  
  lambda.est[1, 1] = mean(lambda.scale[1, ])
  lambda.est[2, 1] = mean(lambda.scale[1, ] + lambda.scale[2, ])
  
  for(l in 1:(n * n) ){
    P.est[l, 1] = mean(P.param[l, ])
  }
  
  P.est.matrix = matrix(data = c(P.est[, 1]), nrow = n, ncol = n, byrow = T)
  
  #plot the estimation runs.
  
  col = c("#120d08", "#bc5356", "#538bbc", "#53bc84")
  
  #lambda
  par(mfrow = c(2,2),
      oma = c(0,0,2,0) + 1,
      mar = c(1,1,1,1) + 3)
  
  
  plot(0,0,xlab="MCMC Runs",
       ylab="Lambda (scaled per 60 seconds)",
       ylim=c(0,max(lambda.scale)), 
       xlim=c(0,n.mcmc), 
       type="n",
       cex.lab = 1)
 
  for(i in 1:n){
    lines(1:n.mcmc, lambda.scale[i, ], col = col[i])
  }
  
  
  #P
  plot(0,0,xlab="MCMC Runs", ylab="P", ylim=c(0, max(P.param)), xlim=c(0,n.mcmc), 
       type="n", cex.lab = 1)
  for(i in 1:(n*n)){
    lines(1:n.mcmc, P.param[i, ], col = col[i])
  }
  
  #Single X
  X = X.param[sample(1:Time, 1), ]
  plot(0, 0, xlab="MCMC Runs", ylab="Single X", ylim=c(0,max(X)), 
       xlim=c(0,n.mcmc), type="n", cex.lab = 1)
  lines(1:n.mcmc, X, col = col[2])
  
  #States over time
  plot(X.est,type = "l",lwd=3, cex.lab = 1, col = col[1])
  
  title(main = title, outer=T)
  
  
  #########################################################
  ##
  ## Fancy Plots with Background Colors
  ##
  #########################################################
  
  par(mfrow = c(1, 1))
  
  
  if(length(unique(location)) == 1){
    
  ##High Density - 4 Hours
  plot(start, 1:int.num, main="High", 
       xlab="Seconds", ylab = "Cumulative Interaction Count", 
       xlim = c(0, maxtime))
  states = X.est #from code above
  rr = rle(states[,1])
  rr$values = round(rr$values, digits = 0)
  embedded.chain = rr$values
  cs = c(0,cumsum(rr$lengths))*delta.t - delta.t
  cols=c('#bc535644','#538bbc44')
  for(j in 1:length(embedded.chain)){
    rect(cs[j],0,cs[j + 1],int.num, 
         col = cols[embedded.chain[j]], density = NA)
    
  }
  points(start, 1:int.num, main="Low", xlab="Seconds",
         ylab = "Cumulative Interaction Count", 
         xlim=c(0,maxtime))
  }else{
    #Low Density - 4 Hours
  
  plot(start, 1:int.num, main="Low", xlab="Seconds",
       ylab = "Cumulative Interaction Count", 
       xlim=c(0,maxtime))
  states = X.est
  rr=rle(states[,1])
  rr$values = round(rr$values, digits = 0)
  embedded.chain=rr$values
  cs=c(0,cumsum(rr$lengths))*delta.t - delta.t
  cols=c('#bc535644','#538bbc44')
  for(j in 1:length(embedded.chain)){
    rect(cs[j],0,cs[j+1],int.num, col=cols[embedded.chain[j]] , density=NA)
  }
  points(start, 1:int.num, main="Low", xlab="Seconds",
         ylab = "Cumulative Interaction Count", 
         xlim=c(0,maxtime))

  # #Low Density - Location 1 
  # 
  # plot(low4.1$start_time, 1:nrow(low4.1), main="Low, Loc 1", xlab="Seconds", xlim=c(0,max(high4$end_time)))
  # ##    plot(one.day,1:length(one.day),main=day,xlab="Minutes")
  # states = out.low1$X.est
  # rr=rle(states[,1])
  # rr$values = round(rr$values, digits = 0)
  # embedded.chain=rr$values
  # cs=c(0,cumsum(rr$lengths))*delta.t - delta.t
  # cols=c('#FF000022','#0000FF22')
  # for(j in 1:length(embedded.chain)){
  #   rect(cs[j],0,cs[j+1],nrow(low4.1), col=cols[embedded.chain[j]] , density=NA)
  # }
  # points(cov.low$time, rep(0, nrow(cov.low)), 
  #        pch=8, col="forestgreen")
  # 
  # #Low Density - Location 4
  # plot(low4.4$start_time, 1:nrow(low4.4), main="Low, Loc 4", xlab="Seconds", xlim=c(0,max(high4$end_time)))
  # ##    plot(one.day,1:length(one.day),main=day,xlab="Minutes")
  # states = out.low4$X.est
  # rr=rle(states[,1])
  # rr$values = round(rr$values, digits = 0)
  # embedded.chain=rr$values
  # cs=c(0,cumsum(rr$lengths))*delta.t - delta.t
  # cols=c('#FF000022','#0000FF22')
  # for(j in 1:length(embedded.chain)){
  #   rect(cs[j],0,cs[j+1],nrow(low4.4), col=cols[embedded.chain[j]] , density=NA)
  # }
  # points(cov.low$time, rep(0, nrow(cov.low)), 
  #        pch=8, col="forestgreen")
  # 
}
  

  
  
  list(X.est = X.est, lambda.est = lambda.est, P.est = P.est.matrix, P.run = P.param)
  
}
 