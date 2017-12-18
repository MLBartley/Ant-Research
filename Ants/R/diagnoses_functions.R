
#' Diagnoses of Penalized HMM Model 
#'
#' @param mcmc_matrix 
#' @param ant_file
#' @param chamber
#' @param Time
#' @param fig_path
#' @param fig_name
#' @param penalty
#'
#' 
#'
#' @return
#' @export
#'


penalty_diagnosis <- function(mcmc_matrix, ant_file, chamber, Time, fig_path, fig_name, penalty){
  
  
  source("http://www.stat.psu.edu/~mharan/batchmeans.R")
  
  #batch means (mean, s.error)
  bmeans <- bmmat(t(mcmc_matrix))
  
pdf( file = paste(fig_path, fig_name, round(penalty, 11), ".diagnostics", ".pdf", sep = ""))  
  #estimates vs sample size
  
    #lamdas - trop.rate(L/H)
    
    par(mfrow = c(1, 2))  
    
    estvssamp(samp = mcmc_matrix[1, ], plotname = "Trophallaxis Rate: Low")
    estvssamp(samp = (mcmc_matrix[1, ] + mcmc_matrix[2, ]), plotname = "Trophallaxis Rate: High")
    
    #gammas - st_switch_rate(LH/HL)
    
    estvssamp(samp = mcmc_matrix[4, ], plotname = "State Switching Rate: Low to High")
    estvssamp(samp = mcmc_matrix[5, ], plotname = "State Switching Rate: High to Low")
    
    #Probabilities
    
    par(mfrow = c(2, 2))
    
    estvssamp(samp = mcmc_matrix[6, ], plotname = "Probability Transitions: Low to Low")
    estvssamp(samp = mcmc_matrix[7, ], plotname = "State Switching Rate: Low to High")
    estvssamp(samp = mcmc_matrix[8, ], plotname = "State Switching Rate: High to Low")
    estvssamp(samp = mcmc_matrix[9, ], plotname = "State Switching Rate: High to High")
    
    
    #Random time for state chain - 

    par(mfrow = c(1, 1))
    
    sample_x <- sample(x = 10:nrow(mcmc_matrix), 1)
    # estvssamp(samp = mcmc_matrix[sample_x, ], plotname = "Behavior State: Random Time Chain")
    # 
  
  # visualization
    
  #lambdas - troph.rate  
  plot(coda::mcmc(t(mcmc_matrix[1:2, ])))
  
  #gammas - st_switch_rate
  plot(coda::mcmc(t(mcmc_matrix[4:5, ])))
  
  #ptm 
  plot(coda::mcmc(t(mcmc_matrix[6:9, ])))
  
  #rounded estimates over time for sto process
  
   plot(round(bmeans[-(1:9), 1 ]), type = "l")
    

   #fancy plot
   
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
   
   
   par(mfrow = c(1, 1))
   
  
     ## High Density - 4 Hours
     plot(start, 1:int.num, xlab = "Seconds", ylab = "Cumulative Interaction Count", 
          xlim = c(0, Time))
     states <- bmeans[-(1:9), 1]  #from code above
     rr <- rle(states)
     rr$values <- round(rr$values, digits = 0)
     embedded.chain <- rr$values
     cs <- c(0, cumsum(rr$lengths)) * delta_t - delta_t
     cols <- c("#bc535644", "#538bbc44")
     for (j in 1:length(embedded.chain)) {
       rect(cs[j], 0, cs[j + 1], int.num, col = cols[embedded.chain[j]], 
            density = NA)
       
     }
     points(start, 1:int.num,
            xlim = c(0, maxtime))
  
    dev.off()
  
  return(bmeans)

}