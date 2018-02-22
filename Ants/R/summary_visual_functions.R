
#' Summary and Visualizations for Ant Trophallaxis Data
#' 
#' The purpose of this function is to easily summarize and 
#' visualize ant trophallaxis data. In addition, if available, 
#' entrance times of foragers into the nest will also be 
#' reflected in the plots.
#' 
#' 
#' @param data Ant trophallaxis data file already loaded and 
#'   prepared with appropriate function.
#' @param entrance Data file of (presumably) forager ants 
#'   entering the colony. Used as a covariate in more complex 
#'   models.
#' @param hours Hours of observations reflected in trophallaxis 
#'   and entrance data files.
#' @param density Options are: \{high, low\} density. Summaries
#'   and plots will also devide data by chamber location in the
#'   "low" density case.
#'   
#' @return  (1) - forager - number of forager ants (2) - 
#'   forager.ID - ids and number of entrances of each forager (3)
#'   - allants - total number of ants engaged in trophallaxis (4)
#'   - all.ID - ids and number of interactions for each ant (5a) 
#'   - Plot - Interactions per Ant with red marked mean (5b) - 
#'   Plot - If low density, interactios per ant by location (6) -
#'   Plot -  Number of trophallaxis events over time with 
#'   entrance times marked in red. (6b) - Plot - If low density, 
#'   trophallaxis events over time by location
#' @export
#' @examples sum = sumvis.troph.data(data = troph.high.4, 
#'                  entrance = inout.high.4, 
#'                  hours = 4, 
#'                  density = "high")

sumvis_troph <- function(data, entrance = FALSE, hours, density = "high"){
  
  
  if (entrance != F) {
    
    #Number of unique ants entering (foragers)
    num.inout = unique(entrance$Ant_ID)
    total.inout = length(num.inout) #to be returned at end
    
    for.ids = table(entrance$Ant_ID) #to be returned at end
    
    # only entrance times
    cov = entrance[which(entrance$Action == "enter" | entrance$Action == "Enter"), ]
    
  }
  
  
  #### How many ants interacting in each?
  
  num.all = unique(data$Ant_ID)
  total.all = length(num.all) #to be returned
  
  all.ids = table(data$Ant_ID)
  
  #histogram of interactions per ant
  par(mfrow = c(1, 1))
  hist.combined <<- qplot(as.numeric(table(data$Ant_ID)), 
                          geom = "histogram",
                          main = "Interactions per Ant", 
                          bins = 20,
                          xlab = "Number of Interactions",
                          xlim = c(0, max(50, max(table(data$Ant_ID))))) +
    geom_vline(xintercept = mean(as.numeric(table(data$Ant_ID))),
               color = "#bc5356", size = 3)
  
  
  print(hist.combined)
  
  # hist(table(data$Ant_ID), xlab = "Count",
  #   main = "Interactions per Ant", 
  #   breaks = 20, 
  #   xlim = c(0, max(50, max(table(data$Ant_ID)))) )
  # abline(v = mean(table(data$Ant_ID)), lty = 3, col="#bc5356", lwd = 3)
  
  #need to split by location 
  if(density != "high"){
    # par(mfrow=c(1,2), oma = c(0, 0, 2, 0)) - want individual ggplots
    
    num.low1 = unique(data$Ant_ID[which(data$Location == 1)])
    length(num.low1)
    
    inter.num1 = table(data$Ant_ID[which(data$Location == 1)])
    inter.num1 #table of interactions per ant by ID
    
    num.low4 = unique(data$Ant_ID[which(data$Location == 4)])
    length(num.low4) 
    
    inter.num4 = table(data$Ant_ID[which(data$Location == 4)])
    inter.num4
    
    xlim = max(inter.num1, inter.num4) #useful so both histograms are on same xlim scale
    
    hist.lowqueen <<- qplot(as.numeric(table(data$Ant_ID[which(data$Location == 1)])),
                            xlab = "Number of Interactions",
                            # main = "",
                            main = "Interactions per Ant - Queen's Chamber",
                            bins = 20, 
                            # fill = I("#120d08"),
                            xlim = c(0, xlim)) +
      geom_vline(xintercept = mean(as.numeric(table(data$Ant_ID[which(data$Location == 1)]))),
                 color = "#bc5356", size = 3)
    
    hist.lowenter <<- qplot(as.numeric(table(data$Ant_ID[which(data$Location == 4)])),
                            xlab = "Number of Interactions",
                            main = "Interactions per Ant - Entrance Chamber",
                            # main = "",
                            bins = 15, 
                            # fill = "#538bbc", 
                            xlim = c(0, xlim)) +
      geom_vline(xintercept = mean(as.numeric(table(data$Ant_ID[which(data$Location == 4)]))),
                 color = "#bc5356", size = 3)
    
    # mtext("Interaction per Ant, by Location", outer = TRUE, cex=1.5)
    
    print(hist.lowenter)
    print(hist.lowqueen)
  }
  
  # #########################################################
  # ##
  # ## Get rid of duplicate entries
  # ## All entries recorded twice with each ant in main/partner position
  # ##
  # #########################################################
  # 
  data.change = data[seq(1, nrow(data), by = 2), ]
  # 
  # #########################################################
  # ##
  # ## plot time of trophylaxis events like a counting process
  # ##
  # #########################################################
  # 
  # #order data frame by start time so plot works better
  data.change = data.change[order(data.change$start_time), ]
  # 
  par(mfrow = c(1, 1))
  
  if(density == "high"){
    
    cumul.h <<- ggplot(data = data.change, 
                       aes(data.change$start_time, 1:nrow(data.change))) +
      ggtitle("Single Chamber Trophallaxis") +
      xlab("Start Time (seconds)") +
      ylab("Cumulative Number of Interactions") +
      geom_point(shape = 2, color = "#120d08", fill = "white")
    
    # plot(data.change$start_time,
    #   1:nrow(data.change),
    #   main = "High Density Trophallaxis",
    #   xlab = "Start Time", 
    #   ylab = "Number of Interactions", 
    #   col = "#120d08")
    if (entrance != F) {
      cumul.h <<- cumul.h + 
        geom_point(cov$time, rep(0, nrow(cov)), shape = 8, col =	"#53bc84")
    }
    
    print(cumul.h)
  }
  else{
    
    # ## Separate Low Density by Location
    # 
    low.1 = data.change[which(data.change$Location == 1), ]
    # 
    low.4 = data.change[which(data.change$Location == 4), ]
    low.4 = low.4[order(low.4$start_time), ]
    # 
    
    cumul.l <<- ggplot(data = data.change, 
                       aes(data.change$start_time, 1:nrow(data.change))) +
      ggtitle("Combined Chambers Trophallaxis") +
      xlab("Start Time (seconds)") +
      ylab("Cumulative Number of Interactions") +
      geom_point(shape = 2, aes(color = factor(Location)), fill = "white") +
      labs(color = "Location")
    
    
    cumul.lowqueen <<- ggplot(data = low.1, 
                              aes(low.1$start_time, 1:nrow(low.1))) +
      ggtitle("Queen's Chamber Trophallaxis") +
      xlab("Start Time (seconds)") +
      ylab("Cumulative Number of Interactions") +
      geom_point(shape = 2, color = "#120d08", fill = "white")
    
    cumul.lowenter <<- ggplot(data = low.4, 
                              aes(low.4$start_time, 1:nrow(low.4))) +
      ggtitle("Entrance Chamber Trophallaxis") +
      xlab("Start Time (seconds)") +
      ylab("Cumulative Number of Interactions") +
      geom_point(shape = 2, color = "#538bbc", fill = "white")
    
    # plot(data.change$start_time, 1:nrow(data.change), main = "Low Density Trophallaxis",
    #   xlab = "Start Time", 
    #   ylab = "Number of Interactions", 
    #   col = c("#120d08", "#538bbc"), 
    #   pch = c(1, 2))
    # legend(50, 300, c("Loc1", "Loc4"), 
    #   pch = c(1,2), 
    #   col = c("#120d08", "#538bbc"),
    #   text.width = 1200)
    if (entrance != F) {
      cumul.l <<- cumul.l +
        geom_point(cov$time, rep(0, nrow(cov)), shape = 8, col = "#53bc84")
    }
    print(cumul.l)
    print(cumul.lowenter)
    print(cumul.lowqueen)
    # par(mfrow = c(1, 1))
    # # par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
    # plot(low.1$start_time, 1:nrow(low.1), main = "Queen's Chamber",
    #   xlim = c(0, max(data.change$end_time)),
    #   xlab = "Interaction Time (s)", 
    #   ylab = "Number of Interactions",
    #   col = "#120d08")
    # if (entrance != F) {
    #       points(cov$time, rep(0, nrow(cov)), pch = 8, col = "#53bc84")
    # }
    # 
    # plot(low.4$start_time, 1:nrow(low.4), main = "Entrance Chamber",
    #   xlim = c(0, max(data.change$end_time)),
    #   xlab = "Interaction Time (s)", 
    #   ylab = "Number of Interactions",
    #   col = "#538bbc", 
    #   pch = 2)
    # if (entrance != F) {
    #       points(cov$time, rep(0, nrow(cov)), pch = 8, col = "#53bc84")
    # }
    # 
    # mtext("Low Density Trophallaxis", outer = TRUE, cex = 1.5 )
    
  }
  
  return(list(allants = total.all, all.ID = all.ids))
  
  if (entrance != F) {
    return(forager = total.inout, forager.ID = for.ids)
  }
  # print(h)
}





#' Post MCMC Summary Table and Visuals
#'
#' @param results
#' @param compare
#' @param file_path
#' @param file_name
#'
#' @export
#'
#'


sumtable_model <- function(results, compare, file_path, file_name, model){
  
  
  #Simple Model:  start rates (high, low), and PTM (P)
  #Penalized Model:  start rates (high, low), and switch rate (gamma)
  #Penalized/Covariate Model: start rates (high, low), switch rate parameters (alpha, betas)
  
  #note that states over time (X_t) are summerized already with figure
  
  
  
  #Create Homes for vectors
  
  st_rate_low_est = rep(0, length(compare))
  st_rate_high_est = rep(0, length(compare))

  tpm_11_est = rep(0, length(compare))
  tpm_12_est = rep(0, length(compare))
  tpm_21_est = rep(0, length(compare))
  tpm_22_est = rep(0, length(compare))

  if (model == "penalized") {
    sw_rate_low_est =  rep(0, length(compare))
    sw_rate_high_est = rep(0, length(compare))
  }

  if (model == "covariate"){

  }

  
  MSPE_est = rep(0, length(compare))
  
  accept = rep(0, length(compare))
  
  # 
  # 
  # for(i in 1:length(compare)){
  #   st_rate_low_est[i] =  results[[i]]$st_rates_est[[1]]$est
  # }
  # 
  # 
  # for(i in 1:length(compare)){
  #   st_rate_high_est[i] = results[[i]]$st_rates_est[[2]]$est
  # }
  # 
  # 
  # if (model != "simple") {
    
  #   for(i in 1:length(compare)){
  #     sw_rate_low_est[i] = results[[i]]$sw_rates_est[[1]]$est
  #   }
  #   
  #   for(i in 1:length(compare)){
  #     sw_rate_high_est[i] = results[[i]]$sw_rates_est[[2]]$est
  #   }
  # }
  # 
  # 
  # 
  # for(i in 1:length(compare)){
  #   tpm_11_est[i] =  results[[i]]$st_ptm_est[[1]]$est
  # }
  # 
  # 
  # for(i in 1:length(compare)){
  #   tpm_12_est[i] = results[[i]]$st_ptm_est[[2]]$est
  # }
  # 
  # 
  # for(i in 1:length(compare)){
  #   tpm_21_est[i] =  results[[i]]$st_ptm_est[[3]]$est
  # }
  # 
  # 
  # for(i in 1:length(compare)){
  #   tpm_22_est[i] = results[[i]]$st_ptm_est[[4]]$est
  # }
  # #
  # #
  # #
    for(i in 1:length(compare)){
      MSPE_est[i] =  results[[i]]$MSPE
    }
  #
  #
  if (model != "simple") {
    
    
    for(i in 1:length(compare)){
      accept[i] = results[[i]]$accept
    }
  }
  
  
  if (model == "simple"){
    table = data.frame(compare, st_rate_low_est, st_rate_high_est,
                       tpm_11_est, tpm_12_est, tpm_21_est, tpm_22_est,
                       MSPE_est
    )
  }
  else{
    table = data.frame(compare, st_rate_low_est, st_rate_high_est,
                       sw_rate_low_est, sw_rate_high_est,
                       tpm_11_est, tpm_12_est, tpm_21_est, tpm_22_est,
                       MSPE_est, accept)
  }
  #
  #
  #
  #
  #
  # plot(log(compare), MSPE_est, col = ifelse(accept<=1000,"red","black"),
  #      ylim = c(min(MSPE_est) - 2*sd(MSPE_est), max(MSPE_est) + 2 * sd(MSPE_est)))
  
  if (model == "simple"){
    compare_plot <<- ggplot(data = table, aes((compare), MSPE_est)) +
      geom_point() +
      ylim(min(MSPE_est) - 2*sd(MSPE_est), max(MSPE_est) + 2 * sd(MSPE_est))
  }
  
  else{
    compare_plot <<- ggplot(data = table, aes(log(compare), MSPE_est)) +
      geom_point(aes( x = log(compare), colour = accept <= 1000))
  }
  
  
  write.csv(x = table, file = paste(file_path, file_name, ".csv", sep = "") )
  #
}


#' Pst MCMC State Switching Visualization
#'
#' @param mcmc_matrix 
#' @param fig_path Path needed to sent plot figures to folder.
#' @param fig_name Base name of plot files to be saved.
#'
#' @export
#'
plot_states  <- function(mcmc_matrix, fig_path, fig_name){

  
  ###### Fancy Plots with Background Colors

    jpeg( file = paste(fig_path, fig_name, round(penalty, 11), ".states", ".jpg", sep = ""))
  
  
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

  
    dev.off()
  
}