#' Summary and Visualizations for Ant Trophallaxis Data
#'
#' The purpose of this function is to easily summarize and visualize 
#' ant trophallaxis data along with entrance times of foragers into 
#' the nest. 
#' 
#'
#' @param data, entrance, hours, density
#' @return  (1) - forager - number of forager ants
#'          (2) - forager.ID - ids and number of entrances of each forager
#'          (3) - allants - total number of ants engaged in trophallaxis
#'          (4) - all.ID - ids and number of interactions for each ant
#'          (5a) - Plot - Interactions per Ant with red marked mean
#'          (5b) - Plot - If low density, interactios per ant by location
#'          (6) - Plot -  Number of trophallaxis events over time with entrance
#'            times marked in red.
#'          (6b) - Plot - If low density, trophallaxis events over time by location
#' @export
#' @examples sum = summary.troph.data(data = troph.high.4, 
#'                  entrance = inout.high.4, 
#'                  hours = 4, 
#'                  density = "high")

summary.troph.data = function(data, entrance, hours, density = "high"){
  
#Number of unique ants entering (foragers)
num.inout = unique(entrance$Ant_ID)
total.inout = length(num.inout) #to be returned at end

for.ids = table(entrance$Ant_ID) #to be returned at end


# only entrance times
cov = entrance[which(entrance$Action == "enter" | entrance$Action == "Enter"), ]


#### How many ants interacting in each?

num.all = unique(data$Ant_ID)
total.all = length(num.all) #to be returned

all.ids = table(data$Ant_ID)

par(mfrow = c(1, 1))
hist(table(data$Ant_ID), xlab = "Count",
     main = "Interactions per Ant", 
     breaks = 20)
abline(v=mean(table(data$Ant_ID)), lty = 3, col="red", lwd = 3)

if(density != "high"){
  par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
  
  num.low1 = unique(data$Ant_ID[which(data$Location == 1)])
  length(num.low1)
  
  inter.num1 = table(data$Ant_ID[which(data$Location == 1)])
  inter.num1 #table of interactions per ant by ID
  
  num.low4 = unique(data$Ant_ID[which(data$Location == 4)])
  length(num.low4) 
  
  inter.num4 = table(data$Ant_ID[which(data$Location == 4)])
  inter.num4
  
  xlim = max(inter.num1, inter.num4) #useful so both histograms are on same xlim scale
  
  hist(table(data$Ant_ID[which(data$Location == 1)]), xlab = "Count",
       main = "",
       #main = "Interactions per Ant:Low Density, Loc 1", 
       breaks = 20, col = "black",
       xlim = c(0, xlim))
  
  hist(table(data$Ant_ID[which(data$Location == 4)]), xlab = "Count",
       #main = "Interactions per Ant:Low Density, Loc 4", 
       main = "",
       breaks = 15, col = "blue", 
       xlim = c(0, xlim))
  
  mtext("Interaction per Ant, by Location", outer = TRUE, cex=1.5)}

#########################################################
##
## Get rid of duplicate entries
## All entries recorded twice with each ant in main/partner position
##
#########################################################

data.change = data[seq(1, nrow(data), by = 2), ]

#########################################################
##
## plot time of trophylaxis events like a counting process
##
#########################################################

#order data frame by start time so plot works better
data.change = data.change[order(data.change$start_time), ]

par(mfrow = c(1, 1))
plot(data.change$start_time,1:nrow(data.change),main="High Density Trophallaxis",
     xlab = "Start Time", 
     ylab = "Number of Interactions")
points(cov$time, rep(0, nrow(cov)), pch=8, col="red")


if(density != "high"){
  
    ## Separate Low Density by Location
    
    low.1 = data.change[which(data.change$Location == 1), ]
    
    low.4 = data.change[which(data.change$Location == 4), ]
    low.4 = low.4[order(low.4$start_time), ]
    
    
    plot(data.change$start_time, 1:nrow(data.change), main="Low Density Trophallaxis",
         xlab = "Start Time", 
         ylab = "Number of Interactions", 
         col=data.change$Location)
    legend(5000, 100, c("Loc1", "Loc4"), lty = c(1,1), col = c("black", "blue"))
    points(cov$time, rep(0, nrow(cov)), pch=8, col="red")
    
    
    par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
    plot(low.1$start_time, 1:nrow(low.1), main="Location 1",
         xlim = c(0, max(data.change$end_time)),
         xlab = "Start Time", 
         ylab = "Number of Interactions")
    points(cov$time, rep(0, nrow(cov)), pch=8, col="red")
    
    plot(low.4$start_time, 1:nrow(low.4), main="Location 4",
         xlim = c(0, max(data.change$end_time)),
         xlab = "Start Time", 
         ylab = "Number of Interactions",
         col = "blue")
    points(cov$time, rep(0, nrow(cov)), pch=8, col="red")
    
    mtext("Low Density Trophallaxis", outer = TRUE, cex = 1.5 )
    
    }

list(forager = total.inout, forager.ID = for.ids, 
       allants = total.all, all.ID = all.ids)

}