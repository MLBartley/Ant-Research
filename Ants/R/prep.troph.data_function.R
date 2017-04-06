#' Data Preparation for Ant Trophallaxis Data File
#'
#' The purpose of this function is to take in .csv files 
#' and apply necissary changes to the data format. Changes include 
#' removing redundant entries, sorting by entrance times, and bining data
#' into smaller chunks.
#' 
#' @param data, delta.t
#' @return (1) One or three (depending on high or low density, 
#'             respectively), file lists with amount of interaction
#'             within binned time chunk (delta.t).
#' @export
#' @examples
#'  prep.troph.data(high4, 60)
#' 
#' 


prep.troph.data = function(data, delta.t) {
    
    
    
    # remove duplicate entries
    data.change = data[seq(1, nrow(data), by = 2), ]
    
    # order data frame by start time so plot works better
    data.change = data.change[order(data.change$start_time), ]
    
    
    # combining the (start time) data into delta.t increments
    if (length(unique(data.change$Location)) == 1) {
        
        # High Density Data
        for (i in 1:nrow(data.change)) {
            tmp = data.change$start_time
            y = rep(0, max(data.change$end_time)/delta.t)
            mint = 0
            for (t in 1:length(y)) {
                y[t] = length(which(tmp > mint & tmp <= mint + delta.t))  #how many interactions start in delta.t increment
                mint = mint + delta.t
            }
        }
        list(high.y = y)
    } else {
        ## Separate Low Density by Location
        low.1 = data.change[which(data.change$Location == 1), ]
        
        low.4 = data.change[which(data.change$Location == 4), ]
        low.4 = low.4[order(low.4$start_time), ]
        
        # Low Density Data: Both Locations
        for (i in 1:nrow(data.change)) {
            tmp = data.change$start_time
            y = rep(0, max(data.change$end_time)/delta.t)  #max time same ~7200
            mint = 0
            for (t in 1:length(y)) {
                y[t] = length(which(tmp > mint & tmp <= mint + delta.t))
                mint = mint + delta.t
            }
        }
        low.y = y
        
        # Low Density Data: Location 1
        for (i in 1:nrow(data.change)) {
            tmp = low.1$start_time
            y = rep(0, max(data.change$end_time)/delta.t)
            mint = 0
            for (t in 1:length(y)) {
                y[t] = length(which(tmp > mint & tmp <= mint + delta.t))
                mint = mint + delta.t
            }
        }
        low1.y = y
        
        
        # Low Density Data: Location 4
        for (i in 1:nrow(data.change)) {
            tmp = low.4$start_time
            y = rep(0, max(data.change$end_time)/delta.t)
            mint = 0
            for (t in 1:length(y)) {
                y[t] = length(which(tmp > mint & tmp <= mint + delta.t))
                mint = mint + delta.t
            }
        }
        low4.y = y
        
        list(low.y = low.y, low1.y = low1.y, low4.y = low4.y)
    }
    
}
