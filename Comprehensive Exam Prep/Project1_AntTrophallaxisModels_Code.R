####
####
#### Project One Final* Code
#### February 2017
####


# *Wishful thinking; who am I kidding? 


###Outline
    # Ant data visualization - 4 hour High/Low density
            #histograms of events per ant
            #scatterplots of each density (and each location in low density)
            #scatter plots of each density with entrance times     
    # Ant data simple model - want to show motivation - doesn't work!
    # Simulated data - visualization 
            # scatter plots of each simulated density?   
    # Simulated data penalized model - show that it can be effective
    # Simulated data penalized model with covariates - include the biology
  

#Call in Trophallaxis Data
    # Currently the .cvs files also load a bunch of columns that are empty 
    # In low density, chamber 4 is entrance and chamber 1 has queen

col1.high4 <- read.csv("./Data/Colony1_trophallaxis_high_density_4hr.csv")
col1.low4 <- read.csv("./Data/Colony1_trophallaxis_low_density_4hr.csv")

col2.high4 <- read.csv("./Data/Colony2_trophallaxis_high_density_4hr.csv")
col2.low4 <- read.csv("./Data/Colony2_trophallaxis_low_density_4hr.csv")

col3.high4 <- read.csv("./Data/Colony3_trophallaxis_high_density_4hr.csv")
col3.low4 <- read.csv("./Data/Colony3_trophallaxis_low_density_4hr.csv")


#Call in Foraging Data
    #



#Prep Data - note decided to keep prep.torph.data function 
            #output to be just the binned data, otherwise function 
            #takes too long to do 1 sec every time

col1.high4.1 = prep.troph.data(col1.high4, 1)
col1.low4.1 = prep.troph.data(col1.low4, 1)

col1.high4.5 = prep.troph.data(col1.high4, 5)
col1.low4.5 = prep.troph.data(col1.low4, 5)

col1.high4.30 = prep.troph.data(col1.high4, 30)
col1.low4.30 = prep.troph.data(col1.low4, 30)

col2.high4.1 = prep.troph.data(col2.high4, 1)
col2.low4.1 = prep.troph.data(col2.low4, 1)

col2.high4.5 = prep.troph.data(col2.high4, 5)
col2.low4.5 = prep.troph.data(col2.low4, 5)

col2.high4.30 = prep.troph.data(col2.high4, 30)
col2.low4.30 = prep.troph.data(col2.low4, 30)

col3.high4.1 = prep.troph.data(col3.high4, 1)
col3.low4.1 = prep.troph.data(col3.low4, 1)

col3.high4.5 = prep.troph.data(col3.high4, 5)
col3.low4.5 = prep.troph.data(col3.low4, 5)

col3.high4.30 = prep.troph.data(col3.high4, 30)
col3.low4.30 = prep.troph.data(col3.low4, 30)



#Visualize the Trophallaxis Data
    #Want these to save to .pdf (in new folder?)
    

par(mfrow = c(1, 1))
plot(col1.high4.1$start_time,1:nrow(col1.high4.1),main="High Density Trophallaxis, 4 Hours",
     xlab = "Start Time", 
     ylab = "Number of Interactions")


plot(low4$start_time, 1:nrow(low4), main="Low Density Trophallaxis",
     xlab = "Start Time", 
     ylab = "Number of Interactions", 
     col=low4$Location)
legend(5000, 100, c("Loc1", "Loc4"), lty = c(1,1), col = c("black", "blue"))



