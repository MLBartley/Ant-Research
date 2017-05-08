#
#
# Exporing Penalty vs MSPE
#
#
##################################


#Import in 1 second .csv files
#append data all together
# take log(penalty)
# scatter plot vs MSPE
# 

setwd("~/Google Drive/PSU/Projects/Ant-Research/simulation-practice/CSV-Results ")
filenames <- list.files(path = "./")
one.data <- do.call("rbind", lapply(filenames, read.csv, header = TRUE))

#remove 'truth' lines (penalty = 0)
one.data = one.data[-which(one.data$c.0..penalty. == 0), ]

one.data$penalty <- log(one.data$c.0..penalty.)

plot(one.data$penalty, one.data$MSPE.est, 
  col = ifelse(one.data$accept<=1000,"red","black"), 
  ylab = "MSPE",
  xlab = "Penalty (exp{value})",
  ylim = c(.055, .06))
abline(v = -13, col = "blue")
# lines(predict(lm(one.data$MSPE.est~one.data$penalty+I(one.data$penalty^2))))


setwd("~/Google Drive/PSU/Projects/Ant-Research/")


