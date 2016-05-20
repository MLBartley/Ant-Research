install.packages("devtools")
library("devtools")
install.packages("roxygen2")
library(roxygen2)


setwd("~/Google Drive/PSU/Projects/Ant-Research")
#create("Ants")

setwd("./Ants")
document()

install_github("MLBartley/Ant-Research", subdir = "Ants")
