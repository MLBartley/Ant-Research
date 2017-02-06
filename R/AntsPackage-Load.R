install.packages("devtools")
library("devtools")
install.packages("roxygen2")
library(roxygen2)

getwd()

setwd("~/Google Drive/PSU/Projects/Ant-Research")
#create("Ants")

setwd("./Ants")
document()

#go updata github first
install_github("MLBartley/Ant-Research", subdir = "Ants", force=T)
library(Ants)

setwd("~/Google Drive/PSU/Projects/Ant-Research")
