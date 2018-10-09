<<<<<<< HEAD
library(magrittr)
library(dplyr)

load("./data-raw/c2l.list.20170608.Rdata")
load("./data-raw/c1l.list.20170608.Rdata")
# ls()
str(c2l.list[[1]])
str(c1l.list[[1]])
load("./data-raw/nest.poly.ld.Rdata")
plot(nest.poly.ld, type = "l")

# #single ant's points - chamber 2
# points(c2l.list[[34]][, -1])
#
# names(c2l.list)
#
# #single ant's points - chamber 1
#
# points(c1l.list[[2]][, -1])
=======
# library(magrittr)
# library(dplyr)

# load("./data-raw/c2l.list.20170608.Rdata")
# ls()
# str(c2l.list[[1]])
# load("./data-raw/nest.poly.ld.Rdata")
# plot(nest.poly.ld, type = "l")
# 
# #single ant's points
# points(c2l.list[[34]][, -1])
# 
# names(c2l.list)
# 
# ##
# ## code to find when ants from colony 2 enter queen chamber
# ##
# 
# ant.id = integer()
# times.enter.chamber = integer()
# time.since.out = integer()
# 

#
# #single ant's points
# points(c2l.list[[34]][, -1])
#
# names(c2l.list)
>>>>>>> 8cc3341e4c9336c9a157eb5736a8c32b9c22651a
#
# ##
# ## code to find when ants from colony 2 enter queen chamber
# ##
#
# ant.id = integer()
# times.enter.chamber = integer()
# time.since.out = integer()
#
<<<<<<< HEAD
#
#
# #single ant's points
# points(c2l.list[[34]][, -1])
#
# names(c2l.list)

##
## code to find when ants from colony 2 enter queen chamber
##

ant.id.Q = integer()
ant.id.E = integer()
times.enter.Qchamber = integer()
times.enter.Echamber = integer()
time.since.out = integer()

for (antnum in 1:length(c2l.list)) {
  ant = c2l.list[[antnum]]
  # plot(nest.poly.ld, type = "l")
  # points(ant[, -1], type = "l")
  T = nrow(ant)
  times.Q = which(ant$x[-1] < 40 & ant$x[-T] > 40)
  times.Q
  times.E = which(ant$x[-1] < 199 & ant$x[-T] == 199)
  times.E
  times.E2 = which(ant$x[-1] < 160 & ant$x[-T] > 160)
  times.E2
  if (length(times.Q) > 0) {
    points(ant[times.Q:(times.Q + 1), -1], col = "red")

    for (i in 1:length(times.Q)) {
      ant.id.Q = c(ant.id.Q, names(c2l.list)[[antnum]])
      times.enter.Qchamber = c(times.enter.Qchamber, times.Q[i])
    }
  }

  if (length(times.E) > 0) {
    points(ant[times.E:(times.E + 1), -1], col = "blue")

    for (i in 1:length(times.E)) {
      ant.id.E = c(ant.id.E, names(c2l.list)[[antnum]])
      times.enter.Echamber = c(times.enter.Echamber, times.E[i])
    }
  }

  if (length(times.E2) > 0) {
    points(ant[times.E2:(times.E2 + 1), -1], col = "darkblue")

    for (i in 1:length(times.E2)) {
      ant.id.E = c(ant.id.E, names(c2l.list)[[antnum]])
      times.enter.Echamber = c(times.enter.Echamber, times.E2[i])
    }
  }
}

ant.id.Q
times.enter.Qchamber

ant.id.E
times.enter.Echamber

covariate <- data.frame(Ant_ID = ant.id, time = times.enter.chamber) %>%
            mutate(Action = rep("Exit"))



### CURRENTLY USED FOR COVARIATE - NEEDS TO BE UPDATED WITH INFO FROM LAUREN
save(covariate, file = "./data-raw/Colony2_covariate_low_density_4hr.Rda")
###

##
## code to find when ants from colony 1 enter queens chamber
##

ant.id = integer()
times.enter.chamber = integer()
time.since.out = integer()

for (antnum in 1:length(c1l.list)) {
  ant = c1l.list[[antnum]]
  # plot(nest.poly.ld, type = "l")
  # points(ant[, -1], type = "l")
  T = nrow(ant)
  times = which(ant$x[-1] < 40 & ant$x[-T] > 40)
  times
  if (length(times) > 0) {
    points(ant[times:(times + 1), -1], col = "red")

    for (i in 1:length(times)) {
      ant.id = c(ant.id, names(c1l.list)[[antnum]])
      times.enter.chamber = c(times.enter.chamber, times[i])
    }
  }
}

ant.id
times.enter.chamber

covariate <- data.frame(Ant_ID = ant.id, time = times.enter.chamber) %>%
  mutate(Action = rep("Exit"))



### CURRENTLY USED FOR COVARIATE - NEEDS TO BE UPDATED WITH INFO FROM LAUREN
save(covariate, file = "./data-raw/Colony1_covariate_low_density_4hr.Rda")
###


##
## code to fine which ants from colony 2 exited the nest
##

ant.id.out <- integer()
times.exit.nest <- integer()


ant.id
times.enter.chamber

covariate <- data.frame(Ant_ID = ant.id, time = times.enter.chamber) %>%
            mutate(Action = rep("Exit"))



### CURRENTLY USED FOR COVARIATE - NEEDS TO BE UPDATED WITH INFO FROM LAUREN
save(covariate, file = "./data-raw/Colony2_covariate_low_density_4hr.Rda")
###


##
## code to fine which ants from colony 2 exited the nest
##

ant.id.out <- integer()
times.exit.nest <- integer()

for (antnum in 1:length(c2l.list)) {
  ant <- c2l.list[[antnum]]
  T <- nrow(ant)
  times <- which(ant$x[-1] > 198 & ant$x[-T] < 198)
  times
  if (length(times) > 0) {
    for (i in 1:length(times)) {
      ant.id.out <- c(ant.id.out, names(c2l.list)[[antnum]])
      times.exit.nest <- c(times.exit.nest, times[i])
    }
  }
}


ant.id.out
times.exit.nest



ant.id.out
times.exit.nest

unique(ant.id.out)
=======
# for (antnum in 1:length(c2l.list)) {
#   ant = c2l.list[[antnum]]
#   # plot(nest.poly.ld, type = "l")
#   # points(ant[, -1], type = "l")
#   T = nrow(ant)
#   times = which(ant$x[-1] < 40 & ant$x[-T] > 40)
#   times
#   if (length(times) > 0) {
#     # points(ant[times:(times + 1), -1], col = "red")

#     for (i in 1:length(times)) {
#       ant.id = c(ant.id, names(c2l.list)[[antnum]])
#       times.enter.chamber = c(times.enter.chamber, times[i])
#     }
#   }
# }
# 
# ant.id
# times.enter.chamber
# 
# covariate <- data.frame(Ant_ID = ant.id, time = times.enter.chamber) %>%
#             mutate(Action = rep("Exit"))
# 
# 
# 
# ### CURRENTLY USED FOR COVARIATE - NEEDS TO BE UPDATED WITH INFO FROM LAUREN
# save(covariate, file = "./data-raw/Colony2_covariate_low_density_4hr.Rda")
# ###
# 
# 
# ##
# ## code to fine which ants from colony 2 exited the nest
# ##
# 
# ant.id.out <- integer()
# times.exit.nest <- integer()
# 
#
# ant.id
# times.enter.chamber
#
# covariate <- data.frame(Ant_ID = ant.id, time = times.enter.chamber) %>%
#             mutate(Action = rep("Exit"))
#
#
#
# ### CURRENTLY USED FOR COVARIATE - NEEDS TO BE UPDATED WITH INFO FROM LAUREN
# save(covariate, file = "./data-raw/Colony2_covariate_low_density_4hr.Rda")
# ###
#
#
# ##
# ## code to fine which ants from colony 2 exited the nest
# ##
#
# ant.id.out <- integer()
# times.exit.nest <- integer()
#
# for (antnum in 1:length(c2l.list)) {
#   ant <- c2l.list[[antnum]]
#   T <- nrow(ant)
#   times <- which(ant$x[-1] > 198 & ant$x[-T] < 198)
#   times
#   if (length(times) > 0) {
#     for (i in 1:length(times)) {
#       ant.id.out <- c(ant.id.out, names(c2l.list)[[antnum]])
#       times.exit.nest <- c(times.exit.nest, times[i])
#     }
#   }
# }
# 
# 
# ant.id.out
# times.exit.nest
# 
#
#
# ant.id.out
# times.exit.nest
#
# unique(ant.id.out)
>>>>>>> 8cc3341e4c9336c9a157eb5736a8c32b9c22651a
