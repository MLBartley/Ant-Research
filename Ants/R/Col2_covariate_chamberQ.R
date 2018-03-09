# library(magrittr)
# library(dplyr)
#
#
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
# for (antnum in 1:length(c2l.list)) {
#   ant = c2l.list[[antnum]]
#   # plot(nest.poly.ld, type = "l")
#   # points(ant[, -1], type = "l")
#   T = nrow(ant)
#   times = which(ant$x[-1] < 40 & ant$x[-T] > 40)
#   times
#   if (length(times) > 0) {
#     # points(ant[times:(times + 1), -1], col = "red")
#
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
# for (antnum in 1:length(c2l.list)) {
#   ant <- c2l.list[[antnum]]
#   T <- nrow(ant)
#   times <- which(ant$x[-1] > 198 & ant$x[-T] < 198)
#   times
#   if (length(times) > 0) {
#
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
# unique(ant.id.out)
