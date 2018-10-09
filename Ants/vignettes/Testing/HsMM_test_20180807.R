## trying out HsMM via HSMM package

library(dplyr)
library(ggplot2)
library(hsmm)
library(magrittr)

source("./R/data_prep_functions.R")
source("./R/summary_visual_functions.R")


#data directory
dat.dir <- "../data-raw/"
#output directory
out.dir <- "../output/"


#column names - helps when only pulling in those columns, no extra

col_names <- c("Location", "Ant_ID", "Ant_ID_partner", "start_time", "end_time")

col2_low4 <- read.csv("./data-raw/Colony2_trophallaxis_low_density_4hr.csv")

#removes any extra columns, rows, and adds column names - depends on col_names being correct length
col2_low4 <- col2_low4[, 1:length(col_names)]
col2_low4 <- col2_low4 %>%
  tidyr::drop_na()
colnames(col2_low4) <- col_names


col2_low4_5 <- prep_troph_data(col2_low4, hours = 4, delta_t =  5)
qstarts <- col2_low4_5$queen_starts_persec
save(col2_low4_5, file = "./data-new/Col2Qstart.Rdata")

pi.par <- rep(.5, 2)
od.par <- list(lambda = c(.007, .05))
tpm.par <- matrix(c(0, 1, 1, 0), 2, 2)
rd.par <- list(r = c(1, 1), pi = c(.0001, .0001))

start.time <- Sys.time()
fit_hsmm_path <- hsmm.viterbi(x = col2_low4_5$queen_starts_persec, od = "pois",
                 rd = "nbinom", od.par = od.par,
                 pi.par = pi.par, rd.par = rd.par, tpm.par = tpm.par)
end.time <- Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='auto'), dig = 2)


sumvis_low <- sumvis_troph(data = col2_low4, entrance = F, hours = 4, density = "low")
#Note: cumul.lowqueen is ggplot saved to environment


cumul.lowqueen + geom_point(color = fit_hsmm_path$path[cumul.lowqueen$data$start_time])



start.time <- Sys.time()
fit_hsmm <- hsmm(x = col2_low4_5$queen_starts_persec, od = "pois",
                rd = "nbinom", od.par = od.par,
                pi.par = pi.par, rd.par = rd.par, tpm.par = tpm.par,
                detailed = TRUE, Q.max = 5000)
end.time <- Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='auto'), dig = 2)

save(fit_hsmm, file = "./data-new/fit_hsmm.Rdata")

## want to check convergence

lambdas <- matrix(as.numeric(unlist(fit_hsmm$ctrl$details$od)), nrow = 2, ncol = 1449)
n.bion.params <- matrix(as.numeric(unlist(fit_hsmm$ctrl$details$rd)), nrow = 4 , ncol = 1449)

plot(lambdas[1, ], type = "l")
plot(lambdas[2, ], type = "l")

plot(n.bion.params[1, ], type = "l")
plot(n.bion.params[2, ], type = "l")

plot(n.bion.params[3, ], type = "l")
plot(n.bion.params[4, ], type = "l")

# start.time <- Sys.time()
# fit_hsmm_smooth <- hsmm.smooth(x = col2_low4_5$queen_starts_persec, od = "pois",
#                  rd = "nbinom", od.par = od.par,
#                  pi.par = pi.par, rd.par = rd.par, tpm.par = tpm.par)
# end.time <- Sys.time()
# elapsed.time = round(difftime(end.time, start.time, units='auto'), dig = 2)
#
# cumul.lowqueen + geom_point(color = fit_hsmm_smooth$path[cumul.lowqueen$data$start_time])
# #same as without _smooth


##checking predictive ability of HsMM

load("./data-new/fit_hsmm.Rdata")

predict <- vector()
current <- 1
nbin.param.ests <- fit_hsmm$para$rd

for (t in 1:(length(fit_hsmm_path$path) - 1)) {

  r <- nbin.param.ests$r[fit_hsmm_path$path[t]]
  p <- nbin.param.ests$pi[fit_hsmm_path$path[t]]


  predict[t] <- dnbinom(x = current, size = r, prob = p) / sum(dnbinom(x = 0:current, size = r, prob = p))

  #how long has path been in current state?
  if (fit_hsmm_path$path[t] == fit_hsmm_path$path[t + 1]) { #stay

    predict[t] <- dnbinom(x = current + 1, size = r, prob = p) /
      (1 - sum(dnbinom(x = 0:current, size = r, prob = p)))

    current <- current + 1
  }else{ #switch
    predict[t] <- dnbinom(x = current, size = r, prob = p) /
      (1 - sum(dnbinom(x = 0:current, size = r, prob = p)))

    current <- 1}
}

head(predict)
plot(predict, col = fit_hsmm_path$path)

lambda <- fit_hsmm$para$od$lambda
x_t <- fit_hsmm_path$path

M <- as.matrix(cbind(predict, x_t))

osa_hsmm <- (apply(M, 1, function(x) x[1] * lambda[x[2]] + (1 - x[1]) * lambda[-x[2]] ))

Y <- as.matrix(cbind(osa_hsmm, col2_low4_5$queen_starts_persec))

MSPE_hsmm <- sum(apply(Y, 1, function(x) (x[1] - x[2])^2)) / nrow(Y) #.0200 which is good?! smaller is better






## function from zuchinni/langrock

hsmm2hmm<-function(omega,dm,eps=1e-10){
  mv<-sapply(dm,length)
  m<-length(mv)
  G<-matrix(0,0,sum(mv))
  for (i in 1:m){
    mi<-mv[[i]]
    F<-cumsum(c(0,dm[[i]][-mi]))
    ci<-ifelse(abs(1-F)>eps,dm[[i]]/(1-F),1)
    cim<-ifelse(1-ci>0,1-ci,0)
    Gi<-matrix(0,mi,0)
    for (j in 1:m){
      if(i==j) { if(mi==1)
      { Gi<-cbind(Gi,c(rep(0,mv[[j]]-1),cim))} else{
        Gi<-cbind(Gi,rbind(cbind(rep(0,mi-1),diag(cim[-mi],mi-1,mi-1)),
                  c(rep(0,mi-1),cim[[mi]])))}
      } else   { if(mi==1)
      { Gi<-cbind(Gi,matrix(c(omega[[i,j]]*ci,rep(0,mv[[j]]-1)),1))} else
      { Gi<-cbind(Gi,cbind(omega[[i,j]]*ci,matrix(0,mv[[i]],mv[[j]]-1)))}
      }
    }
    G<-rbind(G,Gi)
  }
  G }



