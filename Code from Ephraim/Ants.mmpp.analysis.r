##
## read in data
##

## colony.ID = ant colony studied (1 or 2)
## date = day of observation (each colony was watched for 20minutes each day)
## ant.ID_1 = ID number of first ant
## ant.Type_1 = categorical classification of first ant
## ant.ID_2 = ID number of second ant
## ant.Type_2 = categorical classification of second ant
## time = time in days when trophylaxis started
## x = x coordinate of cell where trophylaxis occurred
## y = y coordinate of cell where trophylaxis occurred
## quadrant = quadrant of nest (1,2,3,4) where trophylaxis occurred
load("antdata.Rdata")
head(antdata)
str(antdata)



#########################################################
##
## Select one day for one colony
##
#########################################################

col1 = antdata[which(antdata$colony.ID == "col1"),]
str(col1)
help(str)
## see all days colony was observed
unique(col1$date)

## select one day
day=60213
one.day=col1[which(col1$date==day),]
str(one.day)
head(one.day)




#########################################################
##
## plot time of trophylaxis events like a counting process
##
#########################################################

plot(one.day$time,1:nrow(one.day),main="day60213,colony 1")


diff(one.day$time)
mean=mean(diff(one.day$time))
var(diff(one.day$time))
mean^2

par(mfrow=c(3,3))
for(day in unique(col1$date)){
one.day=col1[which(col1$date==day),]
plot(one.day$time,1:nrow(one.day),main=day)

}
    


day1=col1[which(col1$date==60113),]$time-60113
day2=col1[which(col1$date==60213),]$time-60213
day3=col1[which(col1$date==60313),]$time-60313
day4=col1[which(col1$date==60413),]$time-60413
day5=col1[which(col1$date==60513),]$time-60513
day6=col1[which(col1$date==60613),]$time-60613
day7=col1[which(col1$date==60713),]$time-60713       
day8=col1[which(col1$date==60813),]$time-60813

daylist=list(day1,day2,day3,day4,day5,day6,day7,day8)
alldays=c(day1,day2,day3,day4,day5,day6,day7,day8)
min(alldays)
max(alldays)
max(alldays)*24*60

par(mfrow=c(3,3))
for(day in 1:length(daylist)){
    one.day=daylist[[day]]
    plot(one.day,1:length(one.day),main=day)
}


########################################################
##
## Particle Filter Code
##
########################################################

#########################################################
##
## Simulation Function
##
## X_t ~ Markov Chain ( X_t-1 , P )
## Y_t ~ Pois (lambda_{X_t})
#########################################################

sim.mmpp <- function(tmax,delta.t,start.state=1,P,lambda){
    T=floor(tmax/delta.t)+1
    x=rep(NA,T)
    y=rep(NA,T)
    x[1]=start.state
    y[1]=0
    for(t in 2:T){
        ## sample latent state
        x[t]=sample(1:2,1,prob=P[x[t-1],])
        ## sample observed events
        y[t]=rpois(1,lambda=lambda[x[t]]*delta.t)
    }
    list(y=y,x=x,N=cumsum(y),delta.t=delta.t,t=(0:(T-1))*delta.t)
}



##
## simulate from the model
##

P = matrix(c(.99, .01, .01, .99), nrow = 2, byrow = T)
lambda = k = c(1, 4)
#theta=c(2,2)

delta.t = 1 #needs to be 1, else observations dependent on time

sim = sim.mmpp(7200, delta.t, start.state = 1, P, lambda)
par(mfrow = c(1, 1))
plot(sim$t, sim$N, type = "p", pch = ".", cex = 2, col = sim$x)
max(sim$y)



#########################################################
##
## Particle Filter 
##
## Goals: 1. estimate likelihood [y|P,k,theta]=E_x[y|x,P,k,theta]
##        2. estimate latent states {X_1,X_2,...X_T}
#########################################################

pfilter.mmpp <- function(y,P,lambda,delta.t,M,start.states=rep(1,M),print.iter=F){
    T=length(y)
    ## matrix to save states
    Xmat=matrix(NA,T,M)
    ## vector to save likelihood values from each time point
    lik.vals=rep(NA,T)
    ## initial state
    t=1
    current.states=start.states
    ##
    ## calculate likelihood of observation given states
    ## w_i=[y_1 | x_i1 , P,lambda]
    ##
    w=dpois(rep(y[t],M),lambda[current.states]*delta.t)
    ##
    ## save estimate of likelihood
    ## E_x[y_1 | x_1 , P,k,theta] = 1/M * sum_i [y_1 | x_i1 , P,k,theta]
    ##
    lik.vals[t]=mean(w)
    ##
    ## Bootstrap resample states with weights w_i/sum_i (w_i)
    ## 
    idx=sample(1:M,replace=TRUE,prob=w)
    Xmat[t,]=current.states[idx]
    ##
    ## loop through time
    ##
    for(t in 2:T){
        if(print.iter) cat(t," ")
        ##
        ## forward sample states
        ## X_it ~ Markov Chain ( X_it-1 , P )
        ## 
        prob.of.state.2=P[cbind(Xmat[t-1,],2)]
        current.states=1+rbinom(M,1,prob=prob.of.state.2)
        ## calculate likelihood of observation given states
        ## w_it=[y_t | x_it , P,k,theta]
        ##
        w=dpois(rep(y[t],M),lambda[current.states]*delta.t)
        ## save estimate of likelihood
        ## E_x[y_t | x_t , P,k,theta] = 1/M * sum_i [y_t | x_it , P,k,theta]
        lik.vals[t]=mean(w)
        ## Bootstrap resample states
        idx=sample(1:M,replace=TRUE,prob=w)
        Xmat[t,]=current.states[idx]
    }
    if(print.iter) cat("\n")
    loglik=sum(log(lik.vals))
    list(Xmat=Xmat,lik.vals=lik.vals,loglik=loglik)
}



## simulate from the model

P = matrix(c(.99, .01, .03, .97), nrow=2, byrow=T)
lambda = k = c(1, 4)
theta=c(2,2)

delta.t=.1

sim=sim.mmpp(200,delta.t,start.state=1,P,lambda)
plot(sim$t,sim$N,type="p",pch=".",cex=2,col=sim$x)
max(sim$y)

##
## particle filter with true parameters
##
M=100
start.states=1+rbinom(M,1,.5)
pf.true=pfilter.mmpp(sim$y,P,lambda,delta.t,M,start.states)

## loglikelihood
pf.true$loglik


##
## examine state recovery
##


true.states=sim$x

estimated.states=apply(pf.true$Xmat,1,mean)
plot(true.states,type="l",lwd=3)
points(estimated.states,type="l",col="red",lwd=2,lty=2)






##################################################################3
##
## Particle MCMC
##
##################################################################
#install.packages("mvtnorm")
library(mvtnorm)
source("mcmc.mmpp.r")

P.start=matrix(c(.75,.25,.25,.75),2)
P.start
lambda.start=c(1,3)

source("mcmc.mmpp.r")
out=mcmc.mmpp(sim$y,P.start,lambda.start,delta.t,M=100,diag(.005,4),n.mcmc=1000,adapt.max=500,adapt.int=100)
out$accept


params.true=c(P[1,1],P[2,2],lambda)

par(mfrow=c(2,2))
for(i in 1:4){
    plot(out$params[,i],type="l",main=colnames(out$params)[i])
    abline(h=params.true[i],col="red")
}



## examine state recovery under MCMC estimation

post.means=apply(out$params,2,mean)

P.est=matrix(c(post.means[1],1-post.means[1],1-post.means[2],post.means[2]),nrow=2,byrow=T)
P.est
lambda.est=post.means[3:4]
M=1000
start.states=1+rbinom(M,1,.5)
pf=pfilter.mmpp(sim$y,P.est,lambda.est,delta.t,M,start.states,print.iter=T)

true.states=sim$x

estimated.states=apply(pf$Xmat,1,mean)
plot(true.states,type="l",lwd=3)
points(estimated.states,type="l",col="red",lwd=2,lty=2)







#############################################3
##
##
## Apply to Ants
##
##
##############################################


delta.t=.1

daylist

ylist=list()
for(i in 1:length(daylist)){
    tmp=daylist[[i]]*24*60
    y=rep(0,20/delta.t)
    mint=0
    for(t in 1:length(y)){
        y[t]=length(which(tmp>mint & tmp<=mint+delta.t))
        mint=mint+delta.t
    }
    ylist[[i]]=y
}
ylist

n=integer()
for(i in 1:length(daylist)){
    n[i]=length(daylist[[i]])
}
n

yn=integer()
for(i in 1:length(daylist)){
    yn[i]=sum(ylist[[i]])
}
yn


P.start=matrix(c(.975,.025,.025,.975),2)
P.start
lambda.start=c(10,3)

source("mcmc.mmpp.ants.r")
out=mcmc.mmpp.ants(ylist,P.start,lambda.start,delta.t,M=100,diag(.005,4),n.mcmc=50,adapt.max=25000,adapt.int=100)
out$accept

out$


par(mfrow=c(2,2))
for(i in 1:4){
    plot(outsave$params[,i],type="l",main=colnames(out$params)[i])
}


outsave=out

mmpp.save=outsave
save(mmpp.save,file="mmpp.save.Rdata")

load("mmpp.save.Rdata")
outsave=mmpp.save

params.hat=apply(outsave$params,2,mean)
params.hat

## Values 8/5/2015
## > params.hat
##       p11       p22   lambda1   lambda2 
## 0.9860881 0.9900410 1.1226051 3.0700684 


lambda.hat=params.hat[3:4]
lambda.hat

P.hat=matrix(c(params.hat[1],1-params.hat[1],1-params.hat[2],params.hat[2]),2)
P.hat

## estimated rate matrix (in 1/minutes)
R.hat=P.hat/delta.t
diag(R.hat) <- 0

R.hat

##########################################
##
## State estimation
##
##########################################

M=10000

state.est=matrix(NA,8,length(ylist[[1]]))
for(i in 1:8){
    pf=pfilter.mmpp(ylist[[i]],P.hat,lambda.hat,delta.t,M=M,start.states=rep(1:2,each=M/2))
    state.est[i,]=apply(pf$Xmat,1,mean)
}

hist(state.est)
states=round(state.est)



par(mfrow=c(2,4))
for(day in 1:length(daylist)){
    one.day=daylist[[day]]*24*60
    plot(one.day,1:length(one.day),main=day,xlab="Minutes",ylim=c(0,80))
##    plot(one.day,1:length(one.day),main=day,xlab="Minutes")
    rr=rle(states[day,])
    embedded.chain=rr$values
    cs=c(0,cumsum(rr$lengths))*delta.t-delta.t
    cols=c('#FF000022','#0000FF22')
    for(j in 1:length(embedded.chain)){
        rect(cs[j],0,cs[j+1],100, col=cols[embedded.chain[j]] , density=NA)
    }
    points((0:199)/10,(state.est[day,]-1)*80,type="l",lty=1,lwd=7,col="#00FF0075")
##    axis(4,pretty(c(0,100)),col="green")
}



savePlot("mmpp.jpg")




#######################################################
#######################################################
#######################################################

P.start
lambda.start=c(1,3)
gamma.start=.10

source("mcmc.mmpp.ants.lasso.r")
out=mcmc.mmpp.ants.lasso(ylist,P.start,lambda.start,gamma.start,delta.t,M=100,Sigma.tune=diag(.0001,dim(P.start)[1]^2+1),n.mcmc=20000,adapt.max=25000,adapt.int=1000)
out$accept



par(mfrow=c(2,3))
for(i in 1:5){
    plot(out$params[,i],type="l",main=colnames(out$params)[i])
}



