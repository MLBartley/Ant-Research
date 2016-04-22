mcmc.mmpp.ants <- function(y.list,P.start,lambda.start,delta.t,M,Sigma.tune,n.mcmc,adapt.max=10000,adapt.int=1000){
    ## preliminaries
    T=length(y[[1]])
    params=c(P.start[1,1],P.start[2,2],lambda.start)

    params.save=matrix(NA,nrow=n.mcmc,ncol=4)
    colnames(params.save) <- c("p11","p22","lambda1","lambda2")

    ##
    ## function to evaluate loglikelihood of data given params
    ##

    mmpp.loglik <- function(params,y,M,delta.t){
        p11=params[1]
        p22=params[2]
        lambda=params[3:4]
        ## create P matrix
        P=matrix(c(p11,1-p11,
            1-p22,p22),byrow=T,nrow=2)
        ## draw a random set of initial states
        start.states=sample(1:2,size=M,replace=T)
        ll=0
        for(i in 1:length(y)){
            pf=pfilter.mmpp(y[[i]],P,lambda,delta.t,M,start.states)
            ll=ll+pf$loglik
        }
        ll
    }

    ##
    ## function to evaluate loglikelihood of params under the prior
    ##
    
    prior.loglik <- function(params){
        ## this assumes improper uniform priors for everything
        ##
        ## pi(params) ~ unif(-inf,inf)
        ##
        ## so the density is a constant for all parameter values (as is the log-density)
        0
    }

    ##
    ## function to sample from proposal distribution
    ##

    rproposal <- function(params,Sigma){
        ## Sigma is a covariance matrix of the proposal distribution

        ## transform all variables to have unbounded support
        params.transformed=rep(NA,4)
        params.transformed[1:2]=qnorm(params[1:2])
        params.transformed[3:4]=log(params[3:4])

        ## random walk proposal
        params.transformed.new=rmvnorm(1,params.transformed,Sigma)
        ## back transform
        params.new=rep(NA,4)
        params.new[1:2]=pnorm(params.transformed.new[1:2])
        params.new[3:4]=exp(params.transformed.new[3:4])

        ## return
        params.new
    }


    log.mh.denom <- mmpp.loglik(params,y.list,M,delta.t)+prior.loglik(params)

    ##
    ## MCMC loop
    ##

    accept=0
    for(iter in 1:n.mcmc){
        cat(iter," ")

        ## propose
        params.new <- rproposal(params,Sigma.tune)

        ## mh step
        ## note: the proposal distribution is symmetric (multivariate normal) so the proposal dist'n cancels out!
        ##
        log.mh.num <- mmpp.loglik(params.new,y.list,M,delta.t)+prior.loglik(params.new)
        
        if(runif(1)<exp(log.mh.num-log.mh.denom)){
            params=params.new
            log.mh.denom=log.mh.num
            accept=accept+1
        }

        ## save params
        params.save[iter,]=params
########################################################
        ##
        ## Adaptive MCMC for beta
        ##
########################################################
        if (iter < adapt.max + 1 & iter / adapt.int == round(iter/adapt.int)){
            params.transformed=params.save
            params.transformed[,1:2]=qnorm(params.save[,1:2])
            params.transformed[,3:4]=log(params.save[,3:4])
            Sigma.tune <- 2.4^2/ncol(params.transformed) * var(params.transformed[1:iter, ])
        }
        
        
    }
    list(params=params.save,accept=accept)
}
        
