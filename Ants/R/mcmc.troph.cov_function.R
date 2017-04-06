#' MCMC Estimation function for Trophallaxis data with Covariate(s)
#'
#' The purpose of this function is to find MCMC generated estimates of 
#' (1) - the state (X_t = high/low troph rates) of the colony at time t
#' (2) - the specific rates of interaction (lambda) of each state 
#' (3) - the probability of moving from one state to another (P matrix)
#' This model will take in observed interaction data and covariates. Currently
#' time since last forager ant entered the chamber is our covariate.
#'
#' @param y.data, ant.file, inout.file, title, a, b, theta, states, n.mcmc, cov, mu.cov, sig.cov, tau, delta.t
#' @return  (1) - estimates of X, lambda, P
#'          (2) - 2x2 visual of estimates over time (runs)
#' @export
#' @examples
#' theta = matrix(data = c(90, 10, 10, 90), nrow = 2, ncol = 2, byrow = T) 
#' mu.all = c(2, -1, -0.000004)
#' sig.all = matrix(data = c(0.2, 0, 0, 
#'                          0, 0.2, 0, 
#'                          0, 0, 0.0002), nrow = 3, ncol = 3, byrow = T)
#'                          
#' in.out =  read.csv('Data/Colony_1_in&out_high_density_4hrs.csv')
#' in.out = in.out[which(in.out$Action == 'enter'),]
#' in.out = in.out[order(in.out$time), ]
#' 
#' #Create vector of covariates - time since last ant entered
#' cov = rep(300, in.out$time[1])
#' for(i in 2:nrow(in.out)){
#'    cov = c(cov, 0:(in.out$time[i] - in.out$time[i - 1] ))
#'  }
#'  
#' cov = c(cov, 0:(14400 - in.out$time[nrow(in.out)]))
#' out.high = mcmc.troph.cov(data = high.y, title = 'High Density', a = 5, b = 2, 
#' theta = theta, states = 2, n.mcmc = 3000, cov = cov, mu.cov = mu.all, sig.cov = sig.all)


mcmc.troph.cov = function(y.data, ant.file, inout.file, title, a = 5, 
    b = 2, theta, states = 2, n.mcmc, cov, mu.cov, sig.cov, tau, delta.t, 
    hours) {
    data = y.data
    Time = length(data)
    n = states
    delta = rep(1/n, n)
    
    # needed for final graphic
    location = ant.file$Location
    start = ant.file$start_time
    start = sort(start)
    cov.time = inout.file$time
    int.num = length(start)
    maxtime = hours * 60 * 60
    
    
    if (length(unique(location)) != 1) {
        # Separate Low Density by Location
        
        low.1 = ant.file[which(ant.file$Location == 1), ]
        start2 = low.1$start_time
        
        low.4 = ant.file[which(ant.file$Location == 4), ]
        low.4 = low.4[order(low.4$start_time), ]
        start3 = low.4$start_time
    }
    
    # homes Build Homes for X(1:T), lambda(1:n), and P(nXn) and gam
    # vectors
    
    X.param = matrix(data = rep(NA, Time * n.mcmc), nrow = Time, ncol = n.mcmc, 
        byrow = T)
    
    lambda.param = matrix(data = rep(NA, n * n.mcmc), nrow = n, ncol = n.mcmc, 
        byrow = T)
    
    # P.param = matrix(data = rep(NA, n * n * n.mcmc), nrow = n * n,
    # ncol = n.mcmc, byrow = T)
    
    ## P matrix now varies over time, need new homes
    
    P.11.param = matrix(NA, nrow = Time, ncol = n.mcmc)
    P.12.param = matrix(NA, nrow = Time, ncol = n.mcmc)
    P.21.param = matrix(NA, nrow = Time, ncol = n.mcmc)
    P.22.param = matrix(NA, nrow = Time, ncol = n.mcmc)
    
    gam = matrix(NA, nrow = Time, ncol = n)
    
    alph.beta.params = matrix(NA, nrow = 3, ncol = n.mcmc)
    
    # alpha.param = matrix(NA, nrow = 1, ncol = n.mcmc) betas.param =
    # matrix(NA, nrow = 2, ncol = n.mcmc)
    
    
    ## Initialize parameters
    
    X.param[, 1] = rep(1, Time)
    
    lambda.param[, 1] = rgamma(n = n, shape = a, rate = b)
    
    P.matrix = matrix(data = theta/100, nrow = n, ncol = n, byrow = T)
    
    # P.param[, 1] = as.vector(t(P.matrix)) #holds all P.parameter
    # values over runs
    
    P.11.param[, 1] = P.matrix[1, 1]
    P.12.param[, 1] = P.matrix[1, 2]
    P.21.param[, 1] = P.matrix[2, 1]
    P.22.param[, 1] = P.matrix[2, 2]
    
    alph.beta.params[, 1] = rnorm(3, mean = mu.cov, sd = sqrt(diag(sig.cov)))
    
    ########### Metropolis Hastings/Gibbs Updates
    
    for (l in 2:n.mcmc) {
        
        m = matrix(data = rep(0, n * n), nrow = n, ncol = n)
        # number states going from i to j, refreshes every run
        
        
        
        P.matrix = matrix(data = c(P.11.param[1, l - 1], P.12.param[1, 
            l - 1], P.21.param[1, l - 1], P.22.param[1, l - 1]), nrow = n, 
            ncol = n, byrow = T)
        
        ## P matrix w/ alpha, betas - MH Algorithm
        
        log.fullcond = function(param) {
            
            P.mult = 0
            for (i in 2:Time) {
                P.mult = P.mult + log(P.matrix[X.param[i - 1, l - 
                  1], X.param[i, l - 1]])
            }
            out = P.mult + dnorm(param[1], mean = mu.cov[1], sd = sqrt(diag(sig.cov)[1]), 
                log = T) + dnorm(param[2], mean = mu.cov[2], sd = sqrt(diag(sig.cov)[2]), 
                log = T) + dnorm(param[3], mean = mu.cov[3], sd = sqrt(diag(sig.cov)[3]), 
                log = T)
            
            ## need code to avoid output of -Inf!!
            if (out == -Inf) {
                return(1e-09)
            } else (return(out))
            
            
            
        }
        
        # proposal
        
        if (l%%100 == 0) {
            sigma = c(0, 0, 0)
            for (a in 1:3) {
                sigma[a] = ((2.38^2)/3) * var(alph.beta.params[a, 
                  1:(l - 1)])
            }
            
            if ((sigma[1] != 0) & (sigma[2] != 0) & (sigma[3] != 0)) {
                tau = sigma
            }
            
        }
        
        proposal = rnorm(3, mean = alph.beta.params[, l - 1], sd = tau)
        
        
        
        # accept/reject all params - Block update
        
        prob.all = exp(log.fullcond(proposal) - log.fullcond(alph.beta.params[, 
            l - 1]))
        
        if (runif(1) < prob.all) {
            alph.beta.params[, l] = proposal
        } else {
            alph.beta.params[, l] = alph.beta.params[, l - 1]
        }
        
        
        # Now use priors to calculate P matrix values
        alpha = alph.beta.params[1, l]
        beta.0 = alph.beta.params[2, l]
        # beta.1 = 0
        beta.1 = alph.beta.params[3, l]
        
        
        for (i in 1:Time) {
            
            P.matrix[1, 2] = (exp(beta.0 + beta.1 * cov[i]))/(1 + 
                exp(beta.0 + beta.1 * cov[i]))
            P.matrix[1, 1] = 1 - P.matrix[1, 2]
            P.matrix[2, 2] = (exp(alpha))/(1 + exp(alpha))
            P.matrix[2, 1] = 1 - P.matrix[2, 2]
            
            P.12.param[i, l] = P.matrix[1, 2]
            P.11.param[i, l] = P.matrix[1, 1]
            P.22.param[i, l] = P.matrix[2, 2]
            P.21.param[i, l] = P.matrix[2, 1]
            
        }
        
        ## X Parameters
        
        ## For X at t = 1
        
        P.matrix = matrix(data = c(P.11.param[1, l - 1], P.12.param[1, 
            l - 1], P.21.param[1, l - 1], P.22.param[1, l - 1]), nrow = n, 
            ncol = n, byrow = T)
        for (k in 1:n) {
            gam[1, k] = lambda.param[k, l - 1]^data[1] * exp(-lambda.param[k, 
                l - 1]) * delta[k] * P.matrix[k, X.param[2, l - 1]]
        }
        
        X.param[1, l] = sample(x = (1:n), size = 1, prob = gam[1, 
            ])
        
        m[X.param[1, l], X.param[1, l]] = m[X.param[1, l], X.param[1, 
            l]] + 1
        
        
        ## For X from t = 2 to t = Time - 1
        
        for (t in 2:(Time - 1)) {
            
            P.matrix = matrix(data = c(P.11.param[t, l - 1], P.12.param[t, 
                l - 1], P.21.param[t, l - 1], P.22.param[t, l - 1]), 
                nrow = n, ncol = n, byrow = T)
            
            for (k in 1:n) {
                gam[t, k] = (lambda.param[k, l - 1]^data[t]) * exp(-lambda.param[k, 
                  l - 1]) * P.matrix[X.param[t - 1, l - 1], k] * P.matrix[k, 
                  X.param[t + 1, l - 1]]
            }
            
            X.param[t, l] = sample(x = (1:n), 1, prob = gam[t, ])
            
            m[X.param[t - 1, l], X.param[t, l]] = m[X.param[t - 1, 
                l], X.param[t, l]] + 1
        }
        
        ## For X at t = Time
        
        P.matrix = matrix(data = c(P.11.param[Time, l - 1], P.12.param[Time, 
            l - 1], P.21.param[Time, l - 1], P.22.param[Time, l - 
            1]), nrow = n, ncol = n, byrow = T)
        
        for (k in 1:n) {
            gam[Time, k] = lambda.param[k, l - 1]^data[Time] * exp(-lambda.param[k, 
                l - 1]) * P.matrix[X.param[(Time - 1), l - 1], k]
        }
        
        X.param[Time, l] = sample(x = 1:n, 1, prob = gam[Time, ])
        
        m[X.param[Time - 1, l], X.param[Time, l]] = m[X.param[Time - 
            1, l], X.param[Time, l]] + 1
        
        
        # Lambda
        for (h in 1:n) {
            lambda.param[h, l] = rgamma(n = 1, shape = sum(data[which(X.param[, 
                l] == h)]) + a, rate = sum(m[h, ]) + b)
        }
        
        
        
    }
    
    
    ## Rescale Lambda parameters into per minute segments (instead of
    ## delta.t time segments)
    
    lambda.scale = lambda.param/delta.t * 60
    
    
    
    ## Compile the Estimates
    
    ## X1:XT, Lambda, Pmatrix, alpha/beta0/beta1
    
    # homes
    X.est = matrix(NA, nrow = Time, ncol = 1)
    lambda.est = matrix(data = rep(NA, n), nrow = n, ncol = 1)
    
    P.11.est = matrix(NA, nrow = Time, ncol = 1)
    P.12.est = matrix(NA, nrow = Time, ncol = 1)
    P.21.est = matrix(NA, nrow = Time, ncol = 1)
    P.22.est = matrix(NA, nrow = Time, ncol = 1)
    
    alph.beta.est = matrix(NA, nrow = 3, ncol = 1)
    
    for (t in 1:Time) {
        X.est[t, 1] = mean(X.param[t, ])
    }
    
    for (i in 1:n) {
        lambda.est[i, 1] = mean(lambda.scale[i, ])
    }
    
    for (t in 1:Time) {
        # P.est[l, 1] = mean(P.param[l, ])
        P.11.est[t, 1] = mean(P.11.param[t, ])
        P.12.est[t, 1] = mean(P.12.param[t, ])
        P.21.est[t, 1] = mean(P.21.param[t, ])
        P.22.est[t, 1] = mean(P.22.param[t, ])
        
    }
    
    
    # P.est.matrix = matrix(data = c(P.est[, 1]), nrow = n, ncol = n,
    # byrow = T)
    P.est.matrix = matrix(data = c(mean(P.11.est), mean(P.12.est), 
        mean(P.21.est), mean(P.22.est)), nrow = n, ncol = n, byrow = T)
    
    
    for (i in 1:3) {
        alph.beta.est[i, ] = mean(alph.beta.params[i, ])
    }
    
    # plot the estimation runs.
    
    # lambda
    par(mfrow = c(2, 2), oma = c(0, 0, 2, 0) + 1, mar = c(1, 1, 1, 
        1) + 3)
    plot(0, 0, xlab = "MCMC Runs", ylab = "Lambda (scaled to per 60 seconds)", 
        ylim = c(0, max(lambda.scale)), xlim = c(0, n.mcmc), type = "n", 
        cex.lab = 1)
    for (i in 1:n) {
        lines(1:n.mcmc, lambda.scale[i, ], col = i)
    }
    
    
    # P
    t = sample(1:Time, 1)
    plot(0, 0, xlab = "MCMC Runs", ylab = "Single P", ylim = c(0, 
        1), xlim = c(0, n.mcmc), type = "n", cex.lab = 1)
    # for(i in 1:(n * n)){
    lines(1:n.mcmc, P.11.param[t, ], col = "red")
    lines(1:n.mcmc, P.12.param[t, ], col = "blue")
    lines(1:n.mcmc, P.21.param[t, ], col = "black")
    lines(1:n.mcmc, P.22.param[t, ], col = "green")
    # }
    
    # Single X
    X = X.param[sample(1:Time, 1), ]
    plot(0, 0, xlab = "MCMC Runs", ylab = "Single X", ylim = c(0, 
        max(X)), xlim = c(0, n.mcmc), type = "n", cex.lab = 1)
    lines(1:n.mcmc, X, col = "red")
    
    # States over time
    plot(X.est, type = "l", lwd = 3, cex.lab = 1)
    
    title(main = title, outer = T)
    
    # alpha beta1 beta0 plot
    
    plot(0, 0, xlab = "MCMC Runs", ylab = "Alpha", ylim = c(min(alph.beta.params[1, 
        ]), max(alph.beta.params[1, ])), xlim = c(0, n.mcmc), type = "n", 
        cex.lab = 1)
    lines(1:n.mcmc, alph.beta.params[1, ], col = "red")
    
    plot(0, 0, xlab = "MCMC Runs", ylab = "Beta_0", ylim = c(min(alph.beta.params[2, 
        ]), max(alph.beta.params[2, ])), xlim = c(0, n.mcmc), type = "n", 
        cex.lab = 1)
    lines(1:n.mcmc, alph.beta.params[2, ], col = "green")
    
    plot(0, 0, xlab = "MCMC Runs", ylab = "Beta_1", ylim = c(min(alph.beta.params[3, 
        ]), max(alph.beta.params[3, ])), xlim = c(0, n.mcmc), type = "n", 
        cex.lab = 1)
    lines(1:n.mcmc, alph.beta.params[3, ], col = "blue")
    
    ######################################################### Fancy Plots with Background Colors
    par(mfrow = c(1, 1))
    
    if (length(unique(location)) == 1) {
        
        ## High Density - 4 Hours
        plot(start, 1:int.num, main = "High", xlab = "Seconds", ylab = "Cumulative Interaction Count", 
            xlim = c(0, maxtime))
        ## plot(one.day,1:length(one.day),main=day,xlab='Minutes')
        states = X.est  #from code above
        rr = rle(states[, 1])
        rr$values = round(rr$values, digits = 0)
        embedded.chain = rr$values
        cs = c(0, cumsum(rr$lengths)) * delta.t - delta.t
        cols = c("#bc535644", "#538bbc44")
        for (j in 1:length(embedded.chain)) {
            rect(cs[j], 0, cs[j + 1], int.num, col = cols[embedded.chain[j]], 
                density = NA)
            
        }
        points(cov.time, rep(0, length(cov.time)), pch = 8, col = "#53bc84")
    } else {
        # Low Density - 4 Hours
        
        plot(start, 1:int.num, main = "Low", xlab = "Seconds", ylab = "Cumulative Interaction Count", 
            xlim = c(0, maxtime))
        states = X.est
        rr = rle(states[, 1])
        rr$values = round(rr$values, digits = 0)
        embedded.chain = rr$values
        cs = c(0, cumsum(rr$lengths)) * delta.t - delta.t
        cols = c("#bc535644", "#538bbc44")
        for (j in 1:length(embedded.chain)) {
            rect(cs[j], 0, cs[j + 1], int.num, col = cols[embedded.chain[j]], 
                density = NA)
        }
        points(cov.time, rep(0, length(cov.time)), pch = 8, col = "#53bc84")
        
        # Low Density - Location 1
        
        plot(start2, 1:length(start2), main = "Low, Loc 1", xlab = "Seconds", 
            ylab = "Cumulative Interaction Count", xlim = c(0, maxtime))
        states = X.est
        rr = rle(states[, 1])
        rr$values = round(rr$values, digits = 0)
        embedded.chain = rr$values
        cs = c(0, cumsum(rr$lengths)) * delta.t - delta.t
        cols = c("#bc535644", "#538bbc44")
        for (j in 1:length(embedded.chain)) {
            rect(cs[j], 0, cs[j + 1], length(start2), col = cols[embedded.chain[j]], 
                density = NA)
        }
        points(cov.time, rep(0, length(cov.time)), pch = 8, col = "#53bc84")
        
        # Low Density - Location 4
        plot(start3, 1:length(start3), main = "Low, Loc 4", xlab = "Seconds", 
            ylab = "Cumulative Interaction Count", xlim = c(0, maxtime))
        states = X.est
        rr = rle(states[, 1])
        rr$values = round(rr$values, digits = 0)
        embedded.chain = rr$values
        cs = c(0, cumsum(rr$lengths)) * delta.t - delta.t
        cols = c("#bc535644", "#538bbc44")
        for (j in 1:length(embedded.chain)) {
            rect(cs[j], 0, cs[j + 1], length(start3), col = cols[embedded.chain[j]], 
                density = NA)
        }
        points(cov.time, rep(0, length(cov.time)), pch = 8, col = "#53bc84")
        
    }
    
    list(X.est = X.est, lambda.est = lambda.est, P.est = P.est.matrix, 
        P11.est = P.11.est, P12.est = P.12.est, P21.est = P.21.est, 
        P22.est = P.22.est, alph.beta.est = alph.beta.est)
    
}
