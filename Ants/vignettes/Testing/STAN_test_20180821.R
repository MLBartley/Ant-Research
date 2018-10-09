library(rstan)
library(bayesplot)
library(ggplot2)
library(coda)
library(circular)
library(moveHMM)

rstan_options(auto_write = TRUE) #A logical scalar that controls whether a
# compiled instance of a stanmodel-class is written to the hard disk in the
# same directory as the .stan program.

options(mc.cores = parallel::detectCores())
pal <- c("firebrick","seagreen","navy") # colour palette
set.seed(1)

# Number of states
N <- 2
# transition probabilities
Gamma <- matrix(c(0.9,0.1,
                  0.1,0.9),2,2)

# initial distribution set to the stationary distribution
delta <- solve(t(diag(N)-Gamma +1), rep(1, N))
# state-dependent Gaussian means
lambda <- c(.05,.2)
nobs <- 1000
S <- rep(NA,nobs)
y <- rep(NA,nobs)
# initialise state and observation
S[1] <- sample(1:N, size=1, prob=delta)
y[1] <- rpois(1, lambda[S[1]])
# simulate state and observation processes forward
for(t in 2:nobs) {
  S[t] <- sample(1:N, size=1, prob=Gamma[S[t-1],])
  y[t] <- rpois(1, lambda[S[t]])
}
plot(y, col=pal[S], type="h")


sink("./vignettes/stan_test.stan")
cat("
//define data
data {
  int<lower=0> N; // number of states
  int<lower=1> T; // length of data set
  int<lower=0> y[T]; // observations
}

//define parameters
parameters {
  simplex[N] theta[N]; // N x N tpm
  ordered[N] lambda[N];
}

//initialize tpm with delta
transformed parameters{
  matrix[N, N] ta; //
  simplex[N] statdist; // stationary distribution
  for(j in 1:N){
      for(i in 1:N){
        ta[i,j]= theta[i,j];
      }
  }
statdist =  to_vector((to_row_vector(rep_vector(1.0, N))/
(diag_matrix(rep_vector(1.0, N)) - ta + rep_matrix(1, N, N)))) ; //

}


# Define model

model {
  vector[N] log_theta_tr[N];
  vector[N] lp;
  vector[N] lp_p1;

  // prior for lambda
  lambda[1] ~ gamma(1, 1);
  lambda[2] ~ gamma(1, 1);


  // transpose the tpm and take natural log of entries
  for (n_from in 1:N)
  for (n in 1:N)
  log_theta_tr[n, n_from] = log(theta[n_from, n]);

  // forward algorithm implementation
  for(n in 1:N) // first observation
  lp[n] = log(statdist[n]) + poisson_lpmf(y[1] | lambda[n]);
  for (t in 2:T) { // looping over observations
  for (n in 1:N) // looping over states
  lp_p1[n] = log_sum_exp(log_theta_tr[n] + lp) +
  poisson_lpmf(y[t] | lambda[n]);
  lp = lp_p1; }
  target += log_sum_exp(lp);
}
// end model
    ",fill = TRUE)
sink()

stan.data <- list(y=y, T=nobs, N=2)
fit <- stan(file="./vignettes/stan_test.stan", data=stan.data,
            refresh=5000, control = list(adapt_delta = 0.99))

lambdas <- extract(fit, pars=c("lambda"))
hist(lambdas[[1]][,1 , ], main="",xlab=expression(lambda[1]))
abline(v=lambda[1], col=pal[1], lwd=2)
hist(lambdas[[1]][,2, ],main="",xlab=expression(lambda[2]))
abline(v=lambda[2], col=pal[2], lwd=2)

## extract posterior draws
psam <- extract(fit, pars = c("theta", "lambda"))
## generate new data sets
n.sims <- dim(psam[[1]])[1]
n <- length(y)

# state sequences
ppstates <- matrix(NA, nrow = n.sims, ncol = n) # observations
ppobs <- matrix(NA, nrow = n.sims, ncol = n)
for (j in 1:n.sims) {
  theta <- psam[[1]][j, , ]
  statdist <- solve(t(diag(N) - theta + 1), rep(1, N))
  ppstates[j, 1] <- sample(1:N, size = 1, prob = statdist)
   ppobs[j, 1] <- rpois(1, lambda = psam[[2]][j, 1, ppstates[j, 1]])
  for (i in 2:length(y)) {
    ppstates[j, i] <- sample(1:N, size = 1, prob = theta[ppstates[j, i -
                                                                    1], ])
     ppobs[j, i] <- rpois(1, lambda = psam[[2]][j, 1,  ppstates[j, i]])
  } }
