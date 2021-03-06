% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.DT.cov.troph_function.R
\name{sim.DT.cov.troph}
\alias{sim.DT.cov.troph}
\title{This function allows you to simulate data from the following model: 
P ~ function of covariates, changes over time
X_t ~ Markov Chain ( X_t-1 , P ) [two states]
Y_t ~ Pois (lambda_{X_t})}
\usage{
sim.DT.cov.troph(tmax, delta.t, start.state = 1, alpha.beta, lambda,
  num.locations = 1)
}
\arguments{
\item{tmax, }{delta.t, start.state, alpha.beta, lambda, num.location}
}
\value{
#Output: (1) - state: (1, 2) unobserved two state process 
       (2) - inter.persec: (0, ...) 'observed' number of interactions per 1 second
       (3) - cumu.inter: (0, ...) cumulative count of interactions
       (4) - delta.t
       (5) - bin.sec: (0, ...) time (in seconds) binned by delta.t
       (6) - 1x1 visual of cumlative counts over time, colored by state, 
             separated by location if applicable
}
\description{
This function allows you to simulate data from the following model: 
P ~ function of covariates, changes over time
X_t ~ Markov Chain ( X_t-1 , P ) [two states]
Y_t ~ Pois (lambda_{X_t})
}
\examples{
alpha.beta = c(.005, .0001, .005, .0001)
lambda = k = c(1, 4)
delta.t = 1 #needs to be 1, else observations dependent on time
sim = sim.DT.troph(7200, delta.t, start.state = 1, alpha.beta, lambda, num.location = 1)

}
\keyword{HMM,}
\keyword{Poisson}
\keyword{Process,}
\keyword{covariate}
\keyword{simulation,}
