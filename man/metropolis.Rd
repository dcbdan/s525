\name{metropolis}
\alias{metropolis}
\title{
  Metropolis MCMC
}
\description{
  The metropolis algorithm is a special case of the metropolis-hastings algorithm,
  namely where the proposal distribution is symmetric.
}
\usage{
  metropolis <- function(
    rproposal,
    prob,
    niter,
    init,
    log_prob = FALSE)
}
\arguments{
  \item{rproposal}{rproposal(previous_val) generates a value from the
    proposal distribution}
  \item{prob}{prob(val_proposed)/prob(val_t_minus_1) forms the acceptance probability}
  \item{niter}{number of iterations to perform}
  \item{init}{vector of initial values}
  \item{log_prob}{whether or not prob function specifies log probabilities}
}
\details{
  This function returns a niter x d matrix of values where d is the dimension of
  init and the dimension of each element from rproposal. The returned matrix
  contains all generated values of the metropolis walk.
}
\examples{
  vals <- metropolis(
    rproposal = function(val){ rnorm(1, mean = val, sd = 1){} },
    prob = posterior_prob,
    niter = 10000,
    init = 0)

  library(coda)
  plot(as.mcmc(vals))
}

