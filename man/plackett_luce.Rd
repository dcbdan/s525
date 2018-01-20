\name{plackett_luce}
\alias{plackett_luce}
\title{
  Gibbs Sampler for posterior plackett luce parameters
}
\description{
  The prior on each parameter is given by a gamma distribution.

  The parameter expansion is given by
  Caron, Francois, and Arnaud Doucet. "Efficient Bayesian inference for
  generalized Bradley¿Terry models." Journal of Computational and
  Graphical Statistics 21.1 (2012): 174-196.


}
\usage{
plackett_luce <- function(
  rank_matrix,
  shape = 1,
  rate = 1,
  niter = 1000)
}
\arguments{
  \item{rank_matrix}{n by p matrix of integers ranging from 1,...,k where
    k = max(rank_matrix, na.rm=TRUE).
    Each row may contain NA values provided that at least 2 values of each row
    are not NA. The same integer can appear in each row of the rank matrix as well.}
  \item{shape}{scalar or k vector specifying prior shape parameter over pl parameters}
  \item{rate}{scalar or k vector specifying prior rate parameter over pl parameters}
  \item{niter}{number of gibbs sampling iterations}
}
\details{
  This function returns a niter x k matrix of the gibbs sampler.
}
\examples{
  library(coda)

  vals <- rplackett_luce(100, vs = c(100, 1, 1))
  post_vals <- plackett_luce(vals)
  plot(as.mcmc(post_vals))
}
