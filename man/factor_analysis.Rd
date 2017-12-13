\name{factor_analysis}
\alias{factor_analysis}
\title{
  Runs gibbs sampler for a factor model.
}
\description{
  The model is as follows:
  \deqn{y_i = \Lambda \eta_i + \epsilon_i}
  \deqn{\epsilon_i \sim N_p(0, \Sigma)}
  where \eqn{\Sigma} = diag(\eqn{\sigma_1^2, ..., \sigma_p^2}).

  See Joyee Ghosh & David B. Dunson (2009) Default Prior Distributions and
  Efficient Posterior Computation in Bayesian Factor Analysis, Journal of
  Computational and Graphical Statistics, 18:2, 306-320, DOI: 10.1198/jcgs.2009.07145
  for conditional posteriors.
}
\usage{
factor_analysis <- function(
  Y,
  k,
  niter = 1000,
  shape_psi = 1/2,
  rate_psi = 1/2,
  shape_sigma2 = 1,
  rate_sigma2 = 0.2,
  nonzero_structure = NULL)
}
\arguments{
  \item{Y}{n by p matrix}
  \item{k}{number of factors}
  \item{niter}{number of iterations for the gibbs sampler to run.}
  \item{shape_psi}{shape parameter for psi. Can be a scalar or a k vector}
  \item{rate_psi}{rate parameter for psi. Can be a k vector}
  \item{shape_sigma2}{shape parameter for sigma2. Can be a p vector}
  \item{rate_sigma2}{rate parameter for sigma2. Can be a p vector}
  \item{nonzero_structure}{A boolean p x k matrix. If the i, jth spot is TRUE,
    then \eqn{\lambda_{ij}} is free. If the i, jth spot is FALSE, then
    \eqn{\lambda_{ij}} is zero. If not set, then a lower triangular matrix
    is used.
}
\value{
  \item{sigma}{An niter x p matrix of posterior values}
  \item{eta}{An niter x n x k array of posterior values}
  \item{lambda}{An niter x p x k array of posterior values}
}
