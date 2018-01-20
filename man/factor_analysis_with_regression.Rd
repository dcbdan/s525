\name{factor_analysis_with_regression}
\alias{factor_analysis_with_regression}
\title{
  Runs gibbs sampler for a factor model with regression on the factors.
}
\description{
  The model is as follows:
  \deqn{y_i = \alpha + \Lambda \eta_i + \epsilon_i,}
  \deqn{\eta_i = Bx_i + \delta_i}
  \deqn{\epsilon_i \sim N_p(0, \tau^{-1})}
  \deqn{\delta_i \sim N_k(0, I)}
  where \eqn{\tau} = diag(\eqn{\tau_1, ..., \tau_p}).

  See Joyee Ghosh & David B. Dunson (2009) Default Prior Distributions and
  Efficient Posterior Computation in Bayesian Factor Analysis, Journal of
  Computational and Graphical Statistics, 18:2, 306-320, DOI: 10.1198/jcgs.2009.07145.
}
\usage{
factor_analysis_with_regression <- function(
  Y,
  X,
  k,
  niter = 1,
  shape_psi = 1/2,
  rate_psi = 1/2,
  shape_tau = 1,
  rate_tau = 0.2,
  coef_multiplier = 10,
  nonzero_structure = NULL)
}
\arguments{
  \item{Y}{n by p matrix}
  \item{X}{n by f matrix}
  \item{k}{number of factors}
  \item{niter}{number of iterations for the gibbs sampler to run.}
  \item{shape_psi}{shape parameter for psi. Can be a scalar or a k vector}
  \item{rate_psi}{rate parameter for psi. Can be a k vector}
  \item{shape_tau}{shape parameter for sigma2. Can be a p vector}
  \item{rate_tau}{rate parameter for sigma2. Can be a p vector}
  \item{nonzero_structure}{A boolean p x k matrix. If the i, jth spot is TRUE,
    then \eqn{\lambda_{ij}} is free. If the i, jth spot is FALSE, then
    \eqn{\lambda_{ij}} is zero. If not set, then a lower triangular matrix
    is used.}
}
\value{
  \item{alpha}{An niter x p matrix of posterior values}
  \item{lambda}{An niter x p x k array of posterior values}
  \item{tau}{An niter x p matrix of posterior values}
  \item{eta}{An niter x n x k array of posterior values}
  \item{B}{An niter x k x f array of posterior values}
}
