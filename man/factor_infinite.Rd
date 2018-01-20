\name{factor_infinite}
\alias{factor_infinite}
\title{
  Runs gibbs sampler for a factor model with potentially infinite k
}
\description{
  The model is as follows:
  \deqn{y_i = \Lambda \eta_i + \epsilon_i}
  \deqn{\epsilon_i \sim N_p(0, \Sigma)}
  where \eqn{\Sigma} = diag(\eqn{\sigma_1^2, ..., \sigma_p^2}).

  See Bhattacharya, Anirban, and David B. Dunson.
  'Sparse Bayesian infinite factor models.' Biometrika (2011): 291-306.
}
\usage{
factor_infinite <- function(
  Y,
  k_star,
  niter = 1000,
  a_sig = 1,
  b_sig = 1,
  rho = 3,
  a1 = 1,
  a2 = 3)
}
\arguments{
  \item{Y}{n by p matrix}
  \item{k_star}{number of factors to cut off at}
  \item{niter}{number of iterations for the gibbs sampler to run.}
  \item{a_sig}{shape hyper parameter for the \eqn{\sigma^2}}
  \item{b_sig}{shape hyper parameter for the \eqn{\sigma^2}}
  \item{rho}{hyper parameter}
  \item{a1}{ hyper parameter}
  \item{a2}{ hyper parameter}
}
\value{
  \item{sigma2}{An niter x p matrix of posterior values}
  \item{lambda}{An niter x p x k_star array of posterior values}
}
