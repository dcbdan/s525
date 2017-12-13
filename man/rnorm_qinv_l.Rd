\name{rnorm_qinv_l}
\alias{rnorm_qinv_l}
\title{
  Sample from multivariate normal distribution }
\description{
  Sample from multivariate normal distribution with
  mean \eqn{Q^{-1}l} and covariance matrix \eqn{Q^{-1}}.
}
\usage{
rnorm_qinv_l <- function(
  n,
  Q,
  l)
}
\arguments{
  \item{n}{number of elements to generate}
  \item{Q}{p by p precision matrix. }
  \item{l}{p by 1 vector}
}
\details{
The algorithm is as follows

1. Cholesky decomposition of \eqn{Q} into \eqn{LL^t}.

2. Sample z from rnorm(p). Let \eqn{y = Lz + l}.

3. Solve for \eqn{x} in \eqn{LL^tx = y} and return.
}
\value{
  \item{beta}{An p x 1 vector if \eqn{n=1} otherwise a n by p matrix}
}
