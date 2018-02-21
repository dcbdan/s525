\name{rnorm_qinv_l_chol}
\alias{rnorm_qinv_l_chol}
\title{
  Sample from multivariate normal distribution }
\description{
  Sample from multivariate normal distribution with
  mean \eqn{Q^{-1}l} and covariance matrix \eqn{Q^{-1}} where
  \eqn{Q = LL^t}.
}
\usage{
rnorm_qinv_l_chol <- function(
  n,
  L,
  l)
}
\arguments{
  \item{n}{number of elements to generate}
  \item{L}{NULL by default. If not null, Q is ignored and assumed to be LL^t.}
  \item{l}{p by 1 vector}
}
\details{
The algorithm is as follows

1. Cholesky decomposition of \eqn{Q} into \eqn{LL^t}. (This step is skipped if L is passed in).

2. Sample z from rnorm(p). Let \eqn{y = Lz + l}.

3. Solve for \eqn{x} in \eqn{LL^tx = y} and return.
}
\value{
  \item{x}{An p x 1 vector if \eqn{n=1} otherwise a n by p matrix}
}
