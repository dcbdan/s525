\name{rnorm_qinv_l_eigen}
\alias{rnorm_qinv_l_eigen}
\title{
  Sample from multivariate normal distribution given symmetric eigen decomposition }
\description{
  Sample from multivariate normal distribution with
  mean \eqn{Q^{-1}l} and covariance matrix \eqn{Q^{-1}} where
  \eqn{Q = UDU^t}
}
\usage{
rnorm_qinv_l_eigen <- function(
  n,
  U,
  d,
  l)
}
\arguments{
  \item{n}{number of elements to generate}
  \item{U}{orthogonal matrix such that \eqn{Q = UDU^t}}
  \item{d}{p by 1 vector, \eqn{D = diag(d)} and \eqn{Q = UDU^t}}
  \item{l}{p by 1 vector}
}
\details{
  This function is useful to sample from \eqn{N((Q + sI)^{-1}l, (Q + sI)^{-1})}
  given the eigen decomposition of \eqn{Q}
}
\value{
  \item{x}{An p x 1 vector if \eqn{n=1} otherwise a n by p matrix}
}
\examples{
  W = cbind(c(10,1), c(1,10))
  ee = eigen(W, symmetric = TRUE)
  U = ee$vectors
  d = ee$values
  crossprod(U)

  solve(W)
  U %*% diag(1/d) %*% t(U)

  l = c(50,100)
  m = U %*% diag(1/d) %*% t(U) %*% l
  ret = rnorm_qinv_l_eigen(50000, U, d, l)
  mean(ret[,1] - m[1])
  mean(ret[,2] - m[2])
  sum(cov(ret) - (U %*% diag(1/d) %*% t(U)))

}
