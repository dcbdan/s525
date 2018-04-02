\name{rstnorm_ineq_cork}
\alias{rstnorm_ineq_cork}
\title{
  Sample from multivariate standard normal distribution with inequality constraints }
\description{
  Sample from multivariate normal distribution with
  mean \eqn{0} and covariance matrix \eqn{I_d} with the
  constraint that \eqn{M^t z \leq m} where \eqn{M} is
  given by ineq.mat and \eqn{m} is given by ineq.vec
}
\usage{
rstnorm_ineq_cork = function(
  n,
  ineq.mat,
  ineq.vec,
  thin = 1,
  burn = 0,
  init = NULL)
}
\arguments{
  \item{n}{number of elements to generate}
  \item{ineq.mat}{Inequality matrix}
  \item{ineq.vec}{Inequality vector}
  \item{thin}{Take every thinth value from the gibbs sampler}
  \item{burn}{Discard the first burn elements generated}
  \item{init}{initial z value. If not given, an SVD of ineq.mat happens.}
}
\details{
  See Schmidt, Mikkel. "Linearly constrained bayesian matrix factorization for blind source separation." Advances in neural information processing systems. 2009.
}
\examples{
# sample on a cube!
M = rbind(c(-1,1,0,0),c(0,0,-1,1))
xmin = 1
xmax = 100
ymin = 1
ymax = 100
m = c(-xmin, xmax, -ymin, ymax)
n = 1000
rr = rstnorm_ineq_cork(n, M, m)
plot(rr[,1], rr[,2])

# sample on a triangle and have the third dimension truncated positive
M = cbind(c(-1, 0, 0), c(0, -1, 0), c(1, 1, 0), c(0, 0, -1))
m = c(0, 0, 3, 0)
n = 10000
rr = rstnorm_ineq_cork(n, M, m)
plot(rr[,1], rr[,2])
hist(rr[,3])
}
