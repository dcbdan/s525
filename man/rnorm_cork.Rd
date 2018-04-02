\name{rnorm_cork}
\alias{rnorm_cork}
\title{
  Sample from multivariate normal distribution with
  potential inequality constraints and potential
  equality constraints}
\description{
  Sample \eqn{X} from a multivariate normal distribution with
  mean \eqn{\mu}, covariance matrix \eqn{\Sigma}, subject to the constraints
  that \eqn{M^tX \leq m} and \eqn{A^tX \leq a}.
}
\usage{
rnorm_cork = function(
  n,
  mean.vec, #   Sample X ~ normal with mean.vec, cov.mat subject
  cov.mat,  #   to the constraint that MX <= m, AX = a
  ineq.mat = NULL, #M
  ineq.vec = NULL, #m
  eq.mat = NULL,   #A
  eq.vec = NULL)   #a
}
\arguments{
  \item{n}{number of elements to generate}
  \item{mean.vec}{Mean before constraints}
  \item{cov.mat}{Covariance matrix before constrainnts}
  \item{ineq.mat}{Inequality matrix. If NULL, no inequality constraints are computed}
  \item{ineq.vec}{Inequality vector. If NULL, no inequality constraints are computed}
  \item{eq.mat}{Equality matrix. If NULL, no equality constraints are computed}
  \item{eq.vec}{Equality vector. If NULL, no equality constraints are computed}
}
\details{
  See Schmidt, Mikkel. "Linearly constrained bayesian matrix factorization for blind source separation." Advances in neural information processing systems. 2009.
}
\examples{
M = rbind(c(-1,1,0,0),c(0,0,-1,1))
xmin = 0
xmax = 1
ymin = 0
ymax = 1
m = c(-xmin, xmax, -ymin, ymax)
n = 5000

A = array(c(1,-1), dim = c(2,1))
a = 0

rho = -0.9
mean.vec = c(2,2)
cov.mat = cbind(c(1, rho), c(rho, 1))

par(mfrow = c(2,2))

# no constraints
rr = rnorm_cork(
  n,
  mean.vec,
  cov.mat)
plot(rr[,1], rr[,2])

# just inequality constraints
rr = rnorm_cork(
  n,
  mean.vec,
  cov.mat,
  ineq.mat = M,
  ineq.vec = m)
plot(rr[,1], rr[,2])

# just equality constraints
rr = rnorm_cork(
  n,
  mean.vec,
  cov.mat,
  eq.mat = A,
  eq.vec = a)
plot(rr[,1], rr[,2])

# both equality and inequality constraints
rr = rnorm_cork(
  n,
  mean.vec,
  cov.mat,
  ineq.mat = M,
  ineq.vec = m,
  eq.mat = A,
  eq.vec = a)
plot(rr[,1], rr[,2])

# Another
nN = 10
mean.vec = rep(0, nN)
cov.mat = diag(nN)

M = array(0, c(nN - 1, nN))
diag(M[,1:(nN-1)]) = 1
diag(M[,2:(nN-0)]) = -1
m = rep(0, nrow(M))

# set first value to zero and last value to one
A = array(0, c(2, nN))
A[1,1] = 1
A[2,nN] = 1
a = c(0, 1)

# case 0:
# simulate null, null
rho = 0
oo = rnorm_cork(100, c(1, 1), cbind(c(1, rho), c(rho, 1)))
plot(oo[,1], oo[,2])

# case 1:
# simulate with inequalities
oo = rnorm_cork(100, mean.vec, cov.mat, t(M), m)
stopifnot(all(sapply(2:nN, function(idx){ all(oo[,idx-1] < oo[,idx]) })))
xs = range(oo)
plot(xs, xs, type = "l")
points(oo[,nN-1], oo[,nN])

# case 2:
# simulate with equalities
oo = rnorm_cork(100, mean.vec, cov.mat, eq.mat = t(A), eq.vec = a)
stopifnot(all(oo[,1] == 0 & oo[,nN] == 1))
par(mfrow = c(2,2))
plot(oo[1,])
plot(oo[2,])
plot(oo[3,])
plot(oo[4,])

# case 3:
# simulate with inequalities and equalities
oo = rnorm_cork(100, mean.vec, cov.mat, t(M), m, t(A), a)
stopifnot(all(sapply(2:nN, function(idx){ all(oo[,idx-1] < oo[,idx]) })))
stopifnot(all(oo[,1] == 0 & oo[,nN] == 1))
par(mfrow = c(2,2))
plot(oo[1,])
plot(oo[2,])
plot(oo[3,])
plot(oo[4,])
}
