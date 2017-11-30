\name{probit_regression}
\alias{probit_regression}
\title{
  Gibbs Sampler for Probit Regression
}
\description{
  The probit regression model is given by
  \deqn{y_i|\pi_i \sim Bernoulli(\pi_i),}
  \deqn{\pi_i = \Phi(x_i^t\beta),}
  \deqn{\beta_j \sim N(0, A^2).}
  where \eqn{\Phi} is given by the gaussian CDF.

  The implemented parameter-expanded model is given by
  \deqn{y_i^* = x_i^t\beta + \epsilon_i,}
  \deqn{\epsilon_i \sim N(0, 1),}
  \deqn{y_i = I(y_i^* > 0).}

  The full conditional distributions are given by
  \deqn{[\beta|y, y^*, X] \sim N(Q^{-1}l, Q^{-1})}
  where \eqn{Q = X^tX + A^{-2}I} and \eqn{l = X^ty^*},
  \deqn{[y_i^*|y_i, \beta, X] \sim sgn(y_i, y_i^*)N(x_i^t\beta, 1)}
  where \eqn{sqn} is 1 if both arguments are of the same sign and zero otherwise.

  If A is a vector, \eqn{Q = X^tX + diag(A)}. A flat prior on the intercept
  at position 1 is thus given by A = c(Inf, rep(3, p)).
}
\usage{
probit_regression <- function(
  Y,
  X,
  niter,
  A = 3,
  init = NULL)
}
\arguments{
  \item{Y}{n by 1 vector of ones and zeros}
  \item{X}{n by p predictor matrix}
  \item{niter}{number of gibbs sampling iterations}
  \item{A}{A parameter. The default value is chosen to provide a reasonable
    range of \eqn{\pi}. If A is a vector, then different variances are given for
    each intercept. }
  \item{init}{Initial starting values for beta. If NULL, beta is set to zero.}
}
\details{
  This function returns a niter x p matrix of values where p is the second dimension
  of the predictor matrix X. The returned matrix contains all generated values of
  the gibbs sampling markov chain.
}
\examples{
  library(LearnBayes)
  library(s525)

  data(donner)
  y = donner$survival
  X = cbind(1, donner$age, donner$male)
  niter <- 2000

  gibbs_results <- probit_regression(y, X, niter)
}
