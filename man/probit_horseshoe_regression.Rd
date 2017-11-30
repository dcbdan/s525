\name{probit_horseshoe_regression}
\alias{probit_horseshoe_regression}
\title{
  Gibbs Sampler for Probit Regression with Horseshoe Prior
}
\description{
  The probit regression model with horseshoe prior is given by
  \deqn{y_i|\pi_i \sim Bernoulli(\pi_i),}
  \deqn{\pi_i = \Phi(x_i^t\beta),}
  \deqn{[\beta_j|\lambda_j] \sim N(0, \lambda_j^2), j=2,...,p,}
  \deqn{p(\beta_1) \propto 1,}
  \deqn{[\lambda_j|A] \sim C^{+}(0, A), j=2,...,p,}
  \deqn{A \sim Uniform(0, 10).}
  where \eqn{\Phi} is given by the gaussian CDF.

  The implemented parameter-expanded model is given by
  \deqn{y_i^* = x_i^t\beta + \epsilon_i,}
  \deqn{\epsilon_i \sim N(0, 1),}
  \deqn{y_i = I(y_i^* > 0).}
  The half-Cauchy parameter expansion is also used; given by
  \deqn{[\eta_j|\gamma_j] \sim Gamma(\frac{1}{2}, \gamma_j),}
  \deqn{[\gamma_j] \sim Gamma(\frac{1}{2}, \frac{1}{A^2})}
  and \eqn{\eta_j = \lambda_j^{-2} },
    \eqn{\tau_A = A^{-2} },
  The full conditionals are given by:
  \deqn{[y_i^*|y_i, \beta, X] \sim sgn(y_i, y_i^*)N(x_i^t\beta, 1)}
  where \eqn{sqn} is 1 if both arguments are of the same sign and zero otherwise,
  \deqn{
    [\beta|Y^*, X, \eta] \sim \mathcal{N}(Q^{-1}l, Q^{-1})}
  where \eqn{Q = X'X+diag(0, 1/\eta_2, ..., 1/\eta_p)} and \eqn{l = X'Y^*},
  \deqn{
    [\eta_j|\beta_j, \gamma_j] \sim
      \mathrm{exp}(\frac{\beta_j^2}{2} + \gamma_j),}
  \deqn{
    [\gamma_j|\eta_j, \tau_A] \sim
      \mathrm{exp}(\eta_j + \tau_A),}
  \deqn{
    [\tau_A|\gamma] \sim
      \mathrm{Gamma}(
        \frac{p - 2}{2}, \sum \gamma_i)\mathrm{I}_{(\frac{1}{100}, \infty)}.}
}
\usage{
probit_horseshoe_regression <- function(
  Y,
  X,
  niter,
  init = NULL)
}
\arguments{
  \item{Y}{n by 1 vector of ones and zeros}
  \item{X}{n by p predictor matrix, where p > 1 and the first column of X is all 1.}
  \item{niter}{number of gibbs sampling iterations}
  \item{init}{Initial starting values for beta. If NULL, beta is set to zero.}
}
\details{
  This function returns a niter x p matrix of values where p is the second dimension
  of the predictor matrix X. The returned matrix contains all generated values of
  the gibbs sampling markov chain.
}
\examples{
 print("TODO")
}
