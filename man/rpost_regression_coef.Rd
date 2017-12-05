\name{rpost_regression_coef}
\alias{rpost_regression_coef}
\title{
  Generate posterior for regression coefecients
}
\description{
  Generate a \eqn{N(\mu, \Sigma)} random variable where
  \deqn{\Sigma = (X^tX + D^{-1})^{-1},}
  \deqn{\mu = \Sigma X^t \alpha.}
  The algorithm is \eqn{O(np^2)}; for large p it performs fast.
}
\usage{
rpost_regression_coef <- function(
  X,
  D,
  alpha,
  u = NULL)
}
\arguments{
  \item{X}{n by p matrix}
  \item{D}{p by p matrix}
  \item{alpha}{n by 1 vector}
  \item{u}{Optional. If specified, don't generate \eqn{u \sim N(0, D)}.}
}
\details{
  The algorithm is from

Bhattacharya, Anirban, Antik Chakraborty, and Bani K. Mallick. "Fast sampling with Gaussian scale mixture priors in high-dimensional regression." Biometrika (2016): asw042.

The algorithm is as follows:

1. Sample \eqn{u \sim N(0, D), d \sim N(0, I_n)}

2. Set \eqn{v = Xu + d}

3. Solve \eqn{(XDX^t + I_n)w = (\alpha - v)}

4. Return \eqn{\beta = u + DX^tw}
}
\value{
  \item{beta}{An p x 1 vector}
}

