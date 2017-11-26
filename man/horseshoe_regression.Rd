\name{horseshoe_regression}
\alias{horseshoe_regression}
\title{
  Gibbs Sampler for Horseshoe Regression
}
\description{
  The horseshoe regression model is given by
  \deqn{[y_i|x_i, \beta, \sigma] \sim N(x_i^t\beta, \sigma^2), i = 1,...,n ,}
  \deqn{[\beta_j|\sigma, \lambda_j] \sim N(0, \sigma^2\lambda_j^2),}
  \deqn{[\lambda_j|A] \sim C^{+}(0, A),}
  \deqn{A \sim Uniform(0, 10),}
  \deqn{p(\sigma^2) \propto \frac{1}{\sigma^2}.}
  The half-Cauchy parameter expansion is used; given by
  \deqn{[\eta_j|\gamma_j] \sim Gamma(\frac{1}{2}, \gamma_j),}
  \deqn{[\gamma_j] \sim Gamma(\frac{1}{2}, \frac{1}{A^2}).}
  Let \eqn{\eta_j = \lambda_j^{-2} }, \eqn{\tau_A = A^{-2} },
    \eqn{\tau = \frac{1}{\sigma^2}} and
    \eqn{\Lambda = diag(\eta_1, ..., \eta_p)}.
  The full conditionals are given by:
  \deqn{
    [\beta|Y, X, \eta, \tau] \sim
      \mathcal{N}(
        (X'X+\Lambda)^{-1}X'Y,
        \tau^{-1}(X'X+\Lambda)^{-1}),}
  \deqn{
    [\eta_j|\beta_j, \gamma_j, \tau] \sim
      \mathrm{exp}(\frac{\tau\beta_j^2}{2} + \gamma_j),}
  \deqn{
    [\gamma_j|\eta_j, \tau_A] \sim
      \mathrm{exp}(\eta_j + \tau_A),}
  \deqn{
    [\tau_A|\gamma] \sim
      \mathrm{Gamma}(
        \frac{p - 1}{2}, \sum \gamma_i)\mathrm{I}_{(\frac{1}{100}, \infty)},}
  \deqn{
    [\tau|Y, X, \beta, \eta] \sim
      \mathrm{Gamma}(
        \frac{n + p}{2}, \frac{(y-X\beta)'(y-X\beta) + \beta'\Lambda\beta}{2}).}
}
\usage{
horseshoe_regression <- function(
  Y,
  X,
  niter = 10000)
}
\arguments{
  \item{Y}{n by 1}
  \item{X}{n by p predictor matrix}
  \item{niter}{number of gibbs sampling iterations}
}
\details{
  This function returns the generated parameters from the gibbs sampling markov
  chain.
}
\value{
  \item{beta}{An niter x p matrix}
  \item{lambda}{An niter x p matrix}
  \item{sigma}{An niter x 1 matrix}
}
\examples{
# Load the data
 prostate.data = "https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data"
prostate = read.table(file = "prostate.data", sep="", header = TRUE)
# Training data:
prostate_train = prostate[which(prostate$train),-10]
# Testing data:
prostate_test = prostate[which(!prostate$train),-10]
# Response:
y = prostate_train$lpsa
# Center and scale the data:
y = scale(y)
# And the predictors
X = scale(prostate_train[,names(prostate_train) != "lpsa"])

gibbs_hs <- horseshoe_regression(y, X, niter=10000)
shrinkage_regression_plot(gibbs_hs$beta, y, X)
}
