\name{ridge_regression}
\alias{ridge_regression}
\title{
  Gibbs Sampler for Ridge Regression
}
\description{
  The ridge regression model coeffecients are given by
    \deqn{\arg\min_{\beta} \parallel Y - X\beta \parallel^2 +
      \lambda \parallel \beta \parallel^2.}
  This is the implementation of the bayesian interpretation.
  Namely,
  \deqn{[y_i|x_i, \beta, \tau] \sim N(x_i^t\beta, \tau^{-1}), i = 1,...,n ,}
  \deqn{[\beta|\eta] \sim N(0, \eta^{-1}I),}
  \deqn{\tau \sim Gamma(\frac{1}{100}, \frac{1}{100}),}
  \deqn{\sigma_\beta \sim unif(0, A)}
  where \eqn{A = 1000} and \eqn{\tau = \frac{1}{\sigma_\beta^2}}.
}
\usage{
ridge_regression <- function(
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
  \item{beta}{An niter by p matrix}
  \item{sigma}{An niter by 1 matrix}
  \item{sigma_beta}{An niter by 1 matrix}
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

gibbs_ridge <- ridge_regression(y, X, niter=10000)
shrinkage_regression_plot(gibbs_ridge$beta, y, X, main = "Ridge Prior")
}
