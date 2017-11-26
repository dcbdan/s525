\name{shrinkage_regression_plot}
\alias{shrinkage_regression_plot}
\title{
  Plot the MC chain of regression coeffecients compared to OLS.
}
\description{
  Plot the MC chain of regression coeffecients comparet to OLS using 95\% HPD
  confidence intervals.
}
\usage{
horseshoe_regression_plot <- function(
  beta,
  Y,
  X,
  main = "Horseshoe Prior",
  ylim = NULL)
}
\arguments{
  \item{beta}{n x p MCMC chain of beta coeffecients.}
  \item{Y}{Y values used to generate beta; used as input to OLS regression}
  \item{X}{X values used to generate beta; used as input to OLS regression}
  \item{main}{Title of the plot}
  \item{ylim}{y range of plot}
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
