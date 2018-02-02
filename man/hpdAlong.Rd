\name{hpdAlong}
\alias{hpdAlong}
\title{
  hpdAlong
}
\description{
  Take the hpd interval using replicates along the along dimension.
}
\usage{
  hpdAlong <- function(x, along, prob = 0.95)
}
\arguments{
  \item{x}{An n dimensional array}
  \item{along}{the dimension of x to use replicates along}
  \item{prob}{the target probability content of the intervals}
}
\details{
  This function returns two n - 1 dimensional arrays that is the result of
  taking the HPDinterval along each of the vectors.
}

