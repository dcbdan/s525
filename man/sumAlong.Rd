\name{sumAlong}
\alias{sumAlong}
\title{
  sumAlong
}
\description{
  Sum along a chosen dimension of an arbitrarily dimensioned array.
}
\usage{
  sumAlong <- function(x, along)
}
\arguments{
  \item{x}{An n dimensional array}
  \item{along}{the dimension of x to sum along}
}
\details{
  # The method of the source code came from Hadley Wickham:
  # https://stackoverflow.com/questions/14500707/select-along-one-of-n-dimensions-in-array

  This function returns an n - 1 dimensional array that is the result of
  summing along the along dimension.
}
\examples{
  nx = 3
  ny = 3
  nz = 3
  x = array(1:(nx*ny*nz), dim = c(nx, ny, nz))

  print(x)
  sumAlong(x, 3)
}
