\name{applyAlong}
\alias{applyAlong}
\title{
  applyAlong
}
\description{
  Apply a function to all the submatrices along an arbitrarily chosen dimension.
}
\usage{
  applyAlong <- function(x, f, along = 1)
}
\arguments{
  \item{x}{An n dimensional array}
  \item{f}{A function that takes in a n-1 dimensional array and returns a vector or an array. }
  \item{along}{the dimension of x to apply f to. So if n = 2, along = 1, then f(x[idx,]) is applied. And if n = 5, along = 2, f(x[,idx,,,]) is applied. }
}
\details{
  This function returns an g x (o) dimensional array that is the result of applying
  the function f. Here, g is the size of the along dimension and (o) is the dimension
  of the output of f.
}
\examples{
}
