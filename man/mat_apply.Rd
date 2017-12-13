\name{mat_apply}
\alias{mat_apply}
\title{
  mat apply
}
\description{
  Apply a function to each element in an input and return an arbitrary dimensioned
  array.
}
\usage{
  mat_apply <- function(vec, fun)
}
\arguments{
  \item{vec}{vector of inputs fed into fun}
  \item{fun}{a function to call for each input of vec}
}
\details{
  This function returns an array with first dimension of length(vec) provided
  length(vec) > 1--the first dimension indexes over the length(vec) outputs.
}
\examples{
  # returns a 10 x 3 x 3 array. ret[i,,] contains i*diag(3).
  ret <- mat_apply(1:10, function(i){ i*diag(3) })
}
