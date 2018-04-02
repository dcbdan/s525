\name{rdirichlet}
\alias{rdirichlet}
\title{
  Generate from a dirichlet distribution
}
\description{
  Generate from a dirichlet distribution
}
\usage{
rdirichlet <- function(n, shape)
}
\arguments{
  \item{n}{number of permutations to generate}
  \item{shape}{shape parameter of length k}
}
\details{
  This function returns a n x k matrix where each row is a generated value
  from the dirichlet distribution with the given shape parameter.
}


