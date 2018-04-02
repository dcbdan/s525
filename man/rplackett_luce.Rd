\name{rplackett_luce}
\alias{rplackett_luce}
\title{
  Generate random permutations from a plackett luce distribution
}
\description{
  The plackett luce model assigns a probability distribution over permutations.
}
\usage{
rplackett_luce <- function(n, vs)
}
\arguments{
  \item{n}{number of permutations to generate}
  \item{vs}{k vector of probabilities; need not be scaled to sum to 1}
}
\details{
  This function returns a n x k matrix where each row is a permutation generated
  from a plackett luce distribution with parameters given by vs.
}
\examples{
  library(coda)

  vals <- rplackett_luce(100, vs = c(100, 1, 1))
  post_vals <- plackett_luce(vals)
  plot(as.mcmc(post_vals))
}
