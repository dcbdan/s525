rplackett_luce <- function(n, vs)
{
  # return a n by k matrix. Each row contains
  # a permutation
  k <- length(vs)
  t(replicate(n, sample.int(k, prob = vs)))
}
