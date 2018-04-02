rdirichlet <- function(n, shape)
{
  len <- length(shape)
  rvs <- matrix(
    rgamma(n = len*n, shape = shape),
    nrow = n,
    ncol = len,
    byrow = TRUE)
  ret <- rvs / apply(rvs, 1, sum)

  if(n == 1)
    return(as.vector(ret))
  else
    return(ret)
}
