rnorm_qinv_l_eigen <- function(
  n,
  U,
  d, # a vector
  l) # sample from N(Q^{-1}l, Q^{-1}) where Q = UDU^t, U is orthogonal
{
  p <- ncol(U)
  stopifnot(nrow(U) == p)
  stopifnot(length(l) == p)
  stopifnot(length(d) == p)

  gen_x = function()
  {
    z <- rnorm(p, sd = sqrt(1/d))
    return(as.vector(U %*% z + U %*% (1/d*(crossprod(U, l)))))
  }

  if(n == 1)
  {
    return(gen_x())
  }

  return(t(replicate(n, gen_x())))
}
