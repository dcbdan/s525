rnorm_qinv_l <- function(
  n,
  Q,
  l,
  L = NULL) # sample from N(Q^{-1}l, Q^{-1})
{
  if(length(l) == 1 & is.null(L))
  {
    Q = as.matrix(Q)
  }

  # just leaving this function here for now
  chol_solve <- function(L, y)
  {
    # solve for x in LL'x = y
    backsolve(
      t(L),
      forwardsolve(L, y))
  }
  ###############################################

  if(is.null(L))
  {
    p <- ncol(Q)
    stopifnot(nrow(Q) == p)
    stopifnot(length(l) == p)

    L <- t(chol(Q))
  }
  else
  {
    p <- ncol(L)
  }

  z <- matrix(rnorm(p*n), nrow = p, ncol = n)

  y <- L %*% z + as.vector(l)

  if(n == 1)
  {
    return(as.numeric(chol_solve(L, y)))
  }
  else
  {
    return(t(chol_solve(L, y)))
  }
}
