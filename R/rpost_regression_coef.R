rpost_regression_coef <- function(
  X,
  D,
  alpha,
  u = NULL)
{
  n <- nrow(X)
  p <- ncol(X)

  stopifnot(nrow(D) == p)
  stopifnot(ncol(D) == p)
  stopifnot(length(alpha) == n)

  if(is.null(u))
  {
    u <- mvrnorm(1, mu = rep(0, p), Sigma = D)
  }

  d <- rnorm(n)

  v <- X%*%u + d

  chq <- chol(X %*% D %*% t(X) + diag(n))
  w <- backsolve(chq, forwardsolve(t(chq), alpha - v))

  ret <- u + D %*% t(X) %*% w
}
