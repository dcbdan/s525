rnorm_qinv_l_chol <- function(
  n,
  L, # Let Q = LL^t
  l) # sample from N(Q^{-1}l, Q^{-1})
{
  rnorm_qinv_l(n, NULL, l, L)
}
