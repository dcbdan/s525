probit_regression <- function(
  Y,
  X,
  niter,
  A = 3)
{
  n <- nrow(X); p <- ncol(X); n0 <- sum(Y==0); n1 <- sum(Y==1)

  betas <- array(0, c(niter, p))
  beta_val <- numeric(p) #rnorm(p, mean = 0, sd = A)

  XTX <- crossprod(X)
  Qinv <- solve(XTX + A^-2*diag(p))

  for(count in 1:niter)
  {
    # sample [Z|Y,beta,x]
    Z <- numeric(n)
    Xbeta <- X%*%beta_val

    #print(sum(abs(Xbeta) > 5)); mostly zero
    # to make rtrunc behave, make sure values of Xbeta are within
    # the exreme of (-5, 5)
    Xbeta[Xbeta < -5] = -5; Xbeta[Xbeta > 5] = 5

    # Bwoah. rtrunc printing warning message in
    # vectorized version
    for(i in 1:n)
    {
      a = ifelse(Y[i] == 0, -Inf, 0)
      b = ifelse(Y[i] == 0, 0, Inf)
      Z[i] = rtrunc(n = 1, 'norm', a=a, b=b, mean = Xbeta[i], sd = 1)
    }

    # sample [beta|Y, Z, X]
    l <- crossprod(X, Z)
    beta_val <- mvrnorm(mu = Qinv %*% l, Sigma = Qinv)

    betas[count,] = beta_val
  }

  return(as.mcmc(betas))
}
