probit_horseshoe_regression <- function(
  Y,
  X,
  niter,
  init = NULL)
{
  Y0 <- Y==0; Y1 <- Y==1; n0 <- sum(Y0); n1 <- sum(Y1)
  n <- nrow(X); p <- ncol(X);
  stopifnot(p > 1)
  stopifnot(X[,1] == rep(1, n))

  betas <- array(0, c(niter, p))
  if(is.null(init))
  {
    beta_val <- numeric(p)
  }
  else
  {
    beta_val <- init
  }

  XTX <- crossprod(X)

  gamma_val <- rgamma(p-1, shape = 1/2, rate = 1/25)
  eta_val <- rgamma(p-1, shape = 1/2, rate = gamma_val)

  for(count in 1:niter)
  {
    # sample [Z|Y,beta,x], Z is Y*
    Z <- numeric(n)
    Xbeta <- X%*%beta_val
    Z[Y0] = rtruncnorm(n = n0, a=-Inf, b = 0, mean = Xbeta[Y0], sd = 1)
    Z[Y1] = rtruncnorm(n = n1, a=0,  b = Inf, mean = Xbeta[Y1], sd = 1)

    # sample [tau_A|gamma]
    tau_A_val = rtrunc(
      n = 1,
      'gamma',
      a = 1/100,
      b = Inf,
      shape = (p-2)/2,
      rate = sum(gamma_val))

    # sample [gamma|eta, tau]
    gamma_val = rexp(p-1, rate = eta_val + tau_A_val)

    # sample [eta|beta, gamma]
    eta_val = rexp(p-1, rate = beta_val[2:p]^2/2 + gamma_val)

    # sample [beta|Y, Z, X, eta, gamma]
    Qinv <- solve(XTX + diag(c(0, eta_val)))
    l = crossprod(X, Z)
    beta_val <- mvrnorm(1, mu = Qinv %*% l, Sigma = Qinv)

    betas[count,] = beta_val
  }

  colnames(betas) <- colnames(X)

  return(as.mcmc(betas))
}
