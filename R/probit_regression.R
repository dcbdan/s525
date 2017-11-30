probit_regression <- function(
  Y,
  X,
  niter,
  A = 3,
  init = NULL)
{
  Y0 = Y==0; Y1 = Y==1; n0 <- sum(Y0); n1 <- sum(Y1)
  n <- nrow(X); p <- ncol(X);

  betas <- array(0, c(niter, p))
  if(is.null(init))
  {
    beta_val <- numeric(p) #rnorm(p, mean = 0, sd = A)
  }
  else
  {
    beta_val <- init
  }

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
    #Xbeta[Xbeta < -5] = -5; Xbeta[Xbeta > 5] = 5

    # Bwoah. rtrunc printing warning message in
    # vectorized version
    ### for(i in 1:n)
    ### {
    ###   a = ifelse(Y[i] == 0, -Inf, 0)
    ###   b = ifelse(Y[i] == 0, 0, Inf)
    ###   Z[i] = rtrunc(n = 1, 'norm', a=a, b=b, mean = Xbeta[i], sd = 1)
    ### }
    # Using rtruncnorm instead
    Z[Y0] = rtruncnorm(n = n0, a=-Inf, b = 0, mean = Xbeta[Y0], sd = 1)
    Z[Y1] = rtruncnorm(n = n1, a=0,  b = Inf, mean = Xbeta[Y1], sd = 1)


    # sample [beta|Y, Z, X]
    l <- crossprod(X, Z)
    beta_val <- mvrnorm(mu = Qinv %*% l, Sigma = Qinv)

    betas[count,] = beta_val
  }

  colnames(betas) <- colnames(X)

  return(as.mcmc(betas))
}
