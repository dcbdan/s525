horseshoe_regression <- function(
  Y,
  X,
  niter = 10000)
{
  chain_length <- niter
  XTX <- crossprod(X)
  XTY <- crossprod(X, Y)

  n <- length(Y)
  p <- ncol(X)

  tau <- runif(1)
  tau_a <- (1/runif(1,max = 10))^2
  gammas <- rgamma(p, shape = 1/2, rate = tau_a)
  etas <- rgamma(p, shape = 1/2, rate = gammas)
  betas <- mvrnorm(1, rep(0, p), Sigma = 1/tau*diag(1/etas))

  # TODO why have all these?
  post_betas <- array(0, c(chain_length, p))
  post_etas <- array(0, c(chain_length, p))
  post_tau <- array(0, c(chain_length,1))
  post_tau_a <- array(0, c(chain_length, 1))

  for(count in 1:chain_length)
  {
    etas = rexp(p, tau*betas^2/2 + gammas)

    Lam <- diag(etas)
    #Qi <- solve(XTX + Lam)
    #mu <- Qi %*% XTY
    #Sig <- 1/tau * Qi
    #betas = mvrnorm(1, mu = mu, Sigma = Sig)
    betas = sqrt(1/tau) *
      rpost_regression_coef(
        X,
        Lam,
        Y,
        u = rnorm(p, sd = sqrt(Lam)))

    gammas = rexp(p, etas + tau_a)

    # tau_a can't be less than 1/100 since p(tau_a) is zero
    # outside of (1/100, infinity)
    tau_a = rtrunc(
      n = 1,
      'gamma',
      a = 1/100,
      b = Inf,
      shape = (p-1)/2,
      rate = sum(gammas))

    ols <- crossprod(Y - X%*%betas)
    rrr <- t(betas) %*% Lam %*% betas
    tau = rgamma(1, shape = (n+p)/2, rate = (ols + rrr)/2)

    post_betas[count, ] = betas
    post_etas[count, ] = etas
    post_tau[count] = tau
    post_tau_a[count] = tau_a
  }

  post_sigma <- 1/sqrt(post_tau)
  post_lambdas <- 1/sqrt(post_etas)
  ret <- list(
    as.mcmc(post_betas),
    as.mcmc(post_lambdas),
    as.mcmc(post_sigma))
  names(ret) <- c("beta", "lambda", "sigma")
  colnames(ret$beta) <- paste("beta_",colnames(X), sep="")
  colnames(ret$lambda) <- paste("lambda_", colnames(X), sep="")
  colnames(ret$sigma) <- "sigma"

  return(ret)
}
