ridge_regression <- function(
  Y,
  X,
  niter = 10000)
{
  S = niter
  y = Y
  p = length(X[1,])

  mod_ols = lm(y ~ X - 1)

  # OLS initialization:
  beta = mod_ols$coefficients
  sigma = sd(y - X%*%beta)
  tau = 1/sigma^2

  A = 1000 # Uniform prior on SD in ridge parameter
  sigma_beta = sd(beta)
  #sigma_beta = 1
  eta_beta = 1/sigma_beta^2

  # Recurring terms:
  XtX = crossprod(X)
  Xty = crossprod(X, y)

  post_beta = array(0, c(S, p))
  post_sigma = post_sigma_beta = array(0, c(S, 1))#numeric(S)
  for(s in 1:S)
  {
    # Block 1: sample sigma (via tau)
      # Note: no longer marginalizing over beta (although we could...)
    tau = rgamma(n = 1,
                 shape = .01 + n/2,
                 rate = .01 + sum((y - X%*%beta)^2)/2)
    sigma = 1/sqrt(tau)

    # Block 2: sample beta
    Q_beta = tau*XtX + diag(eta_beta, p)
    ell_beta = tau*Xty

    # Cholesky, then forward/backsolve:
    ch_Q = chol(Q_beta)
    beta = backsolve(ch_Q,
                     forwardsolve(t(ch_Q), ell_beta) +
                       rnorm(p))

    # Block 3: sample sigma_beta ~ Uniform(0, A)
    eta_beta = rtrunc(n = 1,
                       'gamma',   # Family of distribution
                       a = 1/A^2, # Lower interval
                       b = Inf,   # Upper interval
                       shape = p/2 - 1/2,
                       rate =  sum(beta^2)/2)
    sigma_beta = 1/sqrt(eta_beta)

    post_beta[s,] = beta
    post_sigma[s] = sigma
    post_sigma_beta[s] = sigma_beta
  }

  ret <- list(
    as.mcmc(post_beta),
    as.mcmc(post_sigma),
    as.mcmc(post_sigma_beta))
  names(ret) <- c("beta", "sigma", "sigma_beta")
  colnames(ret$beta) <- paste("beta_",colnames(X), sep="")
  colnames(ret$sigma) <- "sigma"

  return(ret)
}

