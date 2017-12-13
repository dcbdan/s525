factor_analysis_with_regression <- function(
  Y, # n by p matrix
  X, # n by f matrix
  k, # number of factors
  niter = 1,
  shape_psi = 1/2,    # Can be a scalar or k by 1 vector
                      # 1/2,1/2 gives half-cauchy, cauchy prior
                      # for 1st row, lower elements of lambda
  rate_psi = 1/2,     # Can be a scalar or k by 1 vector
  shape_tau = 1,      # Can be a scalar or p by 1 vector
  rate_tau = 0.2,     # Can be a scalar or p by 1 vector
  coef_multiplier = 10, # beta_lm ~ N(0, coef_multiplier)
  nonzero_structure = NULL) # non zero structure of factor matrix
                         # if not set, factor matrix will
                         # be interpreted to be lower tri
{
  X <- as.matrix(X)

  n <- nrow(Y)
  p <- ncol(Y)
  f <- ncol(X)

  stopifnot(nrow(X) == n)
  if(is.null(nonzero_structure))
  {
    nonzero_structure <- !upper.tri(matrix(nrow = p, ncol = k))
    piv_idx <- 1:k
  }
  else
  {
    stopifnot(nrow(nonzero_structure) == p)
    stopifnot(ncol(nonzero_structure) == k)

    piv_idx = apply(nonzero_structure, 2, function(v){ which(v)[1] })
  }

  num_nonzeros <- apply(nonzero_structure, 1, sum)

  # piv_idx is a k long.
  # The kth element is the first row that has a nonzero value
  # In Ghosh and Dunson, they had a lower triangular matrix
  # (piv_idx = 1:k). But if lamba[l,l] is set to zero, then
  # the transformation used fails. So when a different nonzero
  # structure is specified, a different "pivot index" is used.

  psi <- diag(rgamma(k, shape = shape_psi, rate = rate_psi), nrow = k)
  tau <- diag(rgamma(p, shape = shape_tau, rate = rate_tau), nrow = p)
  eta <- t(matrix(rnorm(n*k, sd=sqrt(1/psi)), nrow = k, ncol = n))
  bbb <- t(matrix(rnorm(k*f, sd=sqrt(coef_multiplier)), nrow = f, ncol = k))
  alpha <- rnorm(p)
  mu <- rnorm(k)

  post_psi <- array(dim = c(niter, k))
  post_tau <- array(dim = c(niter, p))
  post_mu <- array(dim = c(niter, k))
  post_star_alpha <- array(dim = c(niter, p))
  post_star_eta <- array(dim = c(niter, n, k))
  post_star_lambda <- array(dim = c(niter, p, k))
  post_star_B <- array(dim = c(niter, k, f))

  for(iteration in 1:niter)
  {
    # sample lambda
    lambda <- matrix(0, nrow = p, ncol = k)
    for(j in 1:p)
    {
      alpha_j <- alpha[j]
      nnz_j <- num_nonzeros[j]
      tau_j <- tau[j,j]
      Zj <- matrix(eta[, nonzero_structure[j,]], nrow = n, ncol = nnz_j)
      Yj <- Y[,j]

      lambda[j,nonzero_structure[j,]] <-
        s525::rnorm_qinv_l(1,
        Q = tau_j * crossprod(Zj) + diag(nnz_j),
        l = tau_j * colSums(Zj * (Y[,j] - alpha_j)))
    }

    # sample etas
    for(i in 1:n)
    {
      eta[i,] <- s525::rnorm_qinv_l(1,
        Q = t(lambda) %*% tau %*% lambda + psi,
        l = t(lambda) %*% tau %*% (Y[i,] - alpha) +
            psi %*% (mu + bbb %*% X[i,]))
    }

    # sample alpha
    alpha <- s525::rnorm_qinv_l(1,
      Q = diag(p) + n*tau,
      l = tau %*% colSums(Y - tcrossprod(eta, lambda)))

    # sample mu
    mu <- s525::rnorm_qinv_l(1,
      Q = diag(k) + n*psi,
      l = psi %*% colSums(eta - tcrossprod(X, bbb)))

    # sample bbb
    for(l in 1:k)
    {
      bbb[l,] <- s525::rnorm_qinv_l(1,
        Q = psi[l,l]*crossprod(X) + 1/coef_multiplier*diag(f),
        l = psi[l,l]*colSums(X * (eta[,l] - mu[l])))
    }

    # sample tau
    for(j in 1:p)
    {
      tau[j,j] <- rgamma(1,
        shape = n/2 + shape_tau,
        rate = rate_tau + 1/2*sum((Y[,j] - alpha[j] - lambda[j,] %*% t(eta))^2))
    }

    # sample psi
    for(l in 1:k)
    {
      psi[l,l] <- rgamma(1,
        shape = n/2 + shape_psi,
        rate = rate_psi + 1/2*sum((eta[,l] - mu[l] - bbb[l,] %*% t(X))^2))
    }

    post_psi[iteration,] = diag(psi)
    post_tau[iteration,] = diag(tau)          # to be returned
    post_mu[iteration,] = mu
    post_star_alpha[iteration,] = alpha       # to be transformed and returned
    post_star_lambda[iteration,,] = lambda    # to be transformed and returned
    post_star_eta[iteration,,] = eta          # to be transformed and returned
    post_star_B[iteration,,] = bbb            # to be transformed and returned
  }

  # transform the following
  post_alpha <- array(dim = c(niter, p))
  post_lambda <- array(dim = c(niter, p, k))
  post_eta <- array(dim = c(niter, n, k))
  post_B <- array(dim = c(niter, k, f))

  # this function uses n, p, k, f, piv_idx from environment
  transform_stuff <- function(
    alpha_star,
    lambda_star,
    eta_star,
    B_star,
    mu_star,
    psi)
  {
    lambda_star <- as.matrix(lambda_star, nrow = p, ncol = k)
    eta_star <- as.matrix(eta_star, nrow = n, ncol = k)
    B_star <- as.matrix(B_star, nrow = k, ncol = f)

    ret_alpha <- numeric(p)
    ret_lambda <- matrix(nrow = p, ncol = k)
    ret_eta <- matrix(nrow = n, ncol = k)
    ret_B <- matrix(nrow = k, ncol = f)

    ret_alpha <- alpha_star + lambda_star %*% mu_star
    ret_B <- B_star * sqrt(psi)

    for(l in 1:k)
    {
      ret_eta[,l]    =
        sign(lambda_star[piv_idx[l],l])*(eta_star[,l] - mu_star[l])*sqrt(psi[l])
      ret_lambda[,l] =
        sign(lambda_star[piv_idx[l],l])*lambda_star[,l]*sqrt(1/psi[l])
    }

    return(list(ret_alpha, ret_lambda, ret_eta, ret_B))
  }

  for(iteration in 1:niter)
  {
    tmp <- transform_stuff(
      post_star_alpha[iteration,],
      post_star_lambda[iteration,,],
      post_star_eta[iteration,,],
      post_star_B[iteration,,],
      post_mu[iteration,],
      post_psi[iteration,])

    post_alpha[iteration,] = tmp[[1]]
    post_lambda[iteration,,] = tmp[[2]]
    post_eta[iteration,,] = tmp[[3]]
    post_B[iteration,,] = tmp[[4]]
  }

  ret <- list(post_alpha, post_lambda, post_tau, post_eta, post_B)
  names(ret) <- c("alpha", "lambda", "tau", "eta", "B")
  return(ret)
}
