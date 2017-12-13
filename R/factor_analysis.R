factor_analysis <- function(
  Y, # n by p matrix
  k, # number of factors
  niter = 1000,
  shape_psi = 1/2,    # Can be a scalar or k by 1 vector
                      # 1/2,1/2 gives half-cauchy, cauchy prior
                      # for 1st column, lower elements of lambda
  rate_psi = 1/2,     # Can be a scalar or k by 1 vector
  shape_sigma2 = 1, # Can be a scalar or p by 1 vector
  rate_sigma2 = 0.2,  # Can be a scalar or p by 1 vector
  nonzero_structure = NULL) # non zero structure of factor matrix
                         # if not set, factor matrix will
                         # be interpreted to be lower tri
{
  n <- nrow(Y)
  p <- ncol(Y)
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

  phis <- rgamma(k, shape = shape_psi, rate = rate_psi)
  taus <- rgamma(p, shape = shape_sigma2, rate = rate_sigma2)
  etas <- t(matrix(rnorm(n*k, sd=sqrt(1/phis)), nrow = k, ncol = n))

  post_psi <- array(dim = c(niter, k))
  post_sigma <- array(dim = c(niter, p))
  post_star_etas <- array(dim = c(niter, n, k))
  post_star_lambdas <- array(dim = c(niter, p, k))
  #post_star_lambdas <- lapply(
  #  num_nonzeros,
  #  function(nnz){ matrix(nrow = niter, ncol = nnz) })

  for(iteration in 1:niter)
  {
    # sample lambdas
    lambdas <- matrix(0, nrow = p, ncol = k)
    for(j in 1:p)
    {
      nnz_j <- num_nonzeros[j]
      tau_j <- taus[j]
      Zj <- etas[, nonzero_structure[j,]]
      Yj <- Y[,j]

      lambdas[j,nonzero_structure[j,]] <-
        s525::rpost_regression_coef(
          X = sqrt(tau_j)*Zj,
          D = diag(nnz_j),
          alpha = sqrt(tau_j)*Yj,
          u = rnorm(nnz_j))
    }

    # sample etas
    for(i in 1:n)
    {
      D_eta <-
        if(k == 1)
          1/phis
        else
          diag(1/phis)
      etas[i,] <- s525::rpost_regression_coef(
        X = diag(sqrt(taus)) %*% lambdas,
        D = D_eta,
        alpha = diag(sqrt(taus)) %*% Y[i,],
        u = rnorm(k, sd = sqrt(1/phis)))
    }

    # sample phis = inverse of psis
    phis <- rgamma(
      k,
      shape = shape_psi + n/2,
      rate = rate_psi + 1/2*apply(etas^2, 2, sum))

    # sample taus = sigmas^(-2)
    tau_rate_tmp <- numeric(p)
    for(i in 1:n)
    {
      tau_rate_tmp <- tau_rate_tmp + crossprod(y[i,] - lambdas %*% etas[i,])
    }
    taus <- rgamma(
      p,
      shape = shape_sigma2 + n/2,
      rate = rate_sigma2 + 1/2*tau_rate_tmp)

    post_psi[iteration,] = 1/phis
    post_sigma[iteration,] = 1/sqrt(taus)
    post_star_etas[iteration,,] = etas
    post_star_lambdas[iteration,,] = lambdas
    #for(j in 1:p)
    #{
    #  post_star_lambdas[j] = lambdas[j, nonzero_structure[j,]]
    #}
  }

  # transform post_star_etas and post_star_lambdas
  post_etas <- array(dim = c(niter, n, k))
  post_lambdas <- array(dim = c(niter, p, k))

  # this function uses n, p, k, piv_idx from environment
  transform_eta_lambda <- function(eta_star, lambda_star, phi)
  {
    ret_eta <- matrix(nrow = n, ncol = k)
    ret_lambda <- matrix(nrow = p, ncol = k)

    if(k == 1)
    {
      lambda_star <- as.matrix(lambda_star)
      eta_star <- as.matrix(eta_star)
    }

    for(l in 1:k)
    {
      ret_eta[,l]    =
        sign(lambda_star[piv_idx[l],l])*eta_star[,l]*sqrt(phi[l])
      ret_lambda[,l] =
        sign(lambda_star[piv_idx[l],l])*lambda_star[,l]*sqrt(1/phi[l])
    }
    return(list(ret_eta, ret_lambda))
  }

  for(iteration in 1:niter)
  {
    tmp <- transform_eta_lambda(
      post_star_etas[iteration,,],
      post_star_lambdas[iteration,,],
      1/post_psi[iteration,])
    post_etas[iteration,,] = tmp[[1]]
    post_lambdas[iteration,,] = tmp[[2]]
  }

  ret <- list(post_sigma, post_etas, post_lambdas)
  names(ret) <- c("sigma", "eta", "lambda")
  return(ret)
}

