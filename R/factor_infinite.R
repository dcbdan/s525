factor_infinite <- function(
  Y, # n by p matrix
  k_star, # number of factors to cut off at
  niter = 1000,
  a_sig = 1,
  b_sig = 1,
  rho = 3,
  a1 = 1,
  a2 = 3,
  prt = FALSE)
{
  n <- nrow(Y)
  p <- ncol(Y)
  k = k_star

  stopifnot(k != 1)
  stopifnot(p > k_star)
  stopifnot(a2 > 2)

  # set initial values
  eta <- array(rnorm(n*k), dim = c(n, k))
  del <- rgamma(k, c(a1, rep(a2, k-1)), 1)
  tau <- sapply(1:k, function(h){ prod(del[1:h]) })
  s2g <- 1/rgamma(p, a_sig, b_sig) # s2g is cool kid slang for \sigma^2
  phi <- array(rgamma(p*k, rho/2, rho/2), dim = c(p, k))

  # will be set in the loop
  lam <- array(dim = c(p, k))

  # return values
  post_lam <- array(dim = c(niter, p, k))
  post_s2g <- array(dim = c(niter, p))
  post_eta <- array(dim = c(niter, n, k))

  # run the gibbs sampler
  for(iteration in 1:niter)
  {
    # lam
    for(j in 1:p)
    {
      Q = diag(phi[j,]*tau) + 1/s2g[j]*crossprod(eta)
      l = 1/s2g[j]*crossprod(eta, Y[,j])
      lam[j,] = rnorm_qinv_l(1, Q, l)
    }

    # s2g
    s2g = 1/rgamma(p,
      shape = a_sig + n/2,
      rate = b_sig + 1/2*colSums((Y - eta %*% t(lam))^2))
    # loop version
    #for(j in 1:p)
    #{
    #  s2g[j] = 1/rgamma(1,
    #    shape = a_sig + n/2,
    #    rate = b_sig + 1/2*sum((Y[,j] - eta %*% lam[j,])^2))
    #}

    # eta
    for(i in 1:n)
    {
      lam_sig = crossprod(lam, diag(1/s2g))
      Q = diag(k) + lam_sig %*% lam
      l = lam_sig %*% Y[i,]
      eta[i,] = rnorm_qinv_l(1, Q, l)
    }

    # phi
    for(j in 1:p)
    {
      phi[j,] = rgamma(k,
        shape = 1/2*(rho + 1),
        rate =  1/2*(rho + tau*(lam[j,]^2)))
    }

    # tau
    delta_helper <- function(h)
    {
      tauify = function(l){ prod(del[(1:l)[-h]]) }
      ret = 0
      for(l in h:k)
      {
        ret = ret + tauify(l)*sum(phi[,l]*(lam[,l]^2))
      }
      return(ret)
    }
    del = rgamma(k,
      shape = c(a1, rep(a2, k-1)) + p/2*(k + 1 - 1:k),
      rate = 1 + 1/2*sapply(1:k, delta_helper))
    tau = sapply(1:k, function(h){ prod(del[1:h]) })

    # store lam and s2g
    post_lam[iteration,,] = lam
    post_s2g[iteration,] = s2g
    post_eta[iteration,,] = eta

    if(prt)
    {
      print(iteration)
    }
  }

  return(list(
    sigma2 = post_s2g,
    lambda = post_lam,
    eta = post_eta))
}
