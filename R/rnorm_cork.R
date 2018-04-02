rnorm_cork = function(
  n,
  mean.vec, #   Sample X ~ normal with mean.vec, cov.mat subject
  cov.mat,  #   to the constraint that M^tX <= m, A^tX = a
  ineq.mat = NULL, #M
  ineq.vec = NULL, #m
  eq.mat = NULL,   #A
  eq.vec = NULL)   #a
{
  #Schmidt, Mikkel. "Linearly constrained bayesian matrix factorization for blind source separation." Advances in neural information processing systems. 2009.
  nd = length(mean.vec)

  # mean.vec is nd x 1
  M = ineq.mat # nMconst x nd
  m = ineq.vec # nMconst x 1
  A = eq.mat   # nAconst x nd
  a = eq.vec   # nAconst x 1

  mu  = mean.vec # nd x 1
  sig = cov.mat  # nd x nd

  if(is.null(A) || is.null(a))
  {
    muy = mu             # nd x 1
    sigychol = chol(sig) # nd x nd
  }
  else
  {
    ############################################################################
    # calculate xo, Tt, Tto
    nc = ncol(A)
    svdinfo = svd(A, nu = nd, nv = nc) # nv does not need to be set to nc
                                       # but VV will need to be square

    UU = svdinfo$u
    VV = svdinfo$v

    len_d = length(svdinfo$d)
    stopifnot(len_d != nd) # make sure there can be more than one solution...

    DDi = array(0, dim = c(nrow(UU), nrow(VV)))
    diag(DDi)[1:len_d] = 1/svdinfo$d

    xo = UU %*% DDi %*% t(VV) %*% a
    Tt  = t(UU[1:nd,1:len_d])
    Tto = t(UU[1:nd,(len_d+1):nd])
    ############################################################################

    ############################################################################
    # calculate muy, sigychol
    tStinv = chol2inv(chol(Tt %*% sig %*% t(Tt)))
    lambda = Tto %*% (diag(nd) - sig %*% t(Tt) %*% tStinv %*% Tt)

    muy = lambda %*% (mu - xo)
    sigychol = chol(lambda %*% sig %*% t(Tto))
    ############################################################################
  }

  ##############################################################################
  # sample z
  nz = nrow(sigychol)
  if(is.null(M) || is.null(m))
  {
    z = array(rnorm(nz*n), dim = c(nz, n))
  }
  else
  {
    ############################################################################
    # calculate Mz, mz
    if(is.null(A) || is.null(a))
    {
      #Mz = crossprod(sigychol, M)
      Mz = sigychol %*% M
      mz = m - crossprod(M, muy)
    }
    else
    {
      #Mz = crossprod(sigychol, Tto) %*% M
      Mz = sigychol %*% Tto %*% M
      mz = m - crossprod(M, xo) - crossprod(M, crossprod(Tto, muy))
    }
    ############################################################################

    if(n < 50)
    {
      burn = 50 # arbitrary # TODO
    }
    else
    {
      burn = n
    }

    z = t(rstnorm_ineq_cork(n, Mz, mz, burn = burn))
  }
  ##############################################################################

  if(is.null(A) || is.null(a))
  {
    muy = as.vector(muy)
#    return(t(sigychol %*% z + muy)) # z is nd x n, muy is nd x 1
    return(t(crossprod(sigychol, z) + muy)) # x is nd x n, muy is nd x 1
  }
  else
  {
    muy = as.vector(muy)
    xo = as.vector(xo)
#    return(t(crossprod(Tto, sigychol %*% z + muy) + xo)) # z is nz x n
    return(t(crossprod(Tto, crossprod(sigychol, z) + muy) + xo))
  }
}


