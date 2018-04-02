rstnorm_ineq_cork = function(
  n,
  ineq.mat, # sample z ~ N(0, I) constrained so that t(ineq_mat)z <= ineq_vec
  ineq.vec,
  thin = 1, # take every thinth value
  burn = 0, # ignore the first burn
  init = NULL)
{
  M = ineq.mat # renaming variables. t(M) %*% z <= m
  m = ineq.vec

  nd = nrow(M) # number of dimensions of z

  # do one step of the gibbs sampler on p
  increment = function(p, rep = 1)
  {
    if(rep <= 0)
    {
      return(p)
    }

    for(idx in 1:nd)
    {
      dvec = M[idx,]

      # calculate nvec ...
      if(nd == 1)
      {
        nvec = m
      }
      else if(nd == 2)
      {
        nvec = m - M[-idx,]*p[-idx]
      }
      else
      {
        nvec = m - crossprod(M[-idx,], p[-idx])
      }

      dless = dvec < 0
      dgrea = dvec > 0

      # assuming rtruncnorm is effecient
      a = max(-Inf, nvec[dless]/dvec[dless])
      b = min( Inf, nvec[dgrea]/dvec[dgrea])

      # a==b might happen if it is far from converging
      # all.equal(a, b) == T is not redundant.........all.equal does
      # not return FALSE, only TRUE...............................
      if(is.finite(a) && is.finite(b) && (all.equal(a, b) == T))
      {
        p[idx] = a
      }
      else
      {
        p[idx] = rtruncnorm(1, a, b)
      }
    }

    return(increment(p, rep - 1))
  }

  # we need an initial starting point that is an acceptable value.
  # If not provided, use z such that t(M)z = m.
  if(is.null(init) || crossprod(M, init) > m)
  {
    # TODO
    # computationally, this should be equivalent to doing an svd...
    Mt_pseudoinv = ginv(t(M))

    z = Mt_pseudoinv %*% m
  }
  else
  {
    z = init
  }

  z = increment(z, burn - (thin - 1))

  ret = array(dim = c(n, nd))
  for(iter in 1:n)
  {
    z = increment(z, thin)
    ret[iter,] = z
  }

  return(ret)
}


