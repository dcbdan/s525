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

  perturb = function(p, countdown = 2)
  {
    if(countdown == 0)
      return(p)

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

#      print("PERTURB")
#      print(c(a,b))

      # a==b might happen if it is far from converging
      # all.equal(a, b) == T is not redundant.........all.equal does
      # not return FALSE, only TRUE...............................
      if(is.finite(a) && is.finite(b) && (all.equal(a, b) == T))
      {
        p[idx] = a
      }
      else if(is.finite(a) && is.finite(b))
      {
        if (a > b)
          p[idx] = mean(c(a,b))
        else
          p[idx] = runif(1, a, b)
      }
      else
      {
        p[idx] = rtruncnorm(1, a, b, sd = 10) # more or less uniform
      }
#      print(p[idx])
    }

    stopifnot(all(!is.nan(p)))
    stopifnot(all(!is.na(p)))

    if(all(crossprod(M, p) <= m))
      return(p)
    else
      return(perturb(p, countdown - 1))
  }

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

#      print("INCREMENT")
#      print(c(a,b))

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

    stopifnot(all(!is.nan(p)))
    stopifnot(all(!is.na(p)))

    return(increment(p, rep - 1))
  }

  # we need an initial starting point that is an acceptable value.
  # If not provided, use z such that t(M)z = m.
  if(is.null(init) || any(crossprod(M, init) > m))
  {
    # TODO
    # computationally, this should be equivalent to doing an svd...
    Mt_pseudoinv = ginv(t(M))


#    print(M)

    z = Mt_pseudoinv %*% m
#    print(range(crossprod(M, z) - m))
#    print(t(M) %*% z <= m)
#    z = perturb(z) # the pseudo inverse solution is ok but it puts the
                   # initial value on the edge of the equality and
                   # that can cause truncnorm to freak out.
                   # perturb does the same thing as increment except
                   # instead of rtruncnorm(1, a, b) it does
                   # runif(1, a, b) (for idxs where a and b are finite)
  }
  else
  {
    z = init
  }

  stopifnot(all(crossprod(M, z) <= m))

  z = increment(z, burn - (thin - 1))

  ret = array(dim = c(n, nd))
  for(iter in 1:n)
  {
    z = increment(z, thin)
    ret[iter,] = z
  }

  return(ret)
}


