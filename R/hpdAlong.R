hpdAlong <- function(x, along, prob = 0.95)
{
  # this function is slow...
  # what could be done to make it better is to do HPDinterval(as.mcmc(...))
  #  over more than 1 vector at a time...

  dx = dim(x)
  dim_along = dx[along]

  offidxer = function(dx)
  {
    nn = prod(dx)
    ndim = length(dx)
    ret = array(dim = c(nn, ndim))

    increment = function(v, d)
    {
      v[d] = v[d] + 1
      if(v[d] > dx[d])
      {
        v[d] = 1
        return(increment(v, d-1))
      }
      return(v)
    }

    val = rep(1, ndim)
    for(i in 1:(nn-1))
    {
      ret[i,] = val
      val = increment(val, ndim)
    }
    ret[nn,] = val

    return(ret)
  }

  grab_vector <- function(off_idx)
  {
    idx = numeric(length(dx))
    idx[-along] = off_idx
    indices <- as.list(idx)
    indices[[along]] <- bquote()

    call <- as.call(c(
      list(as.name("["), quote(x)),
      indices))

    return(eval(call))
  }

  set_call <- function(name, idx, val)
  {
    indices = as.list(idx)
    call <- as.call(list(
      as.name("="),
      as.call(c(
        list(as.name("["), as.name(name)),
        indices)),
      val))
    return(call)
  }

  low = array(dim = dx[-along])
  upp = array(dim = dx[-along])

  nn = prod(dx[-along])
  offidxes = offidxer(dx[-along])
  for(idx in 1:nn)
  {
    offidx = offidxes[idx,]
    upper_lower = HPDinterval(as.mcmc(grab_vector(offidx)), prob = prob)

    eval(set_call("low", offidx, upper_lower[1]))
    eval(set_call("upp", offidx, upper_lower[2]))
  }

  return(list(lower = low, upper = upp))
}
