applyAlong <- function(x, f, along = 1)
{
  dx = dim(x)
  dim_along = dx[along]

  grab_submat <- function(idx)
  {
    indices <- rep(list(bquote()), length(dx))
    indices[[along]] <- idx

    call <- as.call(c(
      list(as.name("["), quote(x)),
      indices))

    return(eval(call))
  }

  # return the call equivalent to
  #   oldval[[nm]][idx,...ndim...] = newval
  set_submat_call_helper <- function(val, idx, ndim, newval)
  {
    indices <- rep(list(bquote()), ndim)
    indices[[1]] <- idx

    return(
      as.call(list(
        as.name("="),
        as.call(c(list(as.name("["), as.name(val)), indices)),
        as.name(newval))))
  }


  first_val = f(grab_submat(1))
  if(is.vector(first_val))
  {
    dim_vals = length(first_val)
  }
  else
  {
    dim_vals = dim(first_val)
  }

  ndim_vals = length(dim_vals)
  ret = array(dim = c(dim_along, dim_vals))

  the_first_call = set_submat_call_helper("ret", 1, 1+ndim_vals, "first_val")
  eval(the_first_call)

  if(dim_along > 1)
  {
    for(idx in 2:dim_along)
    {
      newval = f(grab_submat(idx))

      # Return the following call
      #ret[idx,,,,,,,,,,,,,,,] = newval
      the_call = set_submat_call_helper("ret", idx, 1+ndim_vals, "newval")

      # evaluate the call
      eval(the_call)
    }
  }

  return(ret)
}
