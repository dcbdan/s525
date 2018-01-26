sumAlong <- function(x, along)
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

  ret = grab_submat(1)
  if(dim_along > 1)
  {
    for(idx in 2:dim_along)
    {
      ret = ret + grab_submat(idx)
    }
  }

  return(ret)
}

