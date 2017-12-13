mat_apply <- function(vec, fun)
{
  n <- length(vec)

  if(n == 0)
  {
    return(NULL)
  }

  if(n == 1)
  {
    return(fun(vec[1]))
  }

  val1 <- fun(vec[1]) # fun must return array; not vector
  dim_val <- dim(val1)
  neach <- prod(dim_val)
  ret <- array(dim=c(n, neach))
  ret[1,] = val1
  for(iter in 2:n)
  {
    ret[iter,] = as.vector(fun(vec[iter]))
  }
  return(array(ret, dim=c(n, dim_val)))
}

