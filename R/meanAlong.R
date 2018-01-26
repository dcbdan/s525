meanAlong <- function(x, along)
{
  dx = dim(x)
  dim_along = dx[along]
  sumAlong(x, along) / dim_along
}

