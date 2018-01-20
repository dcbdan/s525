plackett_luce <- function(
  rank_matrix, # a n x p matrix of integers ranging from 1,...,k
      # where k = max(rank_matrix, na.rm=TRUE).
      # values in dd can be NA but must have atleast 2 available values per row
      # note: same integer can appear multiple times in each row of rank matrix
  shape = 1,  # prior on lambda
  rate = 1,   # prior on lambda
  niter = 1000)
{
  n <- nrow(rank_matrix)
  p <- ncol(rank_matrix)
  k <- max(rank_matrix, na.rm = TRUE)

  ########################
  # preliminary stuffs

  # set dd and pp
  dd <- array(NA, dim = c(n, p)) # dd is rank_matrix with NA's at back
  pp <- numeric(p) # pp contains number of rankings in each row
  for(i in 1:n)
  {
    d <- rank_matrix[i, !is.na(rank_matrix[i,])]
    pp[i] <- length(d)

    # arrange dd to make sure NAs are at the front
    dd[i,1:pp[i]] = d
  }

  # set ww, ww contains number of times each ranked item k is
  # equal to dd[i,j] for i = 1,...,n, j = 1,...,(pp[i]-1)
  ww <- numeric(k)
  # ugly nested for loops but it is readable and prob correct is greater
  for(i in 1:n)
  {
    for(j in 1:(pp[i]-1))
    {
      l = dd[i,j]
      ww[l] = ww[l] + 1
    }
  }

  # set delt, delt[i,j,l] is equal to sum_{m=j}^{pp[i]}{dd[i,j] == l}
  delt <- array(0, dim = c(n, p-1, k))
  # ugly nested for loops but it is readable and prob correct is greater
  for(i in 1:n)
  {
    for(j in 1:(pp[i] - 1))
    {
      for(l in 1:k)
      {
        delt[i,j,l] = sum(dd[i,j:pp[i]] == l)
      }
    }
  }

  ########################

  post_lambda <- array(dim = c(niter, k))
  Z <- array(0, dim = c(n, p - 1))

  lambda = rgamma(k, shape = shape, rate = rate)

  for(iter in 1:niter)
  {
    # sample Z
    for(i in 1:n)
    {
      for(j in 1:(pp[i]-1))
      {
        Z[i,j] = rexp(1, rate = sum(lambda[dd[i,j:pp[i]]]))
      }
    }

    # sample lambda
    lambda = rgamma(
      k,
      shape = shape + ww,
      rate = rate + sapply(1:k, function(l){ sum(delt[,,l] * Z) }))

    post_lambda[iter,] = k * lambda / sum(lambda)
  }

  colnames(post_lambda) = paste("lambda", int_string(k), sep = "")
  return(post_lambda)
}


