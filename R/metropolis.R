metropolis <- function(
  rproposal, # rproposal(previous_val) generates a value from
             # the proposal distribution
  prob, # prob(val_proposed)/prob(val_t_minus_1) helps form the acceptance prob
  niter, # number of iterations to perform
  init, # vector of initial values
  log_prob = FALSE) # whether or not prob specifies log probabilities
{
  # helper function to
  # calculate r from either the log probabilities or
  # the probabilities
  calc_r <- function(proposed_val, prev_val)
  {
    if(log_prob)
    {
      return(exp(prob(proposed_val) - prob(prev_val)))
    }
    else
    {
      return(prob(proposed_val)/prob(prev_val))
    }
  }

  d <- length(init) # the number of dimensions
  prev_val <- init

  ret <- array(0, c(niter, d))

  number_accepted <- 0
  for(i in 1:niter)
  {
    proposed_val <- rproposal(prev_val)

    r <- calc_r(proposed_val, prev_val)
    u <- runif(1)

    if(u < r)
    {
      ret[i,] = proposed_val
      number_accepted = number_accepted + 1
    }
    else
    {
      ret[i,] = prev_val
    }

    prev_val = ret[i,]
  }
  #print(c(number_accepted, effectiveSize(ret))/niter)

  return(as.mcmc(ret))
}

