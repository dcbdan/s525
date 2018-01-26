gibbs <- function(
  niter,                # number of iterations to run sampler
  init,                 # initial vars to start. A list with names
  hypers,               # constant hyper parameter values, a list
  known_data,           # data that can be considered constant, a named list also
  conditional_samplers, # conditional samplers.
                        #  A list of functions with names corresponding to
                        #  init. known_data is fed into each function as well.
  iter_argname = "iter", # iteration argname for samplers to accept
  ignore = c()) # vector of parameters to not store post values
{
  # We need init and conditional_samplers to have the same
  # names. They need not be in the same order, though
  var_names = names(init)
  stopifnot(sort(var_names) == sort(names(conditional_samplers)))

  # remove vectors and replace with arrays
  vars = lapply(init, as.array)

  # allocate output
  to_keep_vars = sapply(var_names, function(vnm){ all(vnm != ignore) })
  ret = lapply(vars[to_keep_vars], function(v){ array(dim = c(niter, dim(v))) })

  # for use in grab_submat
  ndims = lapply(ret, function(m){ length(dim(m)) })

  # get list of argnames for each conditional sampler
  argnames_possible = c(var_names, names(known_data), iter_argname, names(hypers))
  get_f_vars = function(var_name)
  {
    f_args = names(formals(conditional_samplers[[var_name]]))
    if(any(f_args == "..."))
    {
      return(argnames_possible)
    }
    else
    {
      ret = c()
      for(arg in argnames_possible)
      {
        if(any(arg == f_args))
        {
          ret = c(ret, arg)
        }
      }
      return(ret)
    }
  }
  f_args = lapply(var_names, get_f_vars)
  names(f_args) = var_names

  # return the call equivalent to
  #   oldval[[nm]][idx,...ndim...] = newval[[nm]]
  submat_call_helper <- function(nm, oldval, idx, ndim, newval)
  {
    indices <- rep(list(bquote()), ndim)
    indices[[1]] <- idx

    old_ = as.call(list(as.name('[['), as.name(oldval), nm))
    new_ = as.call(list(as.name('[['), as.name(newval), nm))

    return(
      as.call(list(
        as.name("="),
        as.call(c(
          list(as.name("["), old_), indices)),
        new_)))
  }

  # Do the gibbs sampling step!
  for(iter in 1:niter)
  {
    # for each variable,
    #   sample conditional on all the other variables,
    #   store in ret
    for(var_name in var_names)
    {
      #vars[var_name] = contional_samplers[var_name](vars, known_data, iter)
      argnames = (c(vars, known_data, hypers, list(iter=iter)))[f_args[[var_name]]]
      vars[[var_name]] = do.call(
        conditional_samplers[[var_name]],
        argnames)

      # the following line returns the following r code
      # to be evaulated in this context:
      #   ret[[var_name]][iter,...,] = vars[[var_name]]
      # The trick is that ret[var_name] can have arbitrary number of
      # dimensions
      if(to_keep_vars[var_name])
      {
        set_ret_var_name_to = submat_call_helper(
          var_name,
          "ret",
          iter,
          ndims[[var_name]],
          "vars")

        # evaluate in this context
        eval(set_ret_var_name_to)
      }
    }
  }

  return(lapply(ret, as.mcmc))
}

