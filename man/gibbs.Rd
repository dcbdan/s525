\name{gibbs}
\alias{gibbs}
\title{
  Gibbs Sampler Skeleton
}
\description{
  This function provides the common skeleton used in all gibbs sampler:
  set initial values and then for as many iterations as needed,
  sample conditional posteriors of each of the parameters
}
\usage{
gibbs <- function(
  niter,
  init,
  hypers,
  known_data,
  conditional_samplers,
  iter_argname = "iter",
  ignore = c())
}
\arguments{
  \item{niter}{number of iterations to run sampler}
  \item{init}{A list of initial values to start the sampler. The list must be
    named and the names must correspond to the names in conditional_samplers.}
  \item{hypers}{A list of hyper values. This is passed into the conditional samplers
    as needed.}
  \item{known_data}{A list of data that is passed into the conditional samplers.}
  \item{conditional_samplers}{A list of functions. Each function is a sampler
    for each of the parameters. The list must be named and the names must correspond
    to the names in init.}
  \item{iter_argname}{The conditional samplers may choose to take iter_argname as
    a function argument. For example, if the conditional sampler is only to
    perform a sample every so many iterations.}
  \item{ignore}{A vector of parameter names to not store and return.}
  \item{asmcmc}{Whether or not to make the return value a mcmc object.}
}
\details{
  This function returns a list of arrays where each element is the set of samples
  generated for the corresponding parameter. In addition, the amount of running time
  used to calculate each variable is returned.
}
\examples{
  # Do a gibbs sampling for an AR(1) model:
  #   y_t = phi*y_{t-1} + eps_t
  #   eps_t ~ normal(0, s2a)
  #   phi ~ I_[-1,1](phi)*N(0,1)
  #   s2a ~ inv_gamma(a, b)

  # doing a grid sample along n_grid points; the prior
  # is a truncated normal(0, 1)
  sample_phi = function(y, s2a, n_grid = 1000)
  {
    T = length(y)
    log_prob = function(phi)
    {
      - 0.5/s2a*(
      sum((y[2:T] - phi*y[1:(T-1)])^2)) - 0.5*phi^2
    }

    xs = seq(-0.999, 0.999, length.out = n_grid)
    prob = exp(sapply(xs, log_prob))
    sample(xs, 1, prob = prob)
  }

  # s2a has a inv gamma (a, b) prior; using conjugacy
  sample_s2a = function(y, phi, a, b)
  {
    T = length(y)

    1/rgamma(1,
      shape = a + T/2,
      rate = b + 0.5*sum((y[2:T] - c(phi)*y[1:(T-1)])^2))
  }

  n = 1000
  niter = 500
  ret = gibbs(
    niter,
    init                 = list(s2a = 1, phi = 0),
    hypers               = list(a = 0, b = 0),
    known_data           = list(y = rnorm(n)),
    conditional_samplers = list(s2a = sample_s2a, phi = sample_phi),
    ignore               = "s2a")

  plot(ret$phi)
}
