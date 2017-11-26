horseshoe_regression_plot <- function(
  beta,
  Y,
  X,
  main = "Horseshoe Prior",
  ylim = NULL)
{
  mod_ols = lm(y ~ X - 1)
  summary(mod_ols)

  p <- length(X[1,])
  ci_beta = HPDinterval(as.mcmc(beta))
  beta_means = colMeans(beta)
  beta_ols = mod_ols$coefficients

  if(is.null(ylim))
  {
    ylim = range(ci_beta, beta_ols, na.rm=TRUE)
  }

  # circles = beta ols coefficients
  plot(1:p, beta_ols,
    ylim = ylim,
    main = main,
    xlab = "beta indices")

  for(j in 1:p)
  {
    # blue line marks 95% hpd
    lines(rep(j, 2), ci_beta[j,], col='blue', lwd=3)
    # red x marks mean of estimate
    lines(j, mean(beta[,j]),
          type='p', pch=4, lwd=3, cex = 2, col='red')
  }
  abline(h = 0, lwd=3, col='green', lty=2)
}
