# Same DGP and RNG call order as the historical manuscript generator.
# mu_form='linear' retains its historical x1*x3 term; it is a label, not
# a claim that the prognostic function is additive in the original columns.
generate_data_medical <- function(n = 250, is_te_hetero = TRUE,
  is_mu_nonlinear = TRUE, seed = 1848, RCT = FALSE, scenario = "default",
  z_diff = FALSE, contrast_binary = TRUE, BCF = FALSE, sigma_sq = 1) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  x4 <- rbinom(n, 1, .5)
  if (contrast_binary) x4 <- 2*x4-1
  x5_raw <- sample(1:3, n, replace = TRUE, prob = rep(1/3, 3))
  g <- c(2,-1,-4)[x5_raw]
  mu <- if (!is_mu_nonlinear) 1+g+x1*x3 else -6+g+6*abs(x3-1)
  tau <- if (!is_te_hetero) rep(3,n) else switch(scenario,
    complex_interaction = 1+4*x1+3*x2+2*x2*x1,
    medical_saturation = .5+4/(1+exp(-2.5*x1)),
    mild_saturation = .5+4/(1+exp(-.4*x1)), 1+2*x2*x4)
  s <- sd(mu); u <- runif(n)
  pi_x <- if (RCT) rep(.5,n) else .8*pnorm(3*mu/s-.5*x1)+.05+u/10
  pi_x <- pmin(pmax(pi_x,0),1)
  z <- rbinom(n,1,pi_x); eps <- rnorm(n,0,sqrt(sigma_sq))
  z_binary <- z
  delta <- if (is.logical(z_diff)) if (z_diff) .5 else 0 else z_diff
  if (z_diff != FALSE) z <- z-delta
  y <- mu+z*tau+eps; y_hat <- mu+z*tau
  x5 <- factor(x5_raw, levels=1:3); contrasts(x5) <- contr.sum(3)
  xd <- model.matrix(~x5)[,-1,drop=FALSE]; colnames(xd) <- c("x5_1","x5_2")
  if (BCF) z <- z_binary
  data.frame(x1,x2,x3,x4,x5_1=xd[,1],x5_2=xd[,2],z,y,mu,pi_x,tau,y_hat)
}
