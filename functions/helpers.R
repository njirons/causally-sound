# library('rjags')
library('abtest')
library('truncnorm')
library('latex2exp')
library('logspline')
# library('doParallel')
library('RColorBrewer')


# Independent Beta Bayes factor
get_bfindep_anal <- function(y1, y2, n1, n2, a = 1, b = 1) {
  post0 <- lbeta(2 * a - 1, 2 * b - 1) - 2 * lbeta(a, b)
  post1 <- (
    lbeta(2 * a + y1 + y2 - 1, 2 * b + n1 - y1 + n2 - y2 - 1) -
    lbeta(a + y1, b + n1 - y1) - lbeta(a + y2, b + n2 - y2)
  )

  post0 - post1
}
# get_bfindep_anal <- function(y1, y2, n1, n2, a = 1, b = 1) {
#   post0 <- lbeta(a,b) - 2 * lbeta(a, b)
#   post1 <- (
#     lbeta(a + y1 + y2, b + n1 - y1 + n2 - y2) -
#       lbeta(a + y1, b + n1 - y1) - lbeta(a + y2, b + n2 - y2)
#   )
#   
#   post0 - post1
# }

get_h0_lml <- function(y1, y2, n1, n2, a = 1, b = 1) {
  lchoose(n1,y1) + lchoose(n2,y2) +
    lbeta(a + y1 + y2, b + n1 - y1 + n2 - y2) - lbeta(a,b)
}

get_bfindep_anal_emp <- function(y1, y2, n1, n2, a = 1,b=1) {
  post0 <- lbeta((y1/n1)*a+(y2/n2)*b - 1, 
                 (1-y1/n1)*a+(1-y2/n2)*b - 1) - 
    lbeta((y1/n1)*a, (1-y1/n1)*a) - lbeta((y2/n2)*b, (1-y2/n2)*b)
  post1 <- (
    lbeta((y1/n1)*a+(y2/n2)*b + y1 + y2 - 1, 
          (1-y1/n1)*a+(1-y2/n2)*b + n1 - y1 + n2 - y2 - 1) -
      lbeta((y1/n1)*a + y1, (1-y1/n1)*a + n1 - y1) - 
      lbeta((y2/n2)*b + y2, (1-y2/n2)*b + n2 - y2)
  )
  
  post0 - post1
}


get_ab <- function(
  y1, y2, n1, n2, sigma_psi = 1, sigma_beta = 1, samples = 10000
  ) {
  dat <- list(y1 = y1, y2 = y2, n1 = n1, n2 = n2)
  
  priors <- list(
    mu_psi = 0, sigma_psi = sigma_psi,
    mu_beta = 0, sigma_beta = sigma_beta
  )
  
  ab_test(dat, prior_par = priors, nsamples = samples)
}


# Logit Transformation Bayes factor
# get_bfab <- function(y1, y2, n1, n2, sigma_psi = 1, sigma_beta = 1) {
#   m <- get_ab(y1, y2, n1, n2, sigma_psi = sigma_psi, sigma_beta = sigma_beta)
#   m$logbf$bf10
# }
get_bfab <- function(y1, y2, n1, n2, sigma_psi = 1, sigma_beta = 1) {
  m <- get_ab(y1, y2, n1, n2, sigma_psi = sigma_psi, sigma_beta = sigma_beta)
  m$logml$logml1 - get_h0_lml(y1, y2, n1, n2,a=1,b=1)
}


# Bayesian credible intervals
get_post <- function(y1, y2, n1, n2, a = 1, b = 1, sigma_psi = 1, samples = 5e4) {
  m <- get_ab(y1, y2, n1, n2, samples = samples, sigma_psi = sigma_psi)
  samples <- get_post_samples(m)
  psi_ab <- samples$psi
  eta_ab <- samples$arisk
  
  theta1 <- rbeta(samples, a + y1, b + n1 - y1)
  theta2 <- rbeta(samples, a + y2, b + n2 - y2)
  logit <- function(x) log(x / (1 - x))
  
  psi_ct <- logit(theta2) - logit(theta1)
  eta_ct <- theta2 - theta1
  
  ci_ab_psi <- quantile(psi_ab, c(0.025, 0.975))
  ci_ab_eta <- quantile(eta_ab, c(0.025, 0.975))
  
  ci_ct_psi <- quantile(psi_ct, c(0.025, 0.975))
  ci_ct_eta <- quantile(eta_ct, c(0.025, 0.975))
  
  cbind(
    ci_ab_psi, ci_ct_psi, ci_ab_eta, ci_ct_eta,
    'pmeans_psi' = c(mean(psi_ab), mean(psi_ct)),
    'pmeans_eta' = c(mean(eta_ab), mean(eta_ct))
  )
}


# Plot joint prior for the LT approach
plot_joint_prior <- function(
  cex_lab,
  prior_par = list(mu_psi = 0, sigma_psi = 1, mu_beta = 0, sigma_beta = 1),
  length = 100, ...
) {
  xx <- seq(0, 1, length.out = length)
  yy <- xx
  zz <- outer(xx, yy, dprior, prior_par = prior_par, what = 'p1p2')
  image(
    xx, yy, zz, xlab = '', ylab = '', las = 1,
    col = hcl.colors(500, 'YlOrRd', rev = TRUE), font.main = 1, ...
  )
  mtext(text = expression(theta[1]), side = 1, line = 2.8, cex = cex_lab)
}


# Plot conditional prior for the LT approach
plot_conditional_prior <- function(
  p1, cex_lab, ylim, yticks,
  prior_par = list(mu_psi = 0, sigma_psi = 1, mu_beta = 0, sigma_beta = 1), ...
) {
  col_trans <- rgb(0, 0, 0, alpha = 0.15)
  xx <- seq(0, 1, length.out = 100)
  yy <- dprior(x1 = xx, x2 = p1, prior_par = prior_par, what = 'p2givenp1')
  
  xticks <- pretty(c(0, 1))
  xlim <- range(xticks)
  
  plot(
    1, axes = FALSE, type = 'n', xlim = c(0, 1),
    xlab = '', ylab = '', font.main = 1, ylim = ylim, ...
  )
  lines(xx, yy, lwd = 0, ...)
  graphics::polygon(
    c(xx, xx[length(xx)], xx[1]), y = c(yy, rep(ylim[1], 2)),
    col = adjustcolor(cols[2], 0.50), border = NA
  )
  axis(1, at = xticks, ...)
  axis(2, at = yticks, las = 1, ...)
  mtext(
    text = bquote(theta[2] ~ ' | ' ~ theta[1] ~ ' = ' ~ .(p1)),
    side = 1, line = 3.5, cex = cex_lab
  )
}


# Plot joint prior for the IB approach
plot_beta_joint <- function(
  a, cex_lab = cex.lab, length = 100,
  main = 'Joint Prior Distribution', ...
  ) {
  xx <- seq(0, 1, length.out = length)
  zz <- outer(xx, xx, prior, a, a)
  image(
    xx, xx, zz, xlab = '', ylab = '', las = 1,
    col = hcl.colors(500, 'YlOrRd', rev = TRUE),
    main = main, font.main = 1, ...
  )
  mtext(text = expression(theta[1]), side = 1, line = 2.8, cex = cex_lab)
}


# Plot conditional prior for the IB approach
plot_beta_conditional <- function(
  a, theta1 = 0.10, cex_lab = cex.lab,
  length = 100, main = 'Conditional Prior Distribution',
  ...) {
  xx <- seq(0, 1, length.out = length)
  plot(
    xx, dbeta(xx, a, a), lwd = 0, type = 'l', axes = FALSE,
    main = main, font.main = 1, xlab = '', ylab = '', ...
  )
  axis(1, cex.axis = cex.axis)
  axis(2, las = 2, cex.axis = cex.axis)
  mtext(
    text = bquote(theta[2] ~ ' | ' ~ theta[1] ~ ' = ' ~ .(theta1)),
    side = 1, line = 3.5, cex = cex_lab
  )
}
