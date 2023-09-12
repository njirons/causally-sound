compute_stats <- function(samps){
  samps <- samps[names(samps) %in% c("theta0", "theta1", "etae", "etas")]
  samps$rr <- with(samps, theta1/theta0)
  samps$ef <- with(samps, 1-theta1/theta0)
  samps$or <- with(samps, (theta1/(1-theta1))/((theta0)/(1-theta0)))
  samps$rd <- with(samps, theta1-theta0)
  return(samps)
}

compute_summaries <- function(samps, med, lw, up){
  f <- function(x) c(mean(x), quantile(x, c(med, lw, up)))
  post.summaries <- t(sapply(samps, f))
  colnames(post.summaries) <- c("mean", "med", "lw", "up")
  return(post.summaries)
}

brease <- function(y0, y1, N0, N1, mu0, n0, mue, ne, mus, ns,
                   mono = F,
                   lw = 0.025, up = 0.975, med = .5,
                   n.iter = 1e4){
  # posterior samples
  if(!mono){
    samps <- sample_brease(1e5, y0 = y0, y1 = y1,N0 =  N0, N1 = N1,
                        mu0 = mu0, n0 = n0, mue = mue, ne = ne, mus = mus, ns = ns)
  } else {
    samps <- sample_brease_noharm(1e5, y0 = y0, y1 = y1,N0 =  N0, N1 = N1,
                               mu0 = mu0, n0 = n0, mue = mue, ne = ne)
  }

  # compute RR, OR, RD, EF
  samps <- compute_stats(samps)

  # posterior summaries
  post.summaries <- compute_summaries(samps, med = med, lw = lw, up = up)

  # Bayes Factors
  if(!mono){
    bf01_brease <- exp(bf_brease(y0,y1,N0,N1, mu0 = mu0, n0 = n0, mue = mue, ne = ne, mus = mus, ns = ns))
  } else {
    bf01_brease <- exp(bf_brease_mono(y0,y1,N0,N1, mu0 = mu0, n0 = n0, mue = mue, ne = ne))
  }
  bf10 <- 1/bf01_brease
  bayes.factors <- list(bf10 = bf10)

  out <- list(post.samples = samps,
              post.summaries = post.summaries,
              bayes.factors = bayes.factors)
  return(out)

}




brease.IB <- function(y0, y1, N0, N1, a= c(1,1), b= c(1,1),
                   lw = 0.025, up = 0.975, med = .5,
                   n.iter = 1e4){

  # posterior samples
  samps <- list(theta0 = NA, theta1 = NA)
  samps$theta0 <- rbeta(1e5, y0 + a[1], N0 - y0 + b[1])
  samps$theta1 <- rbeta(1e5, y1 + a[2], N1 - y1 + b[2])

  # compute RR, OR, RD, EF
  samps <- compute_stats(samps)

  # posterior summaries
  post.summaries <- compute_summaries(samps, med = med, lw = lw, up = up)

  # Bayes Factors
  bf10 <-  exp(get_bfindep_anal(y0,y1,N0,N1,1,1))
  bayes.factors <- list(bf10 = bf10)

  out <- list(post.samples = samps,
              post.summaries = post.summaries,
              bayes.factors = bayes.factors)
  return(out)

}


brease.LT <- function(y0 = y0, y1 = y1,N0 =  N0, N1 = N1,
                      mu_psi = 0, sigma_psi = 1,
                      mu_beta = 0, sigma_beta = 1,
                      lw = 0.025, up = 0.975, med = .5,
                      n.iter = 1e4){
  ab <- ab_test(data = data.frame(y1 = y0, n1= N0, y2=y1, n2= N1),
                nsamples = n.iter,
                prior_par = list(mu_psi = mu_psi, sigma_psi = sigma_psi,
                                 mu_beta= mu_beta, sigma_beta = sigma_beta))
  samps <- get_post_samples(ab)
  samps <- samps[c("p1", "p2")]
  colnames(samps) <- c("theta0", "theta1")

  # compute RR, OR, RD, EF
  samps <- compute_stats(samps)

  # posterior summaries
  post.summaries <- compute_summaries(samps, med = med, lw = lw, up = up)

  # Bayes Factors
  bf10 <-  exp(ab$logbf$bf10)
  bayes.factors <- list(bf10 = bf10)

  out <- list(post.samples = samps,
              post.summaries = post.summaries,
              bayes.factors = bayes.factors)
  return(out)
}
