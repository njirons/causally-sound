

bf_seq_brease <- function(y0, y1, N0, N1, mu0, n0, ne, ns,
                       mus.seq = seq(0.01, 0.99, length = 20),
                       mue.seq = seq(0.01, 0.99, length = 20),
                       parallel = T){
  library(doParallel)
  if(parallel){
    cores <- detectCores()
    registerDoParallel(cores = cores)
  }
  f <- function(x, y) -(bf_brease(y0,y1,N0,N1, mu0 = mu0, n0 = n0, mue = x, ne = ne, mus = y, ns = ns))
  z <- foreach(x = mue.seq, .combine = 'cbind') %:%
    foreach(y = mus.seq, .combine = 'c')   %dopar% {
      f(x, y)
    }
  out <- list(mus.seq = mus.seq,
              mue.seq = mue.seq,
              log.bf.seq  = z)
}

bf_seq_ab <- function(y0, y1, N0, N1,
                      a.seq =  seq(0.5,5, length.out = 100),
                      b.seq =  a.seq
                      ){
  f <- function(x, y) get_bfindep_anal(y0,y1,N0,N1,a=x,b=y)
  f <- Vectorize(f)
  z <- outer(a.seq, b.seq, f)
  out <- list(a.seq = a.seq,
              b.seq = b.seq,
              log.bf.seq  = z)
}

bf_seq_lt <- function(y0, y1, N0, N1,
                      mu_range = c(-5,3.5),
                      sigma_range=c(0.01,2),
                      mu_steps = 100,
                      sigma_steps = 100){
  data <- list(y1 = y0, n1 = N0, y2 = y1, n2 = N1)
  ab <- ab_test(data = data)
  p <- plot_robustness(ab,
                       mu_range = mu_range,
                       sigma_range= sigma_range,
                       cores = parallel::detectCores(),
                       mu_steps = mu_steps,
                       sigma_steps = sigma_steps)
  mu.seq <- seq(range(p$mu_psi)[1],range(p$mu_psi)[2],length.out=mu_steps)
  sigma.seq <- seq(range(p$sigma_psi)[1],range(p$sigma_psi)[2],length.out=sigma_steps)
  z.lt <- matrix(p$bf,nrow=mu_steps,ncol=sigma_steps,byrow=FALSE)
  log.z.lt <- log(z.lt)
  out <- list(mu.seq = mu.seq,
              sigma.seq = sigma.seq,
              log.bf.seq  = log.z.lt)
}


contours <- function(x, y, z, n.levels = 20,
                     xlab = "Prior Expected Risk of Side-Effects",
                     ylab = "Prior Expected Efficacy",
                     cex = 2,
                     cex.lab = .8,
                     thr = c(1, 3, 10),
                     thr.cols = c("red", "orange", "blue")){

  # setting graphic parameters
  oldpar <- par(list(mar = c(4,4,1,1), pty = "s"))
  on.exit(par(oldpar))

  n.levels <- 20
  levels <- pretty(z, n = n.levels)
  for(j in seq_along(thr)){
    levels <- levels[!abs((levels) - log(thr[j]))  < (0.15)]
  }
  contour(x, y, z, axes = F,
          col= "grey40",
          labels = formatC(exp(levels),2),
          levels = levels,
          cex = cex,
          cex.lab = cex.lab,
          xlab = xlab,
          ylab = ylab)
  for(j in seq_along(thr)){
    contour(x, y, z,
            cex = 2,
            levels = log(thr[j]), label = thr[j],
            col = thr.cols[j], lty = 2, add = T, lwd = 2)
  }
  box()
  axis(1, at = pretty(x,n=10), cex.axis = .8)
  axis(2, at = pretty(y,n=10), cex.axis = .8)
}


side_plot <- function(y0, y1, N0, N1, mu0, n0, mue, ne,
                      mus.seq = seq(0.001, 0.5, length = 25),
                      ns = c(1/2, 1, 2),
                      lty = c(2, 1, 3),
                      log = F
){
  bf.seq <- list()
  for (i in seq_along(ns)){
    bf.seq[[i]] <- bf_seq_brease(y0, y1, N0, N1, mu0, n0, ne,
                              mue.seq = mue,
                              ns = ns[i],
                              mus.seq = mus.seq)
  }

  oldpar <- par(list(mar = c(4,4,1,1), pty = "s"))
  on.exit(par(oldpar))
  f  <- function(x) x
  fm1 <- exp
  if(!log) {
    f <- exp
    fm1 <- function(x)x
  }
  range <- range(sapply(bf.seq, function(x) f(x$log.bf.seq)))
  plot(mus.seq, f(bf.seq[[1]]$log.bf.seq),
       ylim = range*c(1, 1.1),
       type = "n", axes = F,
       xlab = "Prior Expected Risk of Side-Effects",
       ylab = ('Evidence against H0'),
       cex.lab = .8)
  axis(1, at = pretty(mus.seq, n=10),
       cex.axis = .8)
  yaxis <- pretty(range, n=10)
  axis(2, at = yaxis,
       labels = formatC(fm1(yaxis),2),
       cex.axis = .8)
  abline(h=f(log(1)), col= "red", lty = 2, lwd = 2)
  abline(h=f(log(3)), col= "orange", lty = 2, lwd = 2)
  abline(h=f(log(10)), col= "blue", lty = 2, lwd = 2)
  for (i in seq_along(ns)){
    lines(mus.seq, f(bf.seq[[i]]$log.bf.seq), type = "l", lty = lty[i])
  }

  legend("topright", xpd = TRUE, bty = "n",  ncol = 3,
         title = "Prior sample size", cex = .8,
         legend = ns, lty = lty,inset = c(0,0))
}



