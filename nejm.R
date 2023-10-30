rm(list=ls())

source('functions/helpers.R')
source('functions/brease_func.R')
source('functions/ib_func.R')
source('functions/lt_func.R')

dat <- read.csv('data/NEJM.csv')

# log marginal likelihood comparison
pdf("plots/nejm-lml-comp.pdf", width = 6, height = 6)

bfindep_dat <- sapply(seq(1, 5, length.out = 1), function(a) {
  c(lml1_IB(dat$y1, dat$y2, dat$n1, dat$n2, a, a))
})


lml_lt_thames <- function(y0, y1, N0, N1, sigma.psi){
  capture.output(lt_samps <- jags.lt(y0,y1,N0,N1, sigma.psi = sigma.psi, n.iter=1e4))
  params <- as.matrix(lt_samps[,1:2])
  lps <- apply(params,1,function(x){lp_lt(y0,y1,N0,N1,x[1],x[2])})
  capture.output(lt_thames <- thames(lps,params))
  -lt_thames$log_zhat_inv
}


bfdep_dat <- sapply( seq(1, 2, length.out = 1), function(sigma) {
  apply(dat, 1, function(x) lml_lt_thames(y0 = x[4], y1 = x[5], N0 = x[6], N1 = x[7], sigma.psi = sigma))
})


bfdir_dat <- sapply(seq(0.3,length.out=1), function(a) {
  apply(dat, 1, function(x) {
    return(lml1_brease(x[4],x[5],x[6],x[7],
                  mu0 = 1/2, n0 = 2,
                  mue = a, ne = 1,
                  mus = a, ns = 1))
  })
})

cols <- brewer.pal(3, 'Set1')
cols[3] <- "darkgreen"
cex.main  <- 1.1
cex.lab   <- 1.1
cex.axis  <- 1.1
ylim <- range(c(bfdep_dat, bfindep_dat, bfdir_dat))

par(mar = c(5.2, 5.2, 4.1, 2.1))
plot(1:nrow(dat), bfdep_dat,
     col = cols[2],
     pch = 20, axes = FALSE, cex = 1.5,
     xlab = "Study ID",
     ylab = "",
     ylim = ylim,
     font.main = 1,
     cex.main = cex.main, cex.lab = cex.lab)
axis(1, cex.axis = cex.axis, at = c(1, seq(5, 35, 5), 39))
axis(2, las = 2, cex.axis = cex.axis, at = round(seq(-19.5, -5, len=7), 1))
mtext(text = "Log Marginal Likelihood", side = 2, line = 3.5, cex = cex.lab)
points(bfindep_dat, pch = 18, cex = 1.4, col = cols[1])
points(bfdir_dat, pch = 17, cex = 1.2, col = cols[3])
legend("topleft",
       col = cols[c(1,2,3)],
       pch = c(18,20, 17),
       legend = c("IB", "LT", "BREASE"),
       horiz = T, bty = "n",
       cex= cex.main)

dev.off()

# bayes factor comparison
pdf("plots/nejm-bf-comp.pdf", width = 6, height = 6)

bfindep_dat <- sapply(seq(1, 5, length.out = 100), function(a) {
  get_bfindep_anal(dat$y1, dat$y2, dat$n1, dat$n2, a, a)
})

bfdep_dat <- sapply(seq(1, 4, length.out = 100), function(sigma) {
  apply(dat, 1, function(x) (get_ab(x[4], x[5], x[6], x[7], sigma_psi = sigma)$bf$bf10))
})


library(parallel)

# detect the number of cores
n_cores <- detectCores()

# create a cluster with n_cores
cl <- makeCluster(n_cores)

# export necessary objects to the cluster
clusterExport(cl, ls())

# apply function in parallel
seq_mu <- pretty(c(0.2,.7),n = 100)
seq_mu <- seq_mu[1:100]
mu_def <- which(abs(seq_mu-.3)<0.001)

bfdir_dat <- parSapply(cl, seq_mu, function(a) {
  apply(dat, 1, function(x) {
    return(bf_brease(x[4],x[5],x[6],x[7],
                     mu0 = 1/2, n0 = 2,
                     mue = a, ne = 1,
                     mus = a, ns = 1))
  })
})

# stop the cluster
stopCluster(cl)


cols <- brewer.pal(3, 'Set1')
cols[3] <- "darkgreen"
cex.main  <- 1.1
cex.lab   <- 1.1
cex.axis  <- 1.1
cexx      <- 1.1

alphas_indep <- c(1, seq(0.05, 0.010, length.out = 99))
cols_indep <- sapply(alphas_indep, function(alpha) {
  adjustcolor(cols[1], alpha.f = alpha)
})

alphas_dep <- c(1, seq(0.05, 0.010, length.out = 99))
cols_dep <- sapply(alphas_dep, function(alpha) {
  adjustcolor(cols[2], alpha.f = alpha)
})

alphas_dir <- c(seq(0.01, 0.07, length.out = mu_def-1),
                1,
                seq(0.07, 0.010, length.out =100-mu_def))
cols_dir <- sapply(alphas_dir, function(alpha) {
  adjustcolor(
    cols[3],
    alpha.f = alpha)
})

par(mar = c(5.1, 4.1, 4.1, 2.1))

range.bf <- range(c(exp(-bfindep_dat),(1/bfdep_dat),exp(bfdir_dat) ))
plot(
  seq(1, 39), exp(-bfindep_dat[, 1]), pch = 18, axes = FALSE, cex = 1.4, log = 'y',
  xlab = 'Study ID', ylab = '', ylim = c(1, 200), col = cols[1],
  font.main = 1,
  cex.main = cex.main, cex.lab = cex.lab
)
axis(1, cex.axis = cex.axis, at = c(1, seq(5, 35, 5), 39))
axis(2, las = 2, cex.axis = cex.axis, at = c(1, 3, 10, 20, 100, 150, 200))

for (i in seq(2, 100)) {
  points(
    exp(-bfindep_dat[, i]), pch = 18, cex = 1.4,
    col = cols_indep[i]
  )
}

for (i in seq(1, 100)) {
  points(
    (1/bfdep_dat[, i]), pch = 20, cex = 1.5,
    col = cols_dep[i]
  )
}


for (i in seq(1, 100)) {
  points(
    exp(bfdir_dat[, i]), pch = 17, cex = 1.2,
    col = cols_dir[i]
  )
}



mtext(text = "Evidence for H0", side = 2, line = 2, cex = cex.lab)

points_x <- seq(33, 37.5, length.out = 100)
points_y <- rep(200, 100)
off <- 70


text(30, 200, 'IB', cex = cexx)
points(points_x, points_y, pch = 18, cex = 1.8, col = cols_indep)
text(31.75, 200, '1', cex = cexx)
text(38.75, 200, '5', cex = cexx)

text(30, 200 - off, 'LT', cex = cexx)
points(points_x, points_y - off, pch = 20, cex = 2, col = cols_dep)
text(31.75, 200 - off, '1', cex = cexx)
text(38.75, 200 - off, '2', cex = cexx)

ydir <- exp(log(200-off)- (log(200) - log(200-off)))
text(28.25, ydir, 'BREASE', cex = cexx)
points(points_x, rep(ydir, 100), pch = 17, cex = 1.4, col = cols_dir)
text(31.7,  ydir, '.2', cex = cexx)
text(38.75, ydir, '.7', cex = cexx)

dev.off()
