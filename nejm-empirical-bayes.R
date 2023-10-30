rm(list=ls())

source('functions/helpers.R')
source('functions/brease_func.R')
source('functions/ib_func.R')
source('functions/lt_func.R')

dat <- read.csv('data/NEJM.csv')

# brease bounds
bounds_etae <- apply(dat,1,function(x){
  bounds <- brease_bounds(x[4],x[5],x[6],x[7])
  return(bounds$etae_m)  
})

bounds_etas <- apply(dat,1,function(x){
  bounds <- brease_bounds(x[4],x[5],x[6],x[7])
  return(bounds$etas_m)  
})

# BFs
bfindep_dat <- sapply(seq(1, 5, length.out = 1), function(a) {
  exp(-get_bfindep_anal(dat$y1, dat$y2, dat$n1, dat$n2, a, a))
})

bfdep_dat <- sapply( seq(1, 2, length.out = 1), function(sigma) {
  apply(dat, 1, function(x) {
    1/get_ab(x[4], x[5], x[6], x[7], sigma_psi = sigma)$bf$bf10
})})

bfdir_dat <- sapply(seq(0.3,length.out=1), function(a) {
  apply(dat, 1, function(x) {
    return(exp(bf_brease(x[4],x[5],x[6],x[7],
                    mu0 = 1/2, n0 = 2,
                    mue = a, ne = 1,
                    mus = a, ns = 1)))
  })
})

bfeb_dat <- apply(dat,1,function(x){
  bounds <- brease_bounds(x[4],x[5],x[6],x[7])
  return(exp(bf_brease(x[4],x[5],x[6],x[7],
                  mu0 = 1/2, n0 = 2,
                  mue = bounds$etae_m, ne = 1,
                  mus = bounds$etas_m, ns = 1)))  
})

# LMLs
lmlindep_dat <- sapply(seq(1, 5, length.out = 1), function(a) {
  c(lml1_IB(dat$y1, dat$y2, dat$n1, dat$n2, a, a))
})


lml_lt_thames <- function(y0, y1, N0, N1, sigma.psi){
  capture.output(lt_samps <- jags.lt(y0,y1,N0,N1, sigma.psi = sigma.psi, n.iter=1e4))
  params <- as.matrix(lt_samps[,1:2])
  lps <- apply(params,1,function(x){lp_lt(y0,y1,N0,N1,x[1],x[2])})
  capture.output(lt_thames <- thames(lps,params))
  -lt_thames$log_zhat_inv
}


lmldep_dat <- sapply( seq(1, 2, length.out = 1), function(sigma) {
  apply(dat, 1, function(x) lml_lt_thames(y0 = x[4], y1 = x[5], N0 = x[6], N1 = x[7], sigma.psi = sigma))
})


lmldir_dat <- sapply(seq(0.3,length.out=1), function(a) {
  apply(dat, 1, function(x) {
    return(lml1_brease(x[4],x[5],x[6],x[7],
                    mu0 = 1/2, n0 = 2,
                    mue = a, ne = 1,
                    mus = a, ns = 1))
  })
})

lmleb_dat <- apply(dat,1,function(x){
  bounds <- brease_bounds(x[4],x[5],x[6],x[7])
  return(lml1_brease(x[4],x[5],x[6],x[7],
                  mu0 = 1/2, n0 = 2,
                  mue = bounds$etae_m, ne = 1,
                  mus = bounds$etas_m, ns = 1))  
})

# plot
pdf("plots/nejm-eb-comp.pdf", width = 6, height = 6)

cols <- brewer.pal(3, 'Set1')
cols[3] <- rgb(255/255,140/255,0,1)
cols[4] <- rgb(255/255,165/255,0,0.5)
cex.main  <- 1.1
cex.lab   <- 1.1
cex.axis  <- 1.1

par(mar = c(5.2, 5.2, 1, 0), mfrow=c(1,2))

# plot BFs
ylim <- range(c(bfdep_dat, bfdir_dat,bfeb_dat))
# ylim <- range(c(bfindep_dat,bfdep_dat, bfdir_dat,bfeb_dat))
plot(1:nrow(dat), bfdep_dat,
     col = cols[2],
     pch = 20, axes = FALSE, cex = 1.5,
     xlab = "Study ID",
     ylab = "",
     # ylim = ylim,
     ylim=c(1,20),
     log = 'y',
     font.main = 1,
     cex.main = cex.main, cex.lab = cex.lab)
axis(1, cex.axis = cex.axis, at = c(1, seq(5, 35, 5), 39))
axis(2, las = 2, cex.axis = cex.axis, at = c(1,3,10,20))
mtext(text = "Evidence for H0", side = 2, line = 2.5, cex = cex.lab)
# points(bfindep_dat, pch = 18, cex = 1.4, col = cols[1])
points(bfdir_dat, pch = 17, cex = 1.2, col = cols[3])
points(bfeb_dat, pch = 17, cex = 1.2, col = cols[4])

# plot LMLs
ylim <- range(c(lmldep_dat, lmlindep_dat, lmldir_dat,lmleb_dat))
plot(1:nrow(dat), lmldep_dat,
     col = cols[2],
     pch = 20, axes = FALSE, cex = 1.5,
     xlab = "Study ID",
     ylab = "",
     # ylim = ylim,
     ylim = c(-19.5,-4.0),
     font.main = 1,
     cex.main = cex.main, cex.lab = cex.lab)
axis(1, cex.axis = cex.axis, at = c(1, seq(5, 35, 5), 39))
axis(2, las = 2, cex.axis = cex.axis, at = round(seq(-19.5, -5, len=7), 1))
mtext(text = "Log Marginal Likelihood", side = 2, line = 3.5, cex = cex.lab)
# points(lmlindep_dat, pch = 18, cex = 1.4, col = cols[1])
points(lmldir_dat, pch = 17, cex = 1.2, col = cols[3])
points(lmleb_dat, pch = 17, cex = 1.2, col = cols[4])
legend("bottomright",
       col = cols[c(2,3,4)],
       pch = c(20, 17,17),
       legend = c("LT", "BREASE",'BREASE-EB'),
       horiz = F, bty = "n",
       cex= cex.main)

dev.off()