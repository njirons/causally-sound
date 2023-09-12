rm(list = ls())
source('functions/ib_func.R')
source('functions/brease_func.R')
source('functions/brease.R')
source('functions/contour_func.R')
source('functions/helpers.R')
library(dplyr)
library(latex2exp)
library(rjags)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

mu0 <- 1/2
n0 <- 2

# eta_e needs to be uniform
mue <- 1/2
ne <- 2

mus <- 1/100
ns <- 1

# data
y0 <- 20
N0 <- 1000
y1 <- 40
N1 <- 1000

# sample
n.iter <- 1e5
jags.fit <- jags.brease(y0=y0, y1=y1, N0=N0, N1=N1, mus = mus, ns = ns, ne = ne, n.iter = n.iter)
our.fit  <- sample_brease(n.iter, y0 = y0,y1 =  y1,N0 = N0,N1 =  N1, mus = mus, ns = ns, ne = ne)
stan.fit <- stan.brease(y0=y0, y1=y1, N0=N0, N1=N1, mus = mus, ns = ns, ne = ne, iter = n.iter/4)

par(
  mfrow=c(1,2),
  cex=1,
  mar = c(4,4.5,2.5,2),
  font=1, font.main=2, font.lab=1,font.axis=1
)
hist(our.fit$theta0, freq = F, breaks = 100, ylim = c(0, 150),
     xlab = TeX('$\\theta_0$',bold=TRUE),
     col= "gray",
     axes = F,
     main=NULL,
     ylab=NULL)
axis(1)
hist(jags.fit$theta0, col= rgb(0,0,1,0.3), add= T, freq = F, breaks = 100)
hist(stan.fit$samples$theta0, col= rgb(1,0,0,0.3), add= T, freq = F, breaks = 100)

par(mar = c(4,2,2.5,0))
hist(our.fit$theta1, freq = F, breaks = 100, ylim = c(0, 150),
     xlab = TeX('$\\theta_1$',bold=TRUE),
     col= "gray", border = T, axes = F,
     main=NULL,
     ylab=NULL)
axis(1)
hist(jags.fit$theta1, col= rgb(0,0,1,0.3), add= T, breaks = 100, freq = F, border = T)
hist(stan.fit$samples$theta1, col= rgb(1,0,0,0.3), add= T, breaks = 100, freq = F, border = T)
legend('topright',
       legend=c('Exact','JAGS','Stan'),
       col=c(col= "gray",rgb(0,0,1,0.3),rgb(1,0,0,0.3)),
       lty = rep(1,3), lwd = rep(5,3))
