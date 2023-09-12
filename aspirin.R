### aspirin example
rm(list=ls())
source('functions/ib_func.R')
source('functions/brease_func.R')
source('functions/brease.R')
source('functions/contour_func.R')
source('functions/helpers.R')
library(abtest)
library(epitools)
pacman::p_load(bridgesampling,rstan)
options(mc.cores = parallel::detectCores())
rstan_options(autowrite=TRUE)

# data
y0 <- 26
N0 <- 11034
y1 <- 10
N1 <- 11037

mu0 <- .5;  n0 <- 2;
mue <- .3;  ne <- mu0*n0;
mus <- .3;  ns <- (1-mu0)*n0;




# BREASE ------------------------------------------------------------------

results <- brease(y0, y1, N0, N1, mu0 = mu0, n0 = n0, mue = mue, ne = ne, mus = mus, ns = ns)
round(results$post.summaries*100, 2)
results$bayes.factors$bf10

# mono
mono <- brease(y0, y1, N0, N1, mu0 = mu0, n0 = n0, mue = mue, ne = ne, mus = mus, ns = ns, mono = T)
round(mono$post.summaries*100, 2)
mono$bayes.factors$bf10

# PNS side-effect plot
side_plot(y0, y1, N0, N1, mu0, n0, mue, ne, log = F)

# PNS contour plot
brease.cont <- bf_seq_brease(y0, y1, N0, N1, mu0, n0, ne, ns)
contours(x = brease.cont$mus.seq,
         y = brease.cont$mue.seq,
         z = brease.cont$log.bf.seq)



# IB and LT ---------------------------------------------------------------

# IB 
ib.results <- brease.IB(y0, y1, N0, N1)
round(ib.results$post.summaries*100, 2)
bf01 <- 1/ib.results$bayes.factors$bf10
bf01

# LT
lt.results <- brease.LT(y0, y1, N0, N1)
round(lt.results$post.summaries*100, 2)
lt.results$bayes.factors$bf10
