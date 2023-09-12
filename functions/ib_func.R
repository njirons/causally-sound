# calculate unnormalized log posterior, independent beta
logpost0_IB <- function(theta0,y0,y1,n0,n1,a0=1,b0=1){
  lp <- lchoose(n0,y0) + lchoose(n1,y1) +
    (y0+y1+a0-1)*log(theta0) + 
    (n0+n1-(y0+y1)+b0-1)*log(1-theta0) - lbeta(a0,b0)
  return(lp)
}

logpost1_IB <- function(theta0,theta1,y0,y1,n0,n1,a0=1,b0=1,a1=1,b1=1){
  lp <- lchoose(n0,y0) + lchoose(n1,y1) +
    (y0+a0-1)*log(theta0) + 
    (n0-y0+b0-1)*log(1-theta0) +
    (y1+a1-1)*log(theta1) + 
    (n1-y1+b1-1)*log(1-theta1) - lbeta(a0,b0) - lbeta(a1,b1)
  return(lp)
}

# independent beta logmarglik, H0
lml0_IB <- function(y0,y1,n0,n1,a0=1,b0=1){
  lml <- lchoose(n0,y0) + lchoose(n1,y1) +
    lbeta(y0+y1+a0,n0+n1-(y0+y1)+b0) - lbeta(a0,b0)
  return(lml)
}

# independent beta logmarglik, H1
lml1_IB <- function(y0,y1,n0,n1,a0=1,b0=1,a1=1,b1=1){
  lml <- lchoose(n0,y0) + lchoose(n1,y1) +
    lbeta(y0+a0,n0-y0+b0) - lbeta(a0,b0) +
    lbeta(y1+a1,n1-y1+b1) - lbeta(a1,b1)
  return(lml)
}

# independent beta bayes factor
get_bfindep <- function(y0,y1,n0,n1,a00=1,b00=1,a0=1,b0=1,a1=1,b1=1){
  return(lml0_IB(y0,y1,n0,n1,a00,b00) - lml1_IB(y0,y1,n0,n1,a0,b0,a1,b1))
}

# dependent IB model (see Dablander Appendix D) 
ib_dep_stan <- 'data {
  int<lower=1> N1; // number treated
  int<lower=1> N0; // number control
  int<lower=0,upper=N1> y1; // outcome, treated 
  int<lower=0,upper=N0> y0; // outcome, control 
  
  real<lower=0> sigma_eta;
  real<lower=0> sigma_zeta;
}
parameters {
  real<lower=-1,upper=1> eta; 
  real<lower=0,upper=1> zeta;
}
transformed parameters {
  real<lower=0,upper=1> theta0 = fmin(fmax(zeta-eta/2,0),1);
  real<lower=0,upper=1> theta1 = fmin(fmax(zeta+eta/2,0),1);
}
model {
  // priors
  eta ~ normal(0,sigma_eta)T[-1,1];
  zeta ~ normal(0,sigma_zeta)T[0,1];
  
  // (log) likelihood
  target += binomial_lpmf(y0 | N0, theta0);
  target += binomial_lpmf(y1 | N1, theta1);
}
'

# function to generate samples from ib_dep_stan and calculate LML
bf01_ib_dep <- function(y0,y1,N0,N1,sigma_eta = 0.2, sigma_zeta = 0.5,
                        model=NULL,iter = 2000){
  if(is.null(model)){
    model <- stan_model(model_code = ib_dep_stan)
  }
  
  # fit model
  ib_dep_data <- list(
    N1 = N1,
    y1 = y1,
    N0 = N0, 
    y0 = y0,
    sigma_eta = sigma_eta,
    sigma_zeta = sigma_zeta
  )
  ib_dep_fit <- sampling(model, data = ib_dep_data, iter = iter)
  ib_dep_bridge <- bridge_sampler(ib_dep_fit)
  return(get_h0_lml(y0, y1, N0, N1, a = 1, b = 1) - ib_dep_bridge$logml)
}