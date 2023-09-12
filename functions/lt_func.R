library(rjags)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source('functions/thames.R')
source('functions/trunc_quantile.R')

# sample from LT posterior using STAN
lt_stan_code <- 'data {
  int<lower=1> N1; // number treated
  int<lower=1> N0; // number control
  int<lower=0,upper=N1> y1; // outcome, treated 
  int<lower=0,upper=N0> y0; // outcome, control 
  
  real<lower=0> sigma_beta; // prior SD
  real<lower=0> sigma_psi; // prior SD
}
parameters {
  real beta;
  real psi;
}
transformed parameters {
  real<lower=0,upper=1> theta0 = exp(beta-psi/2)/(1+exp(beta-psi/2));
  real<lower=0,upper=1> theta1 = exp(beta+psi/2)/(1+exp(beta+psi/2));
}
model {
  // priors
  target += normal_lpdf(beta | 0, sigma_beta);
  target += normal_lpdf(psi | 0, sigma_psi);
  
  // (log) likelihood
  target += binomial_lpmf(y0 | N0, theta0);
  target += binomial_lpmf(y1 | N1, theta1);
}
'

stan.lt <- function(y0, y1, N0, N1, 
                     sigma.beta=1, sigma.psi=1, 
                     n.iter = 1e4){
  lt_data <- list(
    N1 = N1,
    y1 = y1,
    N0 = N0, 
    y0 = y0,
    sigma_beta = sigma.beta,
    sigma_psi = sigma.psi)
  
  lt_model <- stan_model(model_code = lt_stan_code)
  lt_fit <- sampling(lt_model, 
                      data = lt_data,
                      iter = n.iter)
  samples <- extract(lt_fit)
  out <- list(model = lt_fit, samples = samples)
  return(out)
}

# sample from LT posterior using JAGS
jags.lt <- function(y0,y1,N0,N1, 
                     sigma.beta=1, sigma.psi=1,
                     n.iter = 1e4){
  
  model_code <- "model{
  # likelihood
  y1 ~ dbinom(theta1, N1);
  y0 ~ dbinom(theta0, N0);

  # priors
  beta ~ dnorm(0, 1/(sigma.beta^2));
  psi ~ dnorm(0, 1/(sigma.psi^2));
  
  theta0 <- exp(beta-psi/2)/(1+exp(beta-psi/2));
  theta1 <- exp(beta+psi/2)/(1+exp(beta+psi/2));
}"
  data <- list(y1=y1, y0=y0, N1=N1, N0=N0,
               sigma.beta=sigma.beta, sigma.psi=sigma.psi)
  model <- jags.model(file = textConnection(model_code),data = data)
  update(model, n.iter = n.iter)
  samps <- coda.samples(model = model, 
                        variable.names = c("beta","psi","theta0",'theta1'),
                        n.iter = n.iter)
  samps <- data.frame(samps[[1]])
  return(samps)
}

jags.lt.h0 <- function(y0,y1,N0,N1, 
                    sigma.beta=1,
                    n.iter = 1e4){
  
  model_code <- "model{
  # likelihood
  y1 ~ dbinom(theta1, N1);
  y0 ~ dbinom(theta0, N0);

  # priors
  beta ~ dnorm(0, 1/(sigma.beta^2));
  
  theta0 <- exp(beta)/(1+exp(beta));
  theta1 <- exp(beta)/(1+exp(beta));
}"
  data <- list(y1=y1, y0=y0, N1=N1, N0=N0,
               sigma.beta=sigma.beta)
  model <- jags.model(file = textConnection(model_code),data = data)
  update(model, n.iter = n.iter)
  samps <- coda.samples(model = model, 
                        variable.names = c("beta","theta0",'theta1'),
                        n.iter = n.iter)
  samps <- data.frame(samps[[1]])
  return(samps)
}

# calculate LT unnormalized log posterior (for THAMES)
lp_lt <- function(y0,y1,N0,N1,beta,psi,sigma.beta=1,sigma.psi=1){
  theta0 <- exp(beta-psi/2)/(1+exp(beta-psi/2))
  theta1 <- exp(beta+psi/2)/(1+exp(beta+psi/2))
  return(
    dbinom(y0,N0,theta0,log=TRUE) +
      dbinom(y1,N1,theta1,log=TRUE) +
      dnorm(beta,0,sigma.beta,log=TRUE) +
      dnorm(psi,0,sigma.psi,log=TRUE)
    )
}

lp_lt_h0 <- function(y0,y1,N0,N1,beta,sigma.beta=1){
  theta0 <- exp(beta)/(1+exp(beta))
  theta1 <- exp(beta)/(1+exp(beta))
  return(
    dbinom(y0,N0,theta0,log=TRUE) +
      dbinom(y1,N1,theta1,log=TRUE) +
      dnorm(beta,0,sigma.beta,log=TRUE)
  )
}


