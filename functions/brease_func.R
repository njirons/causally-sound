
log_prior_brease <- function(theta0, etae, etas, 
                          mu0 = 1/2, n0 = 2,
                          mue = 1/2, ne = mu0*n0,
                          mus = 1/2, ns = (1-mu0)*n0){
  dbeta(theta0, mu0*n0, (1-mu0)*n0, log = T) + 
    dbeta(etae, mue*ne, (1-mue)*ne, log = T) + 
    dbeta(etas, mus*ns, (1-mus)*ns, log = T) 
}

f_theta1 <- function(theta0, etae, etas){
  theta0*(1-etae) + (1-theta0)*etas
}

# simulates from the prior
sim_brease <- function(n.iter = 1e4, 
                    mu0=1/2, n0=2, 
                    mue = 1/2, ne  = mu0*n0,
                    mus = 1/2, ns = (1-mu0)*n0){
  out <- list()
  out$theta0 <- rbeta(n.iter,   mu0*n0, (1-mu0)*n0)
  out$etas   <- rbeta(n.iter,   mus*ns, (1-mus)*ns)
  out$etae   <- rbeta(n.iter,   mue*ne, (1-mue)*ne)
  out$theta1 <- with(out, theta0*(1-etae) + (1-theta0)*etas)
  out <- as.data.frame(out)
  return(out)
}

sim_IB <- function(n.iter = 1e4, 
                   mu0=1/2, n0=2, 
                   mu1 = mu0, n1= n0){
  out <- list()
  out$theta0 <- rbeta(n.iter,   mu0*n0, (1-mu0)*n0)
  out$theta1 <- rbeta(n.iter,   mu1*n1, (1-mu1)*n1)
  out <- as.data.frame(out)
  return(out)
}


sim_lt <- function(n.iter = 1e4, 
                   mu_psi = 0, var_psi  = 1, 
                   mu_beta = 0, var_beta = 1){
  out    <- list()
  out$psi    <- rnorm(n.iter, mean = mu_psi, sd = sqrt(var_psi))
  out$beta   <- rnorm(n.iter, mean = mu_beta, sd = sqrt(var_beta))
  out$theta0 <- with(out, exp(beta - psi/2)/(1 + exp(beta - psi/2)))
  out$theta1 <- with(out, exp(beta + psi/2)/(1 + exp(beta + psi/2)))
  out <- as.data.frame(out)
  return(out)
}

joint_density <- function(theta0, theta1, h, n=300, ...){
  dens <- kde2d(theta0, theta1, h, n = n)
  filled.contour(dens, ...)
}


# brease log marginal likelihood (null hypothesis)
lml0_brease <- function(y0,y1,N0,N1,mu0=1/2,n0=2){
  a0 <- mu0*n0
  b0 <- (1-mu0)*n0
  return(lml0_IB(y0,y1,N0,N1,a0,b0))
}

# brease log marginal likelihood (alternative hypothesis)
lml1_brease <- function(y0,y1,N0,N1,
                     mu0=1/2,n0=2,
                     mue=1/2,ne=2,
                     mus=1/2,ns=2){
  N <- N0+N1
  
  lml <- lchoose(N0, y0) + lchoose(N1, y1) -
    lbeta(mu0*n0,(1-mu0)*n0) - lbeta(mue*ne,(1-mue)*ne) - lbeta(mus*ns,(1-mus)*ns)
  
  J <- (0:y1) %o% rep(1,N1-y1+1)
  K <- rep(1,y1+1) %o% (0:(N1-y1))
  
  summands <- lchoose(y1,J) + lchoose(N1-y1,K) +
    lbeta(y0 + J + K + mu0*n0, N - (y0 + J + K) + (1-mu0)*n0) +
    lbeta(K + mue*ne, J + (1-mue)*ne) +
    lbeta(y1 - J + mus*ns, N1-y1 - K + (1-mus)*ns)
  
  max_log <- max(summands)
  lml <- lml + log(sum(exp(summands-max_log))) + max_log
  return(lml)
}

lml_H0_brease <- function(y0,y1,N0,N1,mu0=1/2,n0=2,mue=1/2,mus=1/2){
  ne <- mu0*n0
  ns <- (1-mu0)*n0
  
  lml <- lchoose(N0, y0) + lchoose(N1, y1) - 
    lgamma((1-mus)*ns) - lgamma(mue*ne+mus*ns) - lgamma((1-mue)*ne)
  
  J <- (0:(y0+y1)) %o% rep(1,N0+N1-y0-y1+1)
  K <- rep(1,y0+y1+1) %o% (0:(N0+N1-y0-y1))
  
  a00 <- N0+N1-y0-y1+(1-mu0)*n0-mus*ns-K
  a10 <- J+K+mue*ne+mus*ns
  a11 <- y0+y1+mu0*n0-mue*ne-J

  summands <- 
    -log(2)*(J+K) +
    lchoose(y0+y1,J) +
    lchoose(N0+N1-y0-y1,K) +
    lgamma(a00) + lgamma(a10) +lgamma(a11) - lgamma(a00+a10+a11)
  
  max_log <- max(summands)
  lml <- lml + log(sum(exp(summands-max_log))) + max_log
  return(lml)
}

# brease log Bayes factor
bf_brease <- function(y0,y1,N0,N1,                     
                   mu0=1/2,n0=2,
                   mue=1/2,ne=2,
                   mus=1/2,ns=2){
  return(lml0_brease(y0,y1,N0,N1,mu0,n0) - 
           lml1_brease(y0,y1,N0,N1,mu0,n0,mue,ne,mus,ns))
}

# monotonic brease log marginal likelihood (alternative hypothesis)
lml1_brease_mono <- function(y0,y1,N0,N1,
                     mu0=1/2,n0=2,
                     mue=1/2,ne=2){
  N <- N0+N1
  
  lml <- lchoose(N0, y0) + lchoose(N1, y1) - 
    lbeta(mu0*n0,(1-mu0)*n0) - lbeta(mue*ne,(1-mue)*ne)
  
  # summands <- rep(NA,N1-y1+1)
  summands <- 
    lchoose(N1-y1, 0:(N1-y1)) +
    lbeta(y0+y1+0:(N1-y1)+mu0*n0,N-(y0+y1+0:(N1-y1))+(1-mu0)*n0) +
    lbeta(0:(N1-y1)+mue*ne,y1+(1-mue)*ne)
  
  max_log <- max(summands)
  lml <- lml + log(sum(exp(summands-max_log))) + max_log
  return(lml)
}

# monotonic brease log Bayes factor
bf_brease_mono <- function(y0,y1,N0,N1,                     
                   mu0=1/2,n0=2,
                   mue=1/2,ne=2){
  return(lml0_brease(y0,y1,N0,N1,mu0,n0) - 
           lml1_brease_mono(y0,y1,N0,N1,mu0,n0,mue,ne))
}

sample_brease <- function(n_sample,y0,y1,N0,N1,
                       mu0=1/2,n0=2,
                       mue=1/2,ne=2,
                       mus=1/2,ns=2,
                       log_weights){
  N <- N0+N1
  
  log_weights <- (outer(lchoose(y1,0:y1),lchoose(N1-y1,0:(N1-y1)),'+')) +
    lbeta(y1+y0+mu0*n0+outer(-(0:y1),0:(N1-y1),"+"),
          N-y0-y1+(1-mu0)*n0+outer(0:y1,-(0:(N1-y1)),'+')) +
    lbeta(rep(1,y1+1) %o% 0:(N1-y1)+mue*ne,
          y1+(1-mue)*ne-(0:y1) %o% rep(1,N1-y1+1)) +
    lbeta(0:y1 %o% rep(1,N1-y1+1) + mus*ns, 
          N1-y1-rep(1,y1+1) %o% 0:(N1-y1) + (1-mus)*ns)
  log_weights <- log_weights - max(log_weights)
  weights <- exp(log_weights)
  rm(log_weights)
  
  # sample x10
  k <- sample(x = 0:(N1-y1),
                size = n_sample,
                prob = colSums(weights),
                replace = TRUE
  )
  j <- sapply(1:n_sample, function(n){
    sample(x = 0:y1,size = 1,prob = weights[,k[n]+1])})
  
  # sample parameters
  etae <- rbeta(n_sample, k + mue*ne, y1 - j + (1-mue)*ne)
  etas <- rbeta(n_sample, j + mus*ns, N1 - y1 - k + (1-mus)*ns)
  theta0 <- rbeta(n_sample, y1 + y0 - j + k + mu0*n0, 
                  N - (y0 + y1 - j + k) + (1-mu0)*n0)
  theta1 <- theta0*(1-etae)+(1-theta0)*etas
  
  return(list(
    j = j, k = k, 
    weights = weights,
    theta0 = theta0,
    etae = etae,
    etas = etas,
    theta1 = theta1))
}

# monotonic brease posterior samples (no harm)
sample_brease_noharm <- function(n_sample,y0,y1,N0,N1,
                       mu0=1/2,n0=2,
                       mue=1/2,ne=2){
  N <- N0+N1
  
  log_weights <- lchoose(N1-y1,0:(N1-y1)) +
    lbeta(y0+y1+mu0*n0+0:(N1-y1),
          N-y0-y1+(1-mu0)*n0-(0:(N1-y1))) +
    lbeta(0:(N1-y1)+mue*ne,
          y1+(1-mue)*ne)
  log_weights <- log_weights - max(log_weights)
  weights <- exp(log_weights)
  rm(log_weights)
  
  # sample x11
  k <- sample(x = 0:(N1-y1),
              size = n_sample,
              prob = weights,
              replace = TRUE
  )
  
  # sample parameters
  etae <- rbeta(n_sample, k + mue*ne, y1 + (1-mue)*ne)
  theta0 <- rbeta(n_sample, y1 + y0 + k + mu0*n0, 
                  N - (y0 + y1 + k) + (1-mu0)*n0)
  theta1 <- theta0*(1-etae)
  
  return(list(
    k = k, 
    weights = weights,
    theta0 = theta0,
    etae = etae,
    theta1 = theta1))
}

# brease log posterior (unnormalized) direct calculation
lp_brease <- function(samples,data){
  y0 <- data[1]; y1 <- data[2]; N0 <- data[3]; N1 <- data[4]; N <- N0+N1;
  mu0 <- data[5]; n0 <- data[6]; 
  mue <- data[7]; ne <- data[8];
  mus <- data[9]; ns <- data[10]
  
  theta0 <- samples[['theta0']]
  eta_e <- samples[['etae']]
  eta_s <- samples[['etas']]
  
  theta1 <- theta0*(1-eta_e)+(1-theta0)*eta_s
  
  # likelihood terms
  lp <- lchoose(N0, y0) + lchoose(N1, y1) +
    y0*log(theta0) + (N0-y0)*log(1-theta0) +
    y1*log(theta1) + (N1-y1)*log(1-theta1)
  
  # prior terms
  lp <- lp + dbeta(theta0,mu0*n0,(1-mu0)*n0,log=TRUE) +
    dbeta(eta_e,mue*ne,(1-mue)*ne,log=TRUE) +
    dbeta(eta_s,mus*ns,(1-mus)*ns,log=TRUE)
  
  return(lp)
}

lp_theta0_brease <- function(theta0,y0,y1,N0,N1,
                          mu0=1/2,n0=2,
                          mue=1/2,ne=2,
                          mus=1/2,ns=2){
  y10 <- 0:y1
  x11 <- 0:(N1-y1)
  N <- N0+N1
  Y10 <- (y10 %o% rep(1,N1-y1+1))
  X11 <- (rep(1,y1+1) %o% x11)
  
  lps <- outer(lchoose(y1,y10),lchoose(N1-y1,x11),'+') +
    lbeta(y0+y1-Y10+X11+mu0*n0, N-(y0+y1-Y10+X11)+(1-mu0)*n0) +
    lbeta(X11+mue*ne, y1-Y10+(1-mue)*ne) +
    lbeta(Y10+mus*ns, N1-y1-X11+(1-mus)*ns) +
    dbeta(x = theta0,
          shape1 = y0+y1-Y10+X11+mu0*n0,
          shape2 = N-(y0+y1-Y10+X11)+(1-mu0)*n0,
          log=TRUE)
  
  return(log(mean(exp(lps-max(lps))))+max(lps))
}

# calculate appell's f1 function (log scale)
log_f1 <- function(a,b1,b2,c,x,y,subdiv=100,rel.tol=1e-10){
  log(integrate(
    f = function(u){
      (u^(a-1))*((1-u)^(c-a-1))*((1-u*x)^(-b1))*((1-u*y)^(-b2))
    },
    lower=0,
    upper=1,
    subdivisions = subdiv,
    rel.tol=rel.tol
    )$value) - lbeta(a,c-a)
}

# Calculate brease conditional log prior of theta1 | theta0
library(tolerance)
logprior_10_brease <- function(theta1,theta0,
                            mu0=1/2,n0=2,
                            mue=1/2,ne=2,
                            mus=1/2,ns=2,
                            subdiv = 100,
                            rel.tol = 1e-10){
  if(theta0 <= 0.5){
    if(theta1 <= theta0){
      lp <- ((1-mue)*ne+mus*ns-1)*log(theta1) + 
        (mue*ne-1)*log(theta0-theta1) + 
        lbeta(mus*ns,(1-mue)*ne) -
        (ne-1)*log(theta0) -
        mus*ns*log(1-theta0) -
        lbeta(mus*ns,(1-mus)*ns) -
        lbeta((1-mue)*ne,mue*ne) +
        log_f1(
          mus*ns,
          1-mue*ne,
          1-(1-mus)*ns,
          (1-mue)*ne+mus*ns,
          -theta1/(theta0-theta1),
          theta1/(1-theta0),
          subdiv = subdiv,
          rel.tol=rel.tol
        )
    }else if(theta1 <= 1-theta0){
      lp <- (mus*ns-1)*log(theta1-theta0) +
        ((1-mus)*ns-1)*log(1-theta1) -
        (ns-1)*log(1-theta0) -
        lbeta(mus*ns,(1-mus)*ns) +
        log_f1(
          mue*ne,
          1-mus*ns,
          1-(1-mus)*ns,
          ne,
          -theta0/(theta1-theta0),
          theta0/(1-theta1),
          subdiv = subdiv,
          rel.tol=rel.tol
        )
    }else{
      lp <- (mue*ne+(1-mus)*ns-1)*log(1-theta1) +
        (mus*ns-1)*log(theta1-theta0) +
        lbeta(mue*ne,(1-mus)*ns) - 
        mue*ne*log(theta0) - 
        (ns-1)*log(1-theta0) - 
        lbeta(mus*ns,(1-mus)*ns) -
        lbeta((1-mue)*ne,mue*ne) +
        log_f1(
          mue*ne,
          1-(1-mue)*ne,
          1-mus*ns,
          mue*ne+(1-mus)*ns,
          (1-theta1)/theta0,
          (theta1-1)/(theta1-theta0),
          subdiv = subdiv,
          rel.tol=rel.tol
        )
    }
    
  }else{
    if(theta1 <= 1-theta0){
      lp <- 0
    }else if(theta1 <= theta0){
      lp <- 0
    }else{
      lp <- 0
    }
  }
  return(lp)
}

# integrate joint prior to get marginal for theta1
logprior_theta1_brease <- function(theta1,
                                mu0=1/2,n0=2,
                                mue=1/2,ne=2,
                                mus=1/2,ns=2){
  log(integrate(
    f = function(theta0){
      exp(
        sapply(theta0,function(t0){logprior_10_brease(theta1,t0,mu0,n0,mue,ne,mus,ns)})+
        dbeta(theta0,mu0*n0,(1-mu0)*n0,log=TRUE)
        )
    },
    lower=0, upper = 0.5)$value)
}

# calculate marginal posterior of theta1
logpost_theta1_brease <- function(theta1,y0,y1,N0,N1,
                                mu0=1/2,n0=2,
                                mue=1/2,ne=2,
                                mus=1/2,ns=2,
                               subdiv = 100,
                               rel.tol=1e-10){
  log(integrate(
    f = function(theta0){
      exp(
        sapply(theta0,
               function(t0){logprior_10_brease(
                 theta1,t0,mu0,n0,mue,ne,mus,ns,subdiv=subdiv,rel.tol=rel.tol)})+
          dbeta(theta0,mu0*n0,(1-mu0)*n0,log=TRUE) +
          lchoose(N0,y0) + lchoose(N1,y1) + 
          y0*log(theta0)+(N0-y0)*log(1-theta0) +
          y1*log(theta1)+(N1-y1)*log(1-theta1)
        )
    },
    lower=0, upper = 0.5, subdivisions=subdiv,rel.tol=rel.tol)$value)
}

# calculate bounds on probabilities of causation
brease_bounds <- function(y0,y1,N0,N1){
  
  # add 1 success and 1 failure (e.g., based on uniform prior)
  # for numerical stability near 0/1
  theta0_hat <- (y0+1)/(N0+2)
  theta1_hat <- (y1+1)/(N1+2)
    
  # bounds based on plug-in estimators
  if(theta0_hat == 1){ #etas not partially identified
    etas_l <- 0
    etas_u <- 1
  }else{
    etas_l <- max(0,(theta1_hat-theta0_hat)/(1-theta0_hat))
    etas_u <- min(1,theta1_hat/(1-theta0_hat))
  }
  etas_m <- (etas_l+etas_u)/2 
  
  if(theta0_hat == 0){ #etae not partially identified
    etae_l <- 0
    etae_u <- 1
  }else{
    etae_l <- max(0,(theta0_hat-theta1_hat)/theta0_hat)
    etae_u <- min(1,(1-theta1_hat)/theta0_hat)    
  }
  etae_m <- (etae_l+etae_u)/2  
  
  return(list(
    etas_l = etas_l, etas_u = etas_u, etas_m = etas_m,
    etae_l = etae_l, etae_u = etae_u, etae_m = etae_m
  ))
}

# brease Stan model
brease_stan_code <- 'data {
  int<lower=1> N1; // number treated
  int<lower=1> N0; // number control
  int<lower=0,upper=N1> y1; // outcome, treated 
  int<lower=0,upper=N0> y0; // outcome, control 
  
  real<lower=0,upper=1> mu0; // prior mean
  real<lower=0,upper=1> mue;
  real<lower=0,upper=1> mus;
  real<lower=0> n0; //prior sample size
  real<lower=0> ne;
  real<lower=0> ns;
}
parameters {
  real<lower=0,upper=1> theta0; // baseline risk
  real<lower=0,upper=1> etae; // efficacy  
  real<lower=0,upper=1> etas; // risk of side effects
}
transformed parameters {
  real<lower=0,upper=1> theta1 = (1-etae)*theta0 + etas*(1-theta0);
  // risk in treated group
}
model {
  // priors
  target += beta_lpdf(theta0 | mu0*n0, (1-mu0)*n0);
  target += beta_lpdf(etae | mue*ne, (1-mue)*ne);
  target += beta_lpdf(etas | mus*ns, (1-mus)*ns);
  
  // (log) likelihood
  target += binomial_lpmf(y0 | N0, theta0);
  target += binomial_lpmf(y1 | N1, theta1);
}
'

# Null model Stan code
null_stan_code <- 'data {
  int<lower=1> N1; // number treated
  int<lower=1> N0; // number control
  int<lower=0,upper=N1> y1; // outcome, treated 
  int<lower=0,upper=N0> y0; // outcome, control 
  
  real<lower=0,upper=1> mu0; // prior mean
  real<lower=0> n0; //prior sample size
}
parameters {
  real<lower=0,upper=1> theta0; // baseline risk
}
model {
  // priors
  target += beta_lpdf(theta0 | mu0*n0, (1-mu0)*n0);
  
  // (log) likelihood
  target += binomial_lpmf(y0 | N0, theta0);
  target += binomial_lpmf(y1 | N1, theta0);
}
'

# Neyman's null model Stan code
neyman_stan_code <- 'data {
  int<lower=1> N1; // number treated
  int<lower=1> N0; // number control
  int<lower=0,upper=N1> y1; // outcome, treated 
  int<lower=0,upper=N0> y0; // outcome, control 
  
  vector[3] a; // Dirichlet parameters
}
parameters {
  simplex[3] p;
}
model {
  // priors
  target += dirichlet_lpdf(p | a);
  
  // (log) likelihood
  target += binomial_lpmf(y0 | N0, p[2]/2+p[3]);
  target += binomial_lpmf(y1 | N1, p[2]/2+p[3]);
}
'


# Mono Stan model
mono_stan_code <- 'data {
  int<lower=1> N1; // number treated
  int<lower=1> N0; // number control
  int<lower=0,upper=N1> y1; // outcome, treated 
  int<lower=0,upper=N0> y0; // outcome, control 
  
  real<lower=0,upper=1> mu0; // prior mean
  real<lower=0,upper=1> mue;
  real<lower=0,upper=1> mus;
  real<lower=0> n0; //prior sample size
  real<lower=0> ne;
  real<lower=0> ns;
}
parameters {
  real<lower=0,upper=1> theta0; // baseline risk
  real<lower=0,upper=1> etae; 
}
transformed parameters {
  real<lower=0,upper=1> theta1 = (1-etae)*theta0;
  // risk in treated group
}
model {
  // priors
  target += beta_lpdf(theta0 | mu0*n0, (1-mu0)*n0);
  target += beta_lpdf(etae | mue*ne, (1-mue)*ne);
  
  // (log) likelihood
  target += binomial_lpmf(y0 | N0, theta0);
  target += binomial_lpmf(y1 | N1, theta1);
}
'
# Mono Stan model with Jeffreys prior
mono_jeffreys_stan_code <- 'data {
  int<lower=1> N1; // number treated
  int<lower=1> N0; // number control
  int<lower=0,upper=N1> y1; // outcome, treated 
  int<lower=0,upper=N0> y0; // outcome, control 
}
parameters {
  real<lower=0,upper=1> theta0; // baseline risk
  real<lower=0,upper=1> etae; 
}
transformed parameters {
  real<lower=0,upper=1> theta1 = (1-etae)*theta0;
  // risk in treated group
}
model {
  // priors
  target += (-0.5)*(log(1-theta0)+log(1-etae)+log(1-theta0*(1-etae)));
  
  // (log) likelihood
  target += binomial_lpmf(y0 | N0, theta0);
  target += binomial_lpmf(y1 | N1, theta1);
}
'

stan.brease <- function(y0, y1, N0, N1, 
                     mu0=1/2, n0=2, 
                     mue=1/2, ne =2, 
                     mus = 1/2, ns = 2,
                     iter = 1e4, 
                     code = "h1"){
  brease_data <- list(
    N1 = N1,
    y1 = y1,
    N0 = N0, 
    y0 = y0,
    mu0 = mu0, n0 = n0,
    mue = mue, ne = ne,
    mus = mus, ns = ns)
  
  code <- switch (code,
    h1 = brease_stan_code,
    mono = mono_stan_code)
  
  brease_model <- stan_model(model_code = code)
  brease_fit <- sampling(brease_model, 
                       data = brease_data,
                       iter = iter)
  samples <- extract(brease_fit)
  out <- list(model = brease_fit, samples = samples)
  return(out)
}

# sampling under Neyman's null
stan.brease.H0 <- function(y0, y1, N0, N1, 
                     mu0=1/2, n0=2, 
                     mue=1/2,
                     mus = 1/2,
                     iter = 1e4){
  ne <- mu0*n0
  ns <- (1-mu0)*n0
  
  a <- c(
    (1-mus)*ns,
    mue*ne+mus*ns,
    (1-mue)*ne
  )
  brease_data <- list(
    N1 = N1,
    y1 = y1,
    N0 = N0, 
    y0 = y0,
    a = a
    )
  
  code <- neyman_stan_code
  
  brease_model <- stan_model(model_code = code)
  brease_fit <- sampling(brease_model, 
                      data = brease_data,
                      iter = iter)
  samples <- extract(brease_fit)
  out <- list(model = brease_fit, samples = samples)
  return(out)
}

# jags brease samples
jags.brease <- function(y0,y1,N0,N1, 
                     mu0=1/2, n0=2,
                     mue=1/2, ne=2,mus=1/2,ns=1/2, n.iter = 1e4){
  
  a0 <- mu0*n0
  b0 <- (1-mu0)*n0
  ae <- ne*mue
  be <- ne*(1-mue)
  as <- ns*mus
  bs <- ns*(1-mus)
  
  
  model_code <- "model{
  # likelihood
  y1 ~ dbinom(theta1, N1);
  y0 ~ dbinom(theta0, N0);

  # priors
  theta0  ~ dbeta(a0,b0);
  etae    ~ dbeta(ae,be);
  etas    ~ dbeta(as,bs);
  theta1 <- theta0*(1-etae)+etas*(1-theta0);
}"
  data <- list(y1=y1, y0=y0, N1=N1, N0=N0,a0=a0, b0=b0, ae=ae, be = be, as=as,bs=bs)
  model <- jags.model(file = textConnection(model_code),data = data)
  update(model, n.iter = n.iter)
  samps <- coda.samples(model = model, 
                        variable.names = c("theta1", "theta0",  "etae",'etas'),
                        n.iter = n.iter)
  samps <- data.frame(samps[[1]])
  return(samps)
}

# Plus-minus BF -----------------------------------------------------------

avg_mono <- function(y0, y1, N0, N1, mue, ne){
  (exp(lml1_brease_mono(y0, y1, N0, N1, mue = mue, ne = ne)) + 
     exp(lml1_brease_mono(N1-y1, N0-y0, N1, N0, mue = mue, ne = ne)))/2
}
avg_side <- function(y0, y1, N0, N1, mue, ne){
  (exp(lml1_brease_mono(N0-y0, N1-y1, N0, N1,mue = mue, ne = ne)) + 
     exp(lml1_brease_mono(y1, y0, N1, N0, mue = mue, ne = ne)))/2
}

bf_side <- function(y0, y1, N0, N1, mue, ne){
  exp(lml0_brease(y0,y1, N0, N1))/avg_side(y0, y1, N0, N1, mue = mue, ne =ne)
}
bf_mono <- function(y0, y1, N0, N1,mue, ne){
  exp(lml0_brease(y0,y1, N0, N1))/avg_mono(y0, y1, N0, N1, mue = mue, ne =ne)
}
