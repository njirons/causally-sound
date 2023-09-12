#' THAMES estimator of the (reciprocal) log marginal likelihood
#'
#' This function computes the THAMES estimate of the reciprocal
#'     log marginal likelihood using posterior samples and
#'     unnormalized log posterior values.
#'
#' @param lps vector of unnormalized log posterior values of length n_samples
#' @param params matrix of parameter posterior samples of dimension n_samples * d
#' @param n_samples integer, number of posterior samples
#' @param d integer, dimension of parameter space
#' @param radius positive number, radius to use for defining the ellipsoid A
#' @param p percentile, used for lower bound of confidence interval
#' @param q percentile, used for upper bound of confidence interval
#' @param lp_func function to compute unnormalized log posterior values
#' @param bound function calculating membership of a point in the posterior support
#' @param n_simuls integer, number of Monte Carlo simulations
#'     to use in the bounded parameter correction calculation.
#'
#' @return Returns a named list with the following elements:
#'
#' @examples
#'
#' @seealso
#'
#' @references Metodiev M, Perrot-Dockès M, Ouadah S, Irons N. J., Raftery A. E.
#'     (2023) Easily Computed Marginal Likelihoods from Posterior Simulation
#'     Using the THAMES Estimator. arXiv preprint.
#'
#' @export
thames <- function(lps = NULL,params,n_samples = NULL,d = NULL, radius = NULL,
                   p = 0.025, q = 1-p, lp_func = NULL,
                   bound = NULL, n_simuls = 1e5){

  # dimension of parameter space
  if(is.null(d)){
    d <- dim(params)[2]
    if(is.null(d)){
      d <- 1
    }
  }

  # radius of A
  if(is.null(radius)){
    radius <- sqrt(d+1)
  }

  # number of posterior samples
  if(is.null(n_samples)){
    if(d==1){
      n_samples <- length(params)
    }else{
      n_samples <- dim(params)[1]
    }
  }

  # calculate unnormalized log posterior values
  if(is.null(lps)){
    lps <- lp_func(params)
  }

  if(length(lps) != n_samples){
    return('Error: number of unnormalized log posterior values does not match posterior sample size.')
  }

  # split the sample
  n1 <- n_samples %/% 2
  n2 <- n_samples - n1
  if(d==1){
    params1 <- params[1:n1]
    params2 <- params[(n1+1):n_samples]
  }else{
    params1 <- params[1:n1,]
    params2 <- params[(n1+1):n_samples,]
  }
  lps1 <- lps[1:n1]
  lps2 <- lps[(n1+1):n_samples]

  # calculate posterior mean and covariance from first half of sample
  if(d==1){
    theta_hat <- params1[which.max(lps1)]
    sigma_hat <- var(params1)
    
    log_det_sigma_hat <- log(sigma_hat)
    
    # which samples are in A?
    inA <- sapply(params2,function(theta){
      # calculate distance from theta_hat
      theta_tilde <- (theta-theta_hat)/sqrt(sigma_hat)
      
      # is distance of theta less than the radius?
      return(sum(theta_tilde^2) < radius^2)
    })
  }else{
    # theta_hat <- colMeans(params1)
    theta_hat <- params1[which.max(lps1),]
    sigma_hat <- cov(params1)    
    
    # calculate SVD of sigma_hat
    sigma_svd = svd(sigma_hat)
    
    # calculate log(det(sigma_hat))
    log_det_sigma_hat = sum(log(sigma_svd$d))
    
    # which samples are in A?
    inA <- apply(params2,1,function(theta){
      # calculate distance from theta_hat
      theta_tilde <-  sigma_svd$d^(-1/2) * (t(sigma_svd$v) %*% (theta-theta_hat))
      
      # is distance of theta less than the radius?
      return(sum(theta_tilde^2) < radius^2)
    })
  }

  # log volume of A
  logvolA = d*log(radius)+(d/2)*log(pi)+log_det_sigma_hat/2-lgamma(d/2+1)

  # calculate log(zhat_inv)
  log_zhat_inv  = log(mean(exp(-(lps2-max(lps2)))*inA))-logvolA-max(lps2)

  # calculate bounded parameter correction (if necessary)
  r_bound <- 1
  if(!is.null(bound)){
    r_bound <- bound_par_cor(theta_hat, sigma_svd, bound, radius, n_simuls)
    log_zhat_inv <- log_zhat_inv - log(r_bound)
  }

  # estimate ar(1) model for lps
  lp_ar <- ar(exp(-(lps2-max(lps2)))*inA, order.max=1)
  phi <- lp_ar$partialacf[1]

  # correct standard error for autocorrelation in MCMC samples
  standard_error <- sd(exp(-lps2+max(lps2))*inA)/
    ((1-phi)*sqrt(n2)*mean(exp(-lps2+max(lps2))*inA))

  # lower bound of confidence interval
  log_zhat_inv_L <- log_zhat_inv + log(trunc_quantile(p,standard_error))

  # upper bound of confidence interval
  log_zhat_inv_U <- log_zhat_inv + log(trunc_quantile(q,standard_error))

  return(list(log_zhat_inv = log_zhat_inv,
              log_zhat_inv_L = log_zhat_inv_L,
              log_zhat_inv_U = log_zhat_inv_U,
              theta_hat = theta_hat,
              sigma_hat = sigma_hat,
              # sigma_svd = sigma_svd,
              log_det_sigma_hat = log_det_sigma_hat,
              logvolA = logvolA, inA = inA, r_bound = r_bound,
              se = standard_error, phi = phi, lp_ar = lp_ar,
              lps = lps, params = params, n_samples = n_samples,
              d = d, radius = radius, p = p, q = q, lp_func = lp_func,
              bound = bound, n_simuls = n_simuls
  ))
}
