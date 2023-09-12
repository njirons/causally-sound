#' Truncated normal quantiles
#'
#' This function computes quantiles of the truncated normal distribution
#'    for calculating THAMES confidence intervals.
#'
#' @param p Percentile
#' @param ratio Ratio of standard error to point estimate
#'     (midpoint of confidence interval)
#'
#' @return Truncated normal quantile
#'
#' @keywords internal
trunc_quantile = function(p,ratio){
  alpha = - 1/ratio
  (qnorm(p=pnorm(alpha)+p*(1-pnorm(alpha))) * (ratio))+1
}
