#' Get Continuous CDF
#' 
#' Internal dispatcher to calculate the Cumulative Distribution Function (CDF)
#' based on the numeric family code from the TMB object.
#' 
#' @param family_code Integer representing the distribution family.
#' @param Y Vector of observed responses.
#' @param mu Vector of location parameters.
#' @param phi Vector of dispersion/shape parameters.
#' @param tau Quantile level.
#' @return A numeric vector of CDF values.
#' @noRd
get_cdf_continuous <- function(family_code, Y, mu, phi, tau) {
  switch(as.character(family_code),
         
         "0" = pbeta_mean(Y, mu, phi),
         "1" = psimplex(Y, mu, phi), 
         "2" = pvasicek_mean(Y, mu, phi),         
         "3" = pugamma_mean(Y, mu, phi),
         "4" = pbessel(Y, mu, phi),
         
         "5" = pbeta_mode(Y, mu, phi),
         "6" = pkum_mode(Y, mu, phi),
         "7" = pugamma_mode(Y, mu, phi),
         "8" = pugompertz_mode(Y, mu, phi),
                  
         "9"  = pkum_quantile(Y, mu, phi, tau),
         "10" = pvasicek_quantile(Y, mu, phi, tau),         
         "11" = puweibull(Y, mu, phi, tau),
         "12" = pugompertz_quantile(Y, mu, phi, tau),
         "13" = pjohnsonsb(Y, mu, phi, tau),
         "14" = pashw(Y, mu, phi, tau),
         "15" = pubs(Y, mu, phi, tau), 
         
         stop(paste("CDF for family code", family_code, "is not implemented."))
  )
}

#' Get Random Continuous Deviates
#' 
#' Internal dispatcher to generate random data from the continuous part of the model.
#' 
#' @param family_code Integer representing the distribution family.
#' @param n Number of observations to generate.
#' @param mu Vector of location parameters.
#' @param phi Vector of dispersion/shape parameters.
#' @param tau Quantile level.
#' @return A numeric vector of simulated continuous responses.
#' @noRd
get_random_continuous <- function(family_code, n, mu, phi, tau) {
  switch(as.character(family_code),
         
         "0" = rbeta_mean(n, mu, phi),
         "1" = rsimplex(n, mu, phi),          
         "2" = rvasicek_mean(n, mu, phi),
         "3" = rugamma_mean(n, mu, phi),
         "4" = rbessel(n, mu, phi),
       
         "5" = rbeta_mode(n, mu, phi),
         "6" = rkum_mode(n, mu, phi),         
         "7" = rugamma_mode(n, mu, phi),
         "8" = rugompertz_mode(n, mu, phi),

         "9"  = rkum_quantile(n, mu, phi, tau),
         "10" = rvasicek_quantile(n, mu, phi, tau),  
         "11" = ruweibull(n, mu, phi, tau),
         "12" = rugompertz_quantile(n, mu, phi, tau),
         "13" = rjohnsonsb(n, mu, phi, tau),
         "14" = rashw(n, mu, phi, tau),
         "15" = rubs(n, mu, phi, tau),
         
         stop(paste("Random data generation for family code", family_code, "is not implemented."))
  )
}

#' Correct Continuous CDF for Zero-One Inflation
#' 
#' Internal function to apply Dunn & Smyth (1996) randomization for discrete bounds.
#' @noRd
correct_for_inflation <- function(Y, cdf_cont, p0, p1) {
  n <- length(Y)
  cdf_val <- numeric(n)
  
  if(length(p0) == 1) p0 <- rep(p0, n)
  if(length(p1) == 1) p1 <- rep(p1, n)
  
  p0[is.na(p0)] <- 0
  p1[is.na(p1)] <- 0
  
  sum_p <- p0 + p1
  idx_overflow <- which(sum_p >= 1)
  if (length(idx_overflow) > 0) {
    p0[idx_overflow] <- p0[idx_overflow] / (sum_p[idx_overflow] + 1e-6)
    p1[idx_overflow] <- p1[idx_overflow] / (sum_p[idx_overflow] + 1e-6)
  }
  
  p_cont <- 1 - p0 - p1
  
  idx0 <- which(Y == 0)
  idx1 <- which(Y == 1)
  idx_cont <- which(Y > 0 & Y < 1)
  
  ## 1. Continuous values: Exact mixture (no randomization)
  if (length(idx_cont) > 0) {
    cdf_val[idx_cont] <- p0[idx_cont] + p_cont[idx_cont] * cdf_cont[idx_cont]
  }
  
  ## 2. Exact zeros: Uniform randomization on the step [0, p0]
  if (length(idx0) > 0) {
    cdf_val[idx0] <- stats::runif(length(idx0), min = 0, max = p0[idx0])
  }
  
  ## 3. Exact ones: Uniform randomization on the step [1 - p1, 1]
  if (length(idx1) > 0) {
    cdf_val[idx1] <- stats::runif(length(idx1), min = 1 - p1[idx1], max = 1)
  }
  
  return(cdf_val)
}

#' Compute Residuals from CDF
#' @noRd
compute_residuals <- function(cdf_val, type) {
  
  eps <- 1e-7
  cdf_val <- pmin(pmax(cdf_val, eps), 1 - eps)
  
  if (type == "quantile") {
    return(stats::qnorm(cdf_val))
  } else if (type == "cox-snell") {
    return(-log(1 - cdf_val))
  } else {
    stop("Unknown residual type.")
  }
}

#' @title Extract Model Residuals for unitregTMB
#' 
#' @description 
#' Computes randomized quantile or Cox-Snell residuals for continuous bounded 
#' models, accounting for potential zero and one inflation (Dunn & Smyth, 1996).
#' 
#' @param object A fitted \code{unitregTMB} model.
#' @param type The type of residuals which should be returned. 
#'        Options are \code{"quantile"} (default), \code{"cox-snell"}, or \code{"response"}.
#' @param ... Additional arguments (currently unused).
#' 
#' @return A numeric vector representing the chosen residuals.
#' 
#' @examples
#' \donttest{
#' # Assuming 'da' is your dataset with a zero-inflated response:
#' # fit <- unitregTMB(Y ~ educ, p0.formula = ~ 1, data = da, family = vasicek())
#' # 
#' # # Extract randomized quantile residuals (default)
#' # res_q <- residuals(fit)
#' # head(res_q)
#' #
#' # # Extract response residuals (Observed - Expected)
#' # res_r <- residuals(fit, type = "response")
#' # head(res_r)
#' }
#' 
#' @export
residuals.unitregTMB <- function(object, type = c("quantile", "cox-snell", "response"), ...) {
  
  type <- match.arg(type)
  Y <- object$Y
  
  if (type == "response") {
    mu_hat <- predict(object, type = "response")
    return(Y - mu_hat)
  }
  
  preds <- predict(object, type = "parameter")
  
  cdf_vals <- get_cdf_continuous(object$obj$env$data$family, 
                                 Y, preds$mu, preds$phi, object$tau)
  
  if (object$has_p0 || object$has_p1) {
    cdf_vals <- correct_for_inflation(Y, cdf_vals, preds$p0, preds$p1)
  }
  
  return(compute_residuals(cdf_vals, type))
}
