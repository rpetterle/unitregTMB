#' @title Simulate Correlated or Independent Bounded Data
#' 
#' @description 
#' Simulates responses from the distributions supported by \code{unitregTMB}. 
#' The function allows the specification of fixed effects for all linear predictors 
#' (\code{mu}, \code{phi}, \code{p0}, \code{p1}), as well as correlated random effects 
#' for the location parameter (\code{mu}) to simulate longitudinal/clustered data. 
#' It can also auto-generate covariates based on user-defined statistical distributions.
#' 
#' @param n_id Integer. Number of clusters/individuals. Default is 100.
#' @param n_rep Integer. Number of repeated measures per individual. Default is 1 (independent data).
#' @param family A \code{unitregTMBFamily} object specifying the distribution and links.
#' @param formula_mu Formula for the location parameter. Can include random effects (e.g., \code{~ x1 + (1 + x1 | id)}).
#' @param formula_phi Formula for the precision/dispersion parameter. Default is \code{~ 1}.
#' @param formula_p0 Formula for the zero-inflation parameter. Default is \code{~ 0} (no inflation).
#' @param formula_p1 Formula for the one-inflation parameter. Default is \code{~ 0} (no inflation).
#' @param beta_mu Numeric vector of fixed effect coefficients for \code{mu}.
#' @param beta_phi Numeric vector of fixed effect coefficients for \code{phi}.
#' @param beta_p0 Numeric vector of fixed effect coefficients for \code{p0}. Default is \code{NULL}.
#' @param beta_p1 Numeric vector of fixed effect coefficients for \code{p1}. Default is \code{NULL}.
#' @param Sigma Covariance matrix for the random effects in \code{mu}. Required if \code{formula_mu} contains random terms.
#' @param cov_config A named list defining the distributions of the covariates present in the formulas. 
#'   Supported distributions: \code{"normal"} (mean, sd), \code{"binomial"} (size, prob), 
#'   \code{"uniform"} (min, max), and \code{"poisson"} (lambda). If a covariate is not specified, 
#'   it defaults to a standard normal distribution.
#' @param seed Optional integer to set the random seed for reproducibility.
#' 
#' @return A \code{data.frame} containing the generated cluster indicator (\code{id}), 
#'   all the simulated covariates, and the simulated response variable (\code{Y}).
#' 
#' @examples
#' \donttest{
#' # 1. Simulate Independent Data (No Random Effects)
#' # vasicek mean model with zero-inflation
#' sim_indep <- sim_unitreg(
#'   n_id = 500, n_rep = 1,
#'   family = vasicek(model_for = "mean"),
#'   formula_mu = ~ cont_var + bin_var,
#'   formula_p0 = ~ 1,
#'   beta_mu = c(0.5, -0.2, 0.8), # Intercept, cont_var, bin_var
#'   beta_phi = c(2.0),           # Intercept (log scale)
#'   beta_p0 = c(-1.5),           # Intercept (logit scale, approx 18% zeros)
#'   cov_config = list(
#'     cont_var = list(dist = "normal", mean = 10, sd = 2),
#'     bin_var  = list(dist = "binomial", size = 1, prob = 0.5)
#'   ),
#'   seed = 123
#' )
#' head(sim_indep)
#' 
#' # 2. Simulate Correlated Longitudinal Data (Random Intercept & Slope)
#' # beta mean model
#' Sigma_matrix <- matrix(c(0.50, 0.15, 
#'                          0.15, 0.20), nrow = 2) # Covariance matrix for (Intercept) and time
#' 
#' sim_long <- sim_unitreg(
#'   n_id = 100, n_rep = 4,
#'   family = beta_fam(model_for = "mean"),
#'   formula_mu = ~ trt + time + (1 + time | id),
#'   formula_phi = ~ trt,
#'   beta_mu = c(-0.5, 1.2, -0.3), # Intercept, trt, time
#'   beta_phi = c(1.5, 0.5),       # Intercept, trt
#'   Sigma = Sigma_matrix,
#'   cov_config = list(
#'     trt  = list(dist = "binomial", size = 1, prob = 0.5),
#'     time = list(dist = "poisson", lambda = 2)
#'   ),
#'   seed = 42
#' )
#' head(sim_long)
#' }
#' 
#' @importFrom stats model.matrix rnorm runif rbinom rpois make.link as.formula
#' @export
sim_unitreg <- function(n_id = 100, n_rep = 1, 
                        family,
                        formula_mu = ~ 1, 
                        formula_phi = ~ 1, 
                        formula_p0 = ~ 0, 
                        formula_p1 = ~ 0,
                        beta_mu, 
                        beta_phi, 
                        beta_p0 = NULL, 
                        beta_p1 = NULL,
                        Sigma = NULL,
                        cov_config = list(),
                        seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  if (!inherits(family, "unitregTMBFamily")) {
    stop("'family' must be a unitregTMBFamily object (e.g., vasicek(model_for = 'mean')).", call. = FALSE)
  }
  
  n_total <- n_id * n_rep
  data_sim <- data.frame(id = rep(1:n_id, each = n_rep))
  
  fixed_mu_form <- reformulas::nobars(formula_mu)
  if (is.null(fixed_mu_form)) fixed_mu_form <- ~ 1
  
  all_formulas <- list(fixed_mu_form, formula_phi, formula_p0, formula_p1)
  all_vars <- unique(unlist(lapply(all_formulas, all.vars)))
  all_vars <- setdiff(all_vars, "id") 
  
  for (var in all_vars) {
    conf <- cov_config[[var]]
    if (is.null(conf)) {
      data_sim[[var]] <- stats::rnorm(n_total, mean = 0, sd = 1)
    } else {
      dist_type <- match.arg(tolower(conf$dist), c("normal", "binomial", "uniform", "poisson"))
      data_sim[[var]] <- switch(dist_type,
                                "normal"   = stats::rnorm(n_total, 
                                                          mean = ifelse(!is.null(conf$mean), conf$mean, 0), 
                                                          sd = ifelse(!is.null(conf$sd), conf$sd, 1)),
                                "binomial" = stats::rbinom(n_total, 
                                                           size = ifelse(!is.null(conf$size), conf$size, 1), 
                                                           prob = ifelse(!is.null(conf$prob), conf$prob, 0.5)),
                                "uniform"  = stats::runif(n_total, 
                                                          min = ifelse(!is.null(conf$min), conf$min, 0), 
                                                          max = ifelse(!is.null(conf$max), conf$max, 1)),
                                "poisson"  = stats::rpois(n_total, 
                                                          lambda = ifelse(!is.null(conf$lambda), conf$lambda, 1))
      )
    }
  }
  
  X_mu  <- stats::model.matrix(fixed_mu_form, data_sim)
  X_phi <- stats::model.matrix(formula_phi, data_sim)
  
  if (ncol(X_mu) != length(beta_mu)) stop("Length of 'beta_mu' does not match the model matrix for mu.", call. = FALSE)
  if (ncol(X_phi) != length(beta_phi)) stop("Length of 'beta_phi' does not match the model matrix for phi.", call. = FALSE)
  
  eta_mu  <- as.vector(X_mu %*% beta_mu)
  eta_phi <- as.vector(X_phi %*% beta_phi)
  
  eta_p0 <- rep(-Inf, n_total) 
  eta_p1 <- rep(-Inf, n_total)
  
  if (!identical(formula_p0, ~ 0) && !is.null(beta_p0)) {
    X_p0 <- stats::model.matrix(formula_p0, data_sim)
    if (ncol(X_p0) != length(beta_p0)) stop("Length of 'beta_p0' does not match the formula.", call. = FALSE)
    eta_p0 <- as.vector(X_p0 %*% beta_p0)
  }
  
  if (!identical(formula_p1, ~ 0) && !is.null(beta_p1)) {
    X_p1 <- stats::model.matrix(formula_p1, data_sim)
    if (ncol(X_p1) != length(beta_p1)) stop("Length of 'beta_p1' does not match the formula.", call. = FALSE)
    eta_p1 <- as.vector(X_p1 %*% beta_p1)
  }
  
  re_terms <- reformulas::findbars(formula_mu)
  if (!is.null(re_terms)) {
    if (is.null(Sigma)) stop("Random effects detected in formula_mu, but 'Sigma' matrix was not provided.", call. = FALSE)
    
    for (bar in re_terms) {
      lhs_expr <- bar[[2]]
      lhs_formula <- stats::as.formula(paste("~", paste(deparse(lhs_expr), collapse = " ")))
      
      Z <- stats::model.matrix(lhs_formula, data_sim)
      q <- ncol(Z)
      
      if (!is.matrix(Sigma) || nrow(Sigma) != q || ncol(Sigma) != q) {
        stop(sprintf("'Sigma' must be a %d x %d covariance matrix to match the random effects structure.", q, q), call. = FALSE)
      }
      
      if (!isSymmetric(Sigma)) {
        stop("'Sigma' must be a symmetric matrix.", call. = FALSE)
      }
      
      L <- tryCatch({
        t(chol(Sigma))
      }, error = function(e) {
        stop("'Sigma' is not positive-definite. Please provide a valid covariance matrix.", call. = FALSE)
      })
      
      u_standard <- matrix(stats::rnorm(n_id * q), nrow = q, ncol = n_id)
      u_correlated <- t(L %*% u_standard) # Result is n_id x q
      
      rand_eff_contribution <- rowSums(Z * u_correlated[data_sim$id, , drop = FALSE])
      eta_mu <- eta_mu + rand_eff_contribution
    }
  }
  
  linkobj_mu <- stats::make.link(family$link_r_name)
  mu_val <- linkobj_mu$linkinv(eta_mu)
  
  linkobj_phi <- stats::make.link(family$phi_link_r_name)
  phi_val <- linkobj_phi$linkinv(eta_phi)
  
  has_p0_model <- any(is.finite(eta_p0))
  has_p1_model <- any(is.finite(eta_p1))
  
  if (has_p0_model && has_p1_model) {
    denom <- 1 + exp(eta_p0) + exp(eta_p1)
    p0_val <- exp(eta_p0) / denom
    p1_val <- exp(eta_p1) / denom
  } else if (has_p0_model) {
    p0_val <- stats::plogis(eta_p0)
    p1_val <- rep(0, n_total)
  } else if (has_p1_model) {
    p0_val <- rep(0, n_total)
    p1_val <- stats::plogis(eta_p1)
  } else {
    p0_val <- rep(0, n_total)
    p1_val <- rep(0, n_total)
  }
  
  y_cont <- get_random_continuous(family$family_code, n_total, mu_val, phi_val, family$tau)
  y_cont <- pmax(pmin(y_cont, 1 - 1e-6), 1e-6)
  
  rand_draws <- stats::runif(n_total)
  cum_p0 <- p0_val
  cum_p1 <- p0_val + p1_val
  
  Y <- y_cont
  Y[rand_draws < cum_p0] <- 0
  Y[rand_draws >= cum_p0 & rand_draws < cum_p1] <- 1
  
  data_sim$Y <- Y
  
  return(data_sim)
}
