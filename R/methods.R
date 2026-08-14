#' @title Summary for unitregTMB Models
#' 
#' @description 
#' Computes and extracts summary statistics for a fitted \code{unitregTMB} model, 
#' including fixed effects, standard errors, z-values, and p-values.
#' 
#' @param object An object of class \code{unitregTMB}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return An object of class \code{summary.unitregTMB} containing the coefficient table, 
#' log-likelihood, model family, link functions, and random effects variance components 
#' (if applicable).
#' 
#' @examples
#' \donttest{
#' # Assuming a dataset 'bodyfat' and a fitted model 'fit':
#' # fit <- unitregTMB(y ~ age + bmi + gender + ipaq + (1 | id), 
#' #                   phi.formula = ~ 1,
#' #                   family = unitgamma(model_for = "mean"),  
#' #                   data = bodyfat_long)
#' #
#' # summary(fit)
#' }
#' @export
summary.unitregTMB <- function(object,...) {
  
  ## Coefficients for mu
  estimates_mu <- object$model_coef$mu
  stderror_mu <- object$model_std$std.mu
  if (is.null(stderror_mu)) {
    stderror_mu <- rep(NA, length(estimates_mu))
  }
  
  ## Coefficients for phi
  estimates_phi <- object$model_coef$phi
  stderror_phi <- object$model_std$std.phi
  if (is.null(stderror_phi)) {
    stderror_phi <- rep(NA, length(estimates_phi))
  }
  
  ## Coefficients for p0 
  estimates_p0 <- NULL
  stderror_p0 <- NULL
  if (isTRUE(object$has_p0)) {
    estimates_p0 <- object$model_coef$p0
    stderror_p0 <- object$model_std$std.p0
    if (is.null(stderror_p0)) {
      stderror_p0 <- rep(NA, length(estimates_p0))
    }
  }
  
  ## Coefficients for p1
  estimates_p1 <- NULL
  stderror_p1 <- NULL
  if (isTRUE(object$has_p1)) {
    estimates_p1 <- object$model_coef$p1
    stderror_p1 <- object$model_std$std.p1
    if (is.null(stderror_p1)) {
      stderror_p1 <- rep(NA, length(estimates_p1))
    }
  }
  
  estimates <- estimates_mu
  stderror <- stderror_mu
  
  n_mu_coefs <- length(estimates_mu)
  n_phi_coefs <- length(estimates_phi)
  n_p0_coefs <- length(estimates_p0)
  n_p1_coefs <- length(estimates_p1)
  
  if (n_phi_coefs > 0) {
    estimates <- c(estimates, estimates_phi)
    stderror <- c(stderror, stderror_phi)
  }
  if (n_p0_coefs > 0) {
    estimates <- c(estimates, estimates_p0)
    stderror <- c(stderror, stderror_p0)
  }
  if (n_p1_coefs > 0) {
    estimates <- c(estimates, estimates_p1)
    stderror <- c(stderror, stderror_p1)
  }
  
  if (length(estimates) == 0) {
    stop("No fixed effect estimates were found.")
  }
  
  stderror_safe <- stderror
  stderror_safe[!is.finite(stderror_safe) | stderror_safe == 0] <- NA
  
  zvalue <- estimates / stderror_safe
  pvalue <- 2 * pnorm(-abs(zvalue))
  
  coeftable <- cbind(
    "Estimate" = estimates,
    "Std. Error" = stderror, 
    "Z value" = zvalue,
    "Pr(>|z|)" = pvalue
  )
  
  rownames_list <- list()
  if (!is.null(names(estimates_mu))) rownames_list$mu <- names(estimates_mu)
  if (!is.null(names(estimates_phi))) rownames_list$phi <- names(estimates_phi)
  if (!is.null(names(estimates_p0))) rownames_list$p0 <- names(estimates_p0)
  if (!is.null(names(estimates_p1))) rownames_list$p1 <- names(estimates_p1)
  
  rownames(coeftable) <- unlist(rownames_list)
  
  phi_link_name <- if (!is.null(object$family_object$phi_link_r_name)) {
     object$family_object$phi_link_r_name
   } else { 
     "log"
 } 

  out <- list(
    coeftable = coeftable,
    logLik = object$logLik,
    call = object$call,
    link.mu = object$link.mu$name,
    link.phi = phi_link_name,
    link.p0 = if (object$has_p0 && object$has_p1) "bivariate logit" else "logit",
    link.p1 = if (object$has_p0 && object$has_p1) "bivariate logit" else "logit",
    family = object$family,
    nobs = object$nobs,
    has_random_effects_mu = object$has_random_effects_mu,
    random_effects_variance = object$random_effects_variance,
    re_info_mu = object$re_info_mu,
    npar = object$npar,
    n_mu_coefs = n_mu_coefs,
    n_phi_coefs = n_phi_coefs,
    n_p0_coefs = n_p0_coefs,
    n_p1_coefs = n_p1_coefs,
    tau = object$tau
  )
  class(out) <- "summary.unitregTMB"
  out
}

#' @title Print Summary for unitregTMB Models
#' 
#' @description 
#' Custom print method for \code{summary.unitregTMB} objects.
#' 
#' @param x An object of class \code{summary.unitregTMB}.
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return Invisibly returns the original \code{summary.unitregTMB} object.
#' 
#' @export
print.summary.unitregTMB <- function(x, digits = max(3, getOption("digits") -3), ...) {
  
  if (!is.matrix(x$coeftable) || !is.numeric(x$coeftable)) {
    stop("The coefficient table is not numeric. The model optimization may have failed.")
  }
  
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  cat("---------------------------------------------------------------\n")

  ## mu coefficients
  cat("Mu coefficients (", x$link.mu, " link):\n", sep = "")
  start_idx <- 1
  end_idx <- x$n_mu_coefs
  
  if (end_idx >= start_idx) {  
    printCoefmat(x$coeftable[start_idx:end_idx, , drop = FALSE], digits = digits, has.Pvalue = TRUE)
  } else {
    cat("  (No coefficients for mu)\n")
  }
  cat("---------------------------------------------------------------\n")
  
  ## phi coefficients
  cat("Phi coefficients (", x$link.phi, " link):\n", sep = "")
  start_idx <- end_idx + 1
  end_idx <- start_idx + x$n_phi_coefs - 1
  if (end_idx >= start_idx) {
    printCoefmat(x$coeftable[start_idx:end_idx, , drop = FALSE], digits = digits, has.Pvalue = TRUE)
  } else {
    cat("  (No coefficients for phi)\n")
  }
  cat("---------------------------------------------------------------\n")
  
  ## p0 coefficients (if applicable)
  if (x$n_p0_coefs > 0) {  
    cat("Zero-inflation (p0) coefficients (", x$link.p0, " link):\n", sep = "")
    start_idx <- end_idx + 1
    end_idx <- start_idx + x$n_p0_coefs - 1
    printCoefmat(x$coeftable[start_idx:end_idx, , drop = FALSE], digits = digits, has.Pvalue = TRUE)
    cat("---------------------------------------------------------------\n")
  }
  
  ## p1 coefficients (if applicable)
  if (x$n_p1_coefs > 0) {  
    cat("One-inflation (p1) coefficients (", x$link.p1, " link):\n", sep = "")
    start_idx <- end_idx + 1
    end_idx <- start_idx + x$n_p1_coefs - 1
    printCoefmat(x$coeftable[start_idx:end_idx, , drop = FALSE], digits = digits, has.Pvalue = TRUE)
    # cat("\n")
    cat("---------------------------------------------------------------\n")
  }
  
  ## Random effects summary
  if (isTRUE(x$has_random_effects_mu) && !is.null(x$random_effects_variance) && length(x$random_effects_variance) > 0) {
    cat("Random Effects:\n")
    
    for (i in seq_along(x$random_effects_variance)) {
      group_name <- if (is.null(x$re_info_mu$group_var_names[i])) "Unknown" else x$re_info_mu$group_var_names[i]
      n_levels <- if (is.null(x$re_info_mu$n_re_levels_list[i])) "NA" else x$re_info_mu$n_re_levels_list[i]
      
      cat(paste0("  Groups: '", group_name, "' (", n_levels, " levels)\n"))
      
      current_re_table <- x$random_effects_variance[[i]]
      
      if (!is.null(current_re_table) && nrow(current_re_table) > 0) {
        std_dev_rows <- startsWith(rownames(current_re_table), "Std.Dev.")
        if (any(std_dev_rows)) {
          cat("  Standard Deviations:\n")
          std_dev_table <- current_re_table[std_dev_rows, , drop = FALSE]
          rownames(std_dev_table) <- gsub("^Std\\.Dev\\.", "    ", rownames(std_dev_table))
          printCoefmat(std_dev_table, digits = digits, has.Pvalue = FALSE)
        }
        
        corr_rows <- startsWith(rownames(current_re_table), "Corr(")
        if (any(corr_rows)) {
          cat("\n  Correlations:\n")
          corr_table <- current_re_table[corr_rows, , drop = FALSE]
          rownames(corr_table) <- paste0("    ", rownames(corr_table))
          printCoefmat(corr_table, digits = digits, has.Pvalue = FALSE)
        }
      }
      cat("\n")
    }
    cat("---------------------------------------------------------------\n") 
  }
  family_string <- paste0("Family: ", x$family)
  if (!is.null(x$tau)) {
    family_string <- paste0(family_string, " (tau = ", x$tau, ")")
  }
  cat(family_string, "\n", sep = "")
  if (is.numeric(x$logLik) && !is.na(x$logLik)) {
    cat("Log-Likelihood: ", format(round(x$logLik[1], 2), nsmall = 2), "\n", sep = "")
  } else {
    cat("Log-Likelihood: NA\n")
  }
  cat("Number of observations: ", x$nobs, "\n", sep = "")
  invisible(x)
}

#' @title Print a unitregTMB Model
#' 
#' @description 
#' Prints the basic information and estimated coefficients of a fitted \code{unitregTMB} model.
#' 
#' @param x An object of class \code{unitregTMB}.
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return Invisibly returns the original \code{unitregTMB} object.
#' 
#' @examples
#' \donttest{
#' # fit <- unitregTMB(Y ~ educ, data = da, family = vasicek())
#' # print(fit)
#' }
#' @export
print.unitregTMB <- function(x, digits = 4, ...) {
  cat("\nCall:  ", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Mu coefficients (", as.character(x$link.mu$name)[1], " link):\n", sep = "")
  print.default(formatC(x$model_coef$mu, digits = digits, width = 4, format = "f"), print.gap = 2, quote = FALSE)
  
  phi_link_name <- if (!is.null(x$family_object$phi_link_r_name)) {
      x$family_object$phi_link_r_name
    } else {
      "log"
    }
  cat("\nPhi coefficients (", phi_link_name, " link):\n", sep = "")

    print.default(formatC(x$model_coef$phi, digits = digits, width = 4, format = "f"), print.gap = 2, quote = FALSE)
    cat("\n")
 
  ## Print for p0 (zero-inflation)
  if (x$has_p0) {
    link_p0_text <- if (x$has_p0 && x$has_p1) "bivariate logit" else "logit"
    cat("Zero-inflation (p0) coefficients (", link_p0_text, " link):\n", sep = "")
    print.default(formatC(x$model_coef$p0, digits = digits, width = 4, format = "f"), print.gap = 2, quote = FALSE)
    cat("\n")
  }
  
  ## Print for p1 (one-inflation)
  if (x$has_p1) {
    link_p1_text <- if (x$has_p0 && x$has_p1) "bivariate logit" else "logit"
    cat("One-inflation (p1) coefficients (", link_p1_text, " link):\n", sep = "")
    print.default(formatC(x$model_coef$p1, digits = digits, width = 4, format = "f"), print.gap = 2, quote = FALSE)
    cat("\n")
  }
  
  family_string <- paste0("Family: ", x$family)
  if (!is.null(x$tau)) {
    family_string <- paste0(family_string, " (tau = ", x$tau, ")")
  }
  cat(family_string, "\n", sep = "")

  if (is.numeric(x$logLik) && !is.na(x$logLik)) {
    cat("Log-Likelihood: ", format(round(x$logLik[1], 2), nsmall = 2), "\n", sep = "")
  } else {
    cat("Log-Likelihood: NA\n")
  }
  invisible(x)
}

#' @name gof_tab.unitregTMB
#'
#' @title Goodness-of-Fit Table for unitregTMB Models
#' 
#' @description 
#' Computes and compares log-likelihood, AIC, and BIC for one or more fitted \code{unitregTMB} models.
#' 
#' @param object An object of class \code{unitregTMB}.
#' @param ... Additional \code{unitregTMB} model objects for comparison.
#' @param digits Number of significant digits to be used in the output.
#' 
#' @return An object of class \code{gof_tab.unitregTMB} containing a list of computed metrics for each model.
#' 
#' @examples
#' \donttest{
#' # Assuming fit1 and fit2 are two fitted unitregTMB models:
#' # fit1 <- unitregTMB(arms ~ age + bmi + gender + ipaq, 
#'                      family = vasicek(model_for = "quantile", tau = 0.75), 
#'                      data = bodyfat)
#' # fit2 <- unitregTMB(arms ~ age + bmi + gender + ipaq, 
#'                      family = kumaraswamy(model_for = "quantile", tau = 0.75), 
#'                      data = bodyfat)
#' # gof_tab(fit1, fit2)
#' }
#' @rdname gof_tab_unitregTMB
#' @export
gof_tab.unitregTMB <- function(object, ..., digits = 2) {
  models <- list(object, ...)
  
  results <- lapply(models, function(model) {
    logLik_val <- model$logLik
    
    if (is.null(logLik_val) || !is.finite(logLik_val)) {
      warning(paste0("logLik for model '", 
                     if (!is.null(model$family)) model$family else "unknown", 
                     "' is not finite or NULL. Setting to NA."))
      logLik_val <- NA_real_
    }
    
    n <- nrow(model$data)
    k <- tryCatch(length(model$opt$par), error = function(e) NA_integer_)
    if (is.na(k)) {
      warning(paste0("Could not determine the number of parameters (k) for model '",
                     if (!is.null(model$family)) model$family else "unknown", 
                     "'. Setting k to NA."))
      AIC_val <- NA_real_
      BIC_val <- NA_real_
    } else {
      AIC_val <- -2 * logLik_val + 2 * k
      BIC_val <- -2 * logLik_val + log(n) * k
    }
    
    if (!is.finite(AIC_val)) {
      AIC_val <- NA_real_
    }
    if (!is.finite(BIC_val)) {
      BIC_val <- NA_real_
    }
    list(
      model_name = if (!is.null(model$family)) model$family else "unitregTMB model",
      tau        = model$tau,
      logLik     = logLik_val,
      AIC        = AIC_val,
      BIC        = BIC_val
    )
  })
  class(results) <- "gof_tab.unitregTMB"
  attr(results, "digits") <- digits
  return(results)
}

#' @rdname gof_tab_unitregTMB
#' @param x An object of class \code{gof_tab.unitregTMB}.
#' @export
print.gof_tab.unitregTMB <- function(x, digits = NULL, ...) {
  
  if (is.null(digits)) {
    digits <- attr(x, "digits")
    if (is.null(digits)) {
      digits <- 3
    }
  }
  
  cat("Goodness-of-fit measures for unitregTMB models:\n")
  
  model_names <- sapply(x, function(m) {
    dist_name <- m$model_name
    dist_name <- sub(" \\(quantile\\)", "", dist_name)
    
    if (!is.null(m$tau) && is.numeric(m$tau)) {
      paste0(dist_name, " (tau = ", m$tau, ")")
    } else {
      dist_name
    }
  })
  
  stats_mat <- rbind(
    LogLik = sapply(x, function(m) m$logLik),
    AIC    = sapply(x, function(m) m$AIC),
    BIC    = sapply(x, function(m) m$BIC)
  )
  colnames(stats_mat) <- model_names
  
  col_width_first <- 10
  col_width_other <- max(nchar(model_names), 15) + 2

  cat(format("Criterion", width = col_width_first, justify = "left"))
  for (name in model_names) {
    cat(format(name, width = col_width_other, justify = "right"))
  }
  cat("\n")
  
  for (stat in rownames(stats_mat)) {
    cat(format(stat, width = col_width_first, justify = "left"))
    for (val in stats_mat[stat, ]) {
      if (!is.finite(val) || is.na(val)) {
        display_val <- "---"
      } else {
        display_val <- format(round(val, digits), nsmall = digits)
      }
      cat(format(display_val, width = col_width_other, justify = "right"))
    }
    cat("\n")
  }
  invisible(x)
}


#' @title Extract Model Coefficients
#' 
#' @description 
#' Extracts the estimated coefficients from a fitted \code{unitregTMB} model.
#' 
#' @param object An object of class \code{unitregTMB}.
#' @param type A character string indicating the type of coefficients to extract. 
#'   Available options are \code{"all"} (default), \code{"mu"} (location), 
#'   \code{"phi"} (precision/dispersion), \code{"p0"} (zero-inflation), or \code{"p1"} (one-inflation).
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return A numeric vector of estimated coefficients.
#' 
#' @examples
#' \donttest{
#' # fit <- unitregTMB(arms ~ age + bmi + gender + ipaq, 
#'                     family = simplex(), 
#'                     data = bodyfat)
#' # coef(fit, type = "all")
#' # coef(fit, type = "mu")
#' }
#' @export
coef.unitregTMB <- function(object, type = c("all", "mu", "phi", "p0", "p1"), ...) {
  
  type <- match.arg(type)
  
  mu_coefs <- object$model_coef$mu
  phi_coefs <- object$model_coef$phi
  p0_coefs <- object$model_coef$p0
  p1_coefs <- object$model_coef$p1
  
  out <- switch(type,
                "all" = {
                  all_coefs <- c(mu_coefs, phi_coefs)
                  if (isTRUE(object$has_p0)) all_coefs <- c(all_coefs, p0_coefs)
                  if (isTRUE(object$has_p1)) all_coefs <- c(all_coefs, p1_coefs)
                  all_coefs
                },
                "mu" = mu_coefs,
                "phi" = phi_coefs,
                "p0" = p0_coefs,
                "p1" = p1_coefs,
                stop("Invalid 'type' argument for coef.unitregTMB.")
  )
  
  return(out)
}

#' @title Extract Log-Likelihood
#' 
#' @description 
#' Extracts the log-likelihood of a fitted \code{unitregTMB} model.
#' 
#' @param object An object of class \code{unitregTMB}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return An object of class \code{logLik} containing the log-likelihood value, 
#'   along with attributes for degrees of freedom (\code{df}) and number of observations (\code{nobs}).
#' 
#' @examples
#' \donttest{
#' # logLik(fit)
#' }
#' @export
logLik.unitregTMB <- function(object, ...) {
  val <- object$logLik
  n_params <- object$npar
  n_obs <- object$nobs
  attr(val, "df") <- n_params
  attr(val, "nobs") <- n_obs
  class(val) <- "logLik"
  
  return(val)
}

#' @title Extract Akaike Information Criterion
#' 
#' @description 
#' Extracts the Akaike Information Criterion (AIC) of a fitted \code{unitregTMB} model.
#' 
#' @param object An object of class \code{unitregTMB}.
#' @param ... Further arguments passed to or from other methods.
#' @param k Numeric, the penalty per parameter to be used; the default \code{k = 2} 
#'   is the classical AIC.
#' 
#' @return A numeric value representing the AIC.
#' 
#' @examples
#' \donttest{
#' # AIC(fit)
#' }
#' @export
AIC.unitregTMB <- function(object, ..., k = 2) {
  ll <- logLik(object)
  return(stats::AIC(ll, k = k))
}

#' @title Extract Bayesian Information Criterion
#' 
#' @description 
#' Extracts the Bayesian Information Criterion (BIC) of a fitted \code{unitregTMB} model.
#' 
#' @param object An object of class \code{unitregTMB}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return A numeric value representing the BIC.
#' 
#' @examples
#' \donttest{
#' # BIC(fit)
#' }
#' @export
BIC.unitregTMB <- function(object, ...) {
  ll <- logLik(object)
  return(stats::BIC(ll))
}

#' @title Calculate Confidence Intervals for Model Parameters
#' 
#' @description 
#' Computes confidence intervals for one or more parameters in a fitted \code{unitregTMB} model
#' using the Wald method (based on standard errors).
#' 
#' @param object A fitted \code{unitregTMB} model.
#' @param parm A specification of which parameters are to be given confidence intervals, 
#'         either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level The confidence level required (default is 0.95).
#' @param ... Additional arguments (currently ignored).
#' 
#' @return A matrix with columns giving lower and upper confidence limits for each parameter.
#' 
#' @examples
#' \donttest{
#' # confint(fit, level = 0.95)
#' }
#' @export
confint.unitregTMB <- function(object, parm, level = 0.95, ...) {
  if (is.null(object$sd_report)) {
    stop("The model must be fitted with sd.report = TRUE to calculate confidence intervals.")
  }
  
  sdr <- summary(object$sd_report, "fixed")
  est <- sdr[, "Estimate"]
  se  <- sdr[, "Std. Error"]
  
  alpha <- 1 - level
  z <- stats::qnorm(c(alpha / 2, 1 - alpha / 2))
  ci <- est + se %o% z
  
  nms_mu  <- paste0("mu.", names(object$model_coef$mu))
  nms_phi <- paste0("phi.", names(object$model_coef$phi))
  nms_p0  <- if(object$has_p0) paste0("p0.", names(object$model_coef$p0)) else character(0)
  nms_p1  <- if(object$has_p1) paste0("p1.", names(object$model_coef$p1)) else character(0)
  
  friendly_names <- c(nms_mu, nms_phi, nms_p0, nms_p1)
  
  rownames(ci) <- friendly_names
  colnames(ci) <- paste(format(100 * c(alpha / 2, 1 - alpha / 2), trim = TRUE, digits = 3), "%")
  
  if (!missing(parm)) {
    ci <- ci[parm, , drop = FALSE]
  }
  
  return(ci)
}

#' @title Compute Profile Likelihood for unitregTMB parameters
#' 
#' @description 
#' Computes the profile likelihood for a specific parameter of a fitted \code{unitregTMB} model.
#' This is computationally intensive but provides more accurate confidence intervals 
#' for asymmetric parameters (like dispersion) than the Wald method.
#' 
#' @param fitted A fitted \code{unitregTMB} model.
#' @param parm The name of the parameter to profile (e.g., "mu.x1", "phi.(Intercept)").
#' @param level.max Maximum confidence level to profile (default is 0.99).
#' @param trace Logical; if TRUE, prints progress.
#' @param ... Additional arguments passed to \code{TMB::tmbprofile}.
#' 
#' @return An object of class \code{profile.unitregTMB} (and \code{tmbprofile}), 
#' containing the computed profile likelihood values.
#' 
#' @export
profile.unitregTMB <- function(fitted, parm, level.max = 0.99, trace = FALSE, ...) {
  
  nms_mu  <- paste0("mu.", names(fitted$model_coef$mu))
  nms_phi <- paste0("phi.", names(fitted$model_coef$phi))
  nms_p0  <- if(fitted$has_p0) paste0("p0.", names(fitted$model_coef$p0)) else character(0)
  nms_p1  <- if(fitted$has_p1) paste0("p1.", names(fitted$model_coef$p1)) else character(0)
  
  friendly_names <- c(nms_mu, nms_phi, nms_p0, nms_p1)
  
  if (missing(parm) || length(parm) != 1 || !(parm %in% friendly_names)) {
    stop("Please specify exactly one valid parameter name to profile (e.g., 'mu.x1' or 'phi.(Intercept)').\n",
         "Available parameters: ", paste(friendly_names, collapse = ", "))
  }
  
  idx <- which(friendly_names == parm)
  
  lc <- rep(0, length(fitted$opt$par))
  lc[idx] <- 1
  
  if(trace) cat(sprintf("Profiling parameter '%s'...\n", parm))
  
  prof <- TMB::tmbprofile(fitted$obj, 
                          name = parm, 
                          lincomb = lc, 
                          trace = trace, 
                          ...)
  
  class(prof) <- c("profile.unitregTMB", class(prof))
  return(prof)
}

#' @title Extract the Number of Observations
#' 
#' @description 
#' Extracts the number of observations used to fit a \code{unitregTMB} model.
#' 
#' @param object An object of class \code{unitregTMB}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return An integer representing the number of observations.
#' 
#' @export
nobs.unitregTMB <- function(object, ...) {
  return(object$nobs)
}

#' @title Extract the Residual Degrees of Freedom
#' 
#' @description 
#' Extracts the residual degrees of freedom from a fitted \code{unitregTMB} model.
#' 
#' @param object An object of class \code{unitregTMB}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return An integer representing the residual degrees of freedom.
#' 
#' @export
df.residual.unitregTMB <- function(object, ...) {
  return(object$nobs - object$npar)
}

#' Extract Fixed Effects
#'
#' Generic function to extract fixed effects from fitted models.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments.
#' @export
fixed_effects <- function(object, ...) 
UseMethod("fixed_effects")

#' Extract Fixed Effects from a unitregTMB Object
#'
#' @param object A fitted object of class "unitregTMB".
#' @param ... Additional arguments.
#' @return A numeric vector of fixed effects.
#' @rdname fixed_effects
#' @export
fixed_effects.unitregTMB <- function(object, ...) {
  res <- object$model_coef
  class(res) <- "fixef.unitregTMB"
  return(res)
}

#' @title Print Method for fixef.unitregTMB
#' 
#' @description 
#' Custom print method for \code{fixef.unitregTMB} objects.
#' 
#' @param x An object of class \code{fixef.unitregTMB}.
#' @param digits Number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return Invisibly returns the original \code{fixef.unitregTMB} object.
#' 
#' @export
print.fixef.unitregTMB <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\n--- Location Model (mu) ---\n")
  print.default(format(x$mu, digits = digits), print.gap = 2L, quote = FALSE)
  
  cat("\n--- Precision/Dispersion Model (phi) ---\n")
  print.default(format(x$phi, digits = digits), print.gap = 2L, quote = FALSE)
  
  if (!is.null(x$p0)) {
    cat("\n--- Zero-Inflation Model (p0) ---\n")
    print.default(format(x$p0, digits = digits), print.gap = 2L, quote = FALSE)
  }
  if (!is.null(x$p1)) {
    cat("\n--- One-Inflation Model (p1) ---\n")
    print.default(format(x$p1, digits = digits), print.gap = 2L, quote = FALSE)
  }
  invisible(x)
}

#' @title Extract the Design Matrix
#' 
#' @description 
#' Extracts the design matrix for a specific model component of a \code{unitregTMB} object.
#' 
#' @param object An object of class \code{unitregTMB}.
#' @param component Character string indicating which design matrix to extract: 
#'   \code{"mu"} (location), \code{"phi"} (precision/dispersion), \code{"p0"} (zero-inflation), 
#'   or \code{"p1"} (one-inflation).
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return A design matrix representing the chosen model component.
#' 
#' @examples
#' \donttest{
#' # model.matrix(fit, component = "mu")
#' }
#' @export
model.matrix.unitregTMB <- function(object, component = c("mu", "phi", "p0", "p1"), ...) {
  component <- match.arg(component)
  form <- reformulas::nobars(formula(object, component = component))
  if (is.null(form)) return(NULL)
  mf <- stats::model.frame(form, object$data, na.action = stats::na.pass)
  return(stats::model.matrix(form, mf))
}

#' @title Extract the Model Formula
#' 
#' @description 
#' Extracts the formula for a specific model component of a \code{unitregTMB} object.
#' 
#' @param x An object of class \code{unitregTMB}.
#' @param component Character string indicating which formula to extract: 
#'   \code{"mu"} (location), \code{"phi"} (precision/dispersion), \code{"p0"} (zero-inflation), 
#'   or \code{"p1"} (one-inflation).
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return An object of class \code{formula}.
#' 
#' @examples
#' \donttest{
#' # formula(fit, component = "phi")
#' }
#' @export
formula.unitregTMB <- function(x, component = c("mu", "phi", "p0", "p1"), ...) {
  component <- match.arg(component)
  switch(component,
         "mu"  = x$formula,
         "phi" = x$phi.formula,
         "p0"  = x$p0.formula,
         "p1"  = x$p1.formula)
}

#' @title Extract Model Fitted Values for unitregTMB
#' 
#' @description 
#' Computes the fitted values for a \code{unitregTMB} model. This method reconstructs the linear
#' predictors using the original design matrices stored in the TMB object environment
#' and the optimized parameters (including random effects and offsets).
#' 
#' @param object An object of class \code{unitregTMB}.
#' @param type The type of prediction required. 
#' \itemize{
#'   \item \code{"response"}: The expected mean value of the response variable (default).
#'   \item \code{"mu"}: The location parameter of the continuous component.
#'   \item \code{"phi"}: The precision parameter.
#'   \item \code{"p0"}: The probability of zero inflation.
#'   \item \code{"p1"}: The probability of one inflation.
#' }
#' @param ... Additional arguments (currently unused).
#' 
#' @return A numeric vector of fitted values corresponding to the specified \code{type}.
#' 
#' @examples
#' \donttest{
#' # fitted(fit, type = "response")
#' # fitted(fit, type = "phi")
#' }
#' @export
fitted.unitregTMB <- function(object, type = c("response", "mu", "phi", "p0", "p1"), ...) {
  
  type <- match.arg(type)
  
  env_data <- object$obj$env$data
  
  X_mu  <- env_data$X_mu
  X_phi <- env_data$X_phi
  X_p0  <- env_data$X_p0
  X_p1  <- env_data$X_p1
  Z_mu  <- env_data$Z_mu
  
  beta_mu  <- object$model_coef$mu
  beta_phi <- object$model_coef$phi
  beta_p0  <- object$model_coef$p0
  beta_p1  <- object$model_coef$p1
  
  eta_mu <- as.vector(X_mu %*% beta_mu)
  
  if (!is.null(env_data$offset_mu)) {
    eta_mu <- eta_mu + env_data$offset_mu
  }
  
  if (object$has_random_effects_mu) {
    full_par <- object$obj$env$last.par.best
    u_mu_hat <- full_par[names(full_par) == "u_mu"]
    
    if (length(u_mu_hat) == ncol(Z_mu)) {
      eta_mu <- eta_mu + as.vector(Z_mu %*% u_mu_hat)
    } else {
      warning("Dimension mismatch: 'u_mu' length does not match Z matrix columns. Random effects ignored in fitted values.")
    }
  }
  
  mu_vals <- object$link.mu$linkinv(eta_mu)
  if (type == "mu") return(mu_vals)
  
  eta_phi <- as.vector(X_phi %*% beta_phi)
  linkobj_phi <- stats::make.link(object$family_object$phi_link_r_name)
  phi_vals <- linkobj_phi$linkinv(eta_phi)
  if (type == "phi") return(phi_vals)
  
  n_obs <- nrow(X_mu)
  eta_p0 <- if (object$has_p0) as.vector(X_p0 %*% beta_p0) else rep(0, n_obs)
  eta_p1 <- if (object$has_p1) as.vector(X_p1 %*% beta_p1) else rep(0, n_obs)
  
  if (object$has_p0 && object$has_p1) {
    denom <- 1 + exp(eta_p0) + exp(eta_p1)
    p0_vals <- exp(eta_p0) / denom
    p1_vals <- exp(eta_p1) / denom
  } else if (object$has_p0) {
    p0_vals <- stats::plogis(eta_p0)
    p1_vals <- rep(0, n_obs)
  } else if (object$has_p1) {
    p0_vals <- rep(0, n_obs)
    p1_vals <- stats::plogis(eta_p1)
  } else {
    p0_vals <- rep(0, n_obs)
    p1_vals <- rep(0, n_obs)
  }
  
  if (type == "p0") return(p0_vals)
  if (type == "p1") return(p1_vals)
  
  fitted_resp <- (1 - p0_vals - p1_vals) * mu_vals + p1_vals
  
  return(fitted_resp)
}

#' @title Extract Variance-Covariance Matrix
#' 
#' @description 
#' Extracts the variance-covariance matrix of the main parameters of a fitted \code{unitregTMB} model.
#' 
#' @param object An object of class \code{unitregTMB}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return A numeric matrix representing the estimated covariances between the parameter estimates.
#' 
#' @examples
#' \donttest{
#' # vcov(fit)
#' }
#' @export
vcov.unitregTMB <- function(object, ...) {
  
  if (!is.null(object$sd_report) && !is.null(object$sd_report$cov.fixed)) {
    V <- object$sd_report$cov.fixed
  } 
  else if (!is.null(object$opt$hessian)) {
    V <- tryCatch(solve(object$opt$hessian), error = function(e) NULL)
    if (is.null(V)) {
      stop("Hessian is computationally singular. Variance-covariance matrix cannot be computed.")
    }
  } else {
    stop("Variance-covariance matrix not available. Please fit the model with control = unitregTMB.control(sd.report = TRUE).")
  }
  
  tmb_names <- rownames(V)
  
  idx_mu  <- which(tmb_names == "beta_mu")
  idx_phi <- which(tmb_names == "beta_phi")
  idx_p0  <- if (object$has_p0) which(tmb_names == "beta_p0") else integer(0)
  idx_p1  <- if (object$has_p1) which(tmb_names == "beta_p1") else integer(0)
  
  idx_fixed <- c(idx_mu, idx_phi, idx_p0, idx_p1)
  
  V_fixed <- V[idx_fixed, idx_fixed, drop = FALSE]
  
  coef_names <- names(coef(object, type = "all"))
  
  if (length(coef_names) == nrow(V_fixed)) {
    rownames(V_fixed) <- colnames(V_fixed) <- coef_names
  } else {
    warning("Mismatch between number of coefficients and covariance matrix dimensions.")
  }
  
  return(V_fixed)
}

#' @title Update and Re-fit a unitregTMB Model
#' 
#' @description 
#' \code{update} will update and (by default) re-fit a \code{unitregTMB} model. It 
#' allows you to easily change the formula, data, or other arguments without having 
#' to specify the entire model call again.
#' 
#' @param object A fitted \code{unitregTMB} model object.
#' @param formula. Changes to the formula. This is a two-sided formula where the 
#'   left-hand side is usually empty (e.g., \code{~ . + x3}). See \code{\link[stats]{update.formula}} for details.
#' @param ... Additional arguments to the call, or arguments with changed values 
#'   (e.g., \code{phi.formula}, \code{p0.formula}, \code{data}, \code{family}, \code{control}).
#' @param evaluate Logical. If \code{TRUE} (default), the updated model is evaluated and re-fitted. 
#'   If \code{FALSE}, the updated call is returned without evaluating.
#' 
#' @return A new \code{unitregTMB} model object (if \code{evaluate = TRUE}), 
#'   or the unevaluated model call (if \code{evaluate = FALSE}).
#' 
#' @examples
#' \donttest{
#' # Assuming 'da' is your dataset:
#' # fit1 <- unitregTMB(y ~ bmi + gender + ipaq + (1 | id), 
#' #                   family = unitgamma(model_for = "mean"),  
#' #                   data = bodyfat_long))
#' #
#' # # 1. Update the model by adding a covariate to the mean formula
#' # fit2 <- update(fit1, formula. = ~ . + age)
#' #
#' # # 2. Update the model by adding a dispersion formula and changing the family
#' # fit3 <- update(fit1, phi.formula = ~ bmi, family = kumaraswamy())
#' #
#' # # 3. Just get the updated call without fitting (useful for debugging)
#' # call_only <- update(fit1, phi.formula = ~ age, evaluate = FALSE)
#' # print(call_only)
#' }
#' 
#' @importFrom stats update.formula
#' @export
update.unitregTMB <- function(object, formula., ..., evaluate = TRUE) {
  
  if (is.null(object$call)) {
    stop("The model object does not contain a 'call' component.")
  }
  
  call <- object$call
  
  extras <- match.call(expand.dots = FALSE)$...
  
  if (!missing(formula.)) {
    call$formula <- stats::update.formula(object$formula, formula.)
  }
  
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    
    for (a in names(extras)[existing]) {
      call[[a]] <- extras[[a]]
    }
    
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  
  if (evaluate) {
    eval(call, parent.frame())
  } else {
    call
  }
}
