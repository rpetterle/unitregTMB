#' @importFrom stats AIC BIC anova coef confint df.residual fitted formula
#' @importFrom stats logLik model.matrix nobs predict profile residuals simulate
#' @importFrom stats update vcov printCoefmat setNames nlminb optim symnum
#'
#' @importFrom regcore df_envelopes extract_equations gof_tab
#' @importFrom regcore pairwise_vuong_test plot_envelope plot_ranef
#' @importFrom regcore random_effects stepCriterion summary_coef
#' @importFrom regcore vuong_test write_methods
NULL

.unitreg_components <- c("mu", "phi", "p0", "p1")

.unitreg_labels <- c(mu  = "Mu", phi = "Phi", 
                     p0  = "Zero-inflation (p0)",
                     p1  = "One-inflation (p1)")

.unitreg_active_component <- function(object, component) {
  switch(component, mu = TRUE, phi = TRUE,
         p0 = isTRUE(object$has_p0),
         p1 = isTRUE(object$has_p1))
}

.unitreg_coef_table <- function(object, component) {
  if (!.unitreg_active_component(object, component)) return(NULL)

  est <- object$model_coef[[component]]
  se <- object$model_std[[paste0("std.", component)]]
  if (is.null(est)) return(NULL)
  if (is.null(se)) se <- rep(NA_real_, length(est))

  se_safe <- se
  se_safe[!is.finite(se_safe) | se_safe == 0] <- NA_real_
  z <- est / se_safe

  out <- cbind(Estimate = est, 
               `Std. Error` = se,
               `Z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  rownames(out) <- names(est)
  out
}

.unitreg_coef_names <- function(object) {
  unlist(lapply(.unitreg_components, function(component) {
    if (!.unitreg_active_component(object, component)) return(character())
    nm <- names(object$model_coef[[component]])
    paste0(component, ".", nm)
  }), use.names = FALSE)
}

.unitreg_coef_values <- function(object) {
  unlist(lapply(.unitreg_components, function(component) {
    if (!.unitreg_active_component(object, component)) return(numeric())
    object$model_coef[[component]]
  }), use.names = FALSE)
}

.unitreg_coef_se <- function(object) {
  unlist(lapply(.unitreg_components, function(component) {
    if (!.unitreg_active_component(object, component)) return(numeric())
    se <- object$model_std[[paste0("std.", component)]]
    if (is.null(se)) rep(NA_real_, length(object$model_coef[[component]])) else se
  }), use.names = FALSE)
}

.unitreg_link_names <- function(object) {
  c(mu = object$link.mu$name,
    phi = if (!is.null(object$family_object$phi_link_r_name))
             object$family_object$phi_link_r_name else "log",

    p0 = if (isTRUE(object$has_p0) && isTRUE(object$has_p1))
             "bivariate logit" else "logit",

    p1 = if (isTRUE(object$has_p0) && isTRUE(object$has_p1))
             "bivariate logit" else "logit")
}

.unitreg_model_for <- function(object) {
  out <- object$family_object$model_for
  if (!is.null(out) && length(out)) return(as.character(out[1L]))
  code <- object$family_object$family_code
  if (is.null(code) || !length(code) || !is.finite(code[1L])) return(NA_character_)
  code <- as.integer(code[1L])
  if (code %in% 0:4) return("mean")
  if (code %in% 5:8) return("mode")
  if (code %in% 9:15) return("quantile")
  NA_character_
}

.unitreg_family_label <- function(object, include_tau = TRUE) {
  label <- if (!is.null(object$family) && length(object$family)) as.character(object$family[1L]) else "unitregTMB"
  model_for <- .unitreg_model_for(object)
  if (!is.na(model_for) && !grepl(paste0("(", model_for, ")"), label, fixed = TRUE)) {
    label <- paste0(label, " [", model_for, "]")
  }
  if (isTRUE(include_tau) && !is.null(object$tau) && length(object$tau) == 1L && is.finite(object$tau)) {
    label <- paste0(label, " (tau = ", format(object$tau, trim = TRUE), ")")
  }
  label
}

.unitreg_re_list <- function(object) {
  out <- object$random_effects_variance_by_component
  if (!is.null(out)) return(out)

  list(mu = object$random_effects_variance,
       phi = NULL,
       p0 = NULL,
       p1 = NULL)
}

.unitreg_re_info <- function(object, component) {
  object[[paste0("re_info_", component)]]
}

.unitreg_has_re <- function(object, component) {
  isTRUE(object[[paste0("has_random_effects_", component)]])
}

.unitreg_has_spatial <- function(object, component) {
  isTRUE(object[[paste0("has_spatial_effects_", component)]])
}

.unitreg_spatial_info <- function(object, component) {
  object[[paste0("spatial_info_", component)]]
}

.unitreg_spatial_parameters <- function(object, component) {
  object[[paste0("spatial_parameters_", component)]]
}

.unitreg_print_re <- function(label, variance_list, re_info, digits) {
  if (is.null(variance_list) || !length(variance_list)) return(invisible(NULL))

  cat("Random Effects in ", label, ":\n", sep = "")

  for (i in seq_along(variance_list)) {
    group <- if (!is.null(re_info$group_var_names[i]))
      re_info$group_var_names[i] else "Unknown"
    nlev <- if (!is.null(re_info$n_re_levels_list[i]))
      re_info$n_re_levels_list[i] else NA_integer_

    cat("  Groups: '", group, "' (", nlev, " levels)\n", sep = "")

    tab <- variance_list[[i]]
    if (is.null(tab) || !nrow(tab)) {
      cat("\n")
      next
    }

    sd_rows <- startsWith(rownames(tab), "Std.Dev.")
    if (any(sd_rows)) {
      sd_tab <- tab[sd_rows, , drop = FALSE]
      rownames(sd_tab) <- sub("^Std\\.Dev\\.", "    ", rownames(sd_tab))
      cat("  Standard Deviations:\n")
      printCoefmat(sd_tab, digits = digits, has.Pvalue = FALSE)
    }

    cor_rows <- startsWith(rownames(tab), "Corr(")
    if (any(cor_rows)) {
      cor_tab <- tab[cor_rows, , drop = FALSE]
      rownames(cor_tab) <- paste0("    ", rownames(cor_tab))
      cat("\n  Correlations:\n")
      printCoefmat(cor_tab, digits = digits, has.Pvalue = FALSE)
    }
    cat("\n")
  }

  invisible(NULL)
}

.unitreg_print_spatial <- function(label, info, parameters, digits) {
  if (is.null(info) || !isTRUE(info$has_spatial_effects)) return(invisible(NULL))

  cat("Spatial Effect in ", label, ":\n", sep = "")
  cat("  Structure: ", info$structure, "\n", sep = "")
  cat("  Region variable: ", info$region_var, "\n", sep = "")
  cat("  Number of regions: ", info$n_regions, "\n", sep = "")

  if (!is.null(parameters) && nrow(parameters)) {
    tab <- parameters
    rownames(tab) <- paste0("    ", rownames(tab))
    cat("  Spatial Parameters:\n")
    printCoefmat(tab, digits = digits, has.Pvalue = FALSE)
  }

  cat("\n")
  invisible(NULL)
}

.unitreg_full_par <- function(object) {
  flags <- c(vapply(.unitreg_components, function(x) .unitreg_has_re(object, x), logical(1)),
             vapply(.unitreg_components, function(x) .unitreg_has_spatial(object, x), logical(1)))

  if (!any(flags)) return(NULL)

  object$obj$fn(object$opt$par)
  object$obj$env$parList( x = object$opt$par, par = object$obj$env$last.par)
}

.unitreg_compact_model_label <- function(model) {
  label <- .unitreg_family_label(model, include_tau = FALSE)

  role <- NULL

  if (!is.null(model$family_object$model_for) &&
      length(model$family_object$model_for)) {
    role <- as.character(model$family_object$model_for)[1L]
  }

  if (is.null(role) || !role %in% c("mean", "mode", "quantile")) {
    code <- model$family_object$family_code

    if (!is.null(code) && length(code)) {
      code <- as.integer(code)[1L]

      if (!is.na(code) && code %in% 0:4) role <- "mean"
      if (!is.na(code) && code %in% 5:8) role <- "mode"
      if (!is.na(code) && code %in% 9:15) role <- "quantile"
    }
  }

  label <- gsub( "\\s*(\\(quantile\\)|\\[quantile\\])", "", label, ignore.case = TRUE)

  label <- trimws(label)

  if (!is.null(role) && role %in% c("mean", "mode") &&
      !grepl(paste0("(", role, ")"), label, fixed = TRUE)) {
    label <- paste0(label, " (", role, ")")
  }

  if (identical(role, "quantile") && !is.null(model$tau) &&
      length(model$tau) == 1L && is.finite(model$tau)) {
    label <- paste0(label, " (tau = ", format(model$tau, trim = TRUE), ")")
  }
  label
}


.unitreg_print_blocks <- function(tab, row_labels, first_header, 
                                  min_col_width = 12L, 
                                  max_width = getOption("width", 80L)) {

  if (!is.matrix(tab)) tab <- as.matrix(tab)

  if (!nrow(tab) || !ncol(tab)) return(invisible(NULL))

  width1 <- max(nchar(c(first_header, row_labels), type = "width"),
                na.rm = TRUE) + 2L

  col_widths <- vapply(seq_len(ncol(tab)), function(j) {
    as.integer(max(min_col_width, nchar(colnames(tab)[j], type = "width"),
               nchar(tab[, j], type = "width"), na.rm = TRUE) + 2L)
  }, integer(1))

  max_width <- max(40L, as.integer(max_width))

  blocks <- list()
  current <- integer()
  used <- width1

  for (j in seq_len(ncol(tab))) {
    w <- col_widths[j]

    if (length(current) && used + w > max_width) {
      blocks[[length(blocks) + 1L]] <- current
      current <- j
      used <- width1 + w
    } else {
      current <- c(current, j)
      used <- used + w
    }
  }

  if (length(current)) {
    blocks[[length(blocks) + 1L]] <- current
  }

  for (b in seq_along(blocks)) {
    idx <- blocks[[b]]

    cat(format(first_header, width = width1, justify = "left"))

    for (j in idx) {
      cat(format(colnames(tab)[j], width = col_widths[j], justify = "right"))
    }
    cat("\n")

    for (i in seq_len(nrow(tab))) {
      cat(format(row_labels[i], width = width1, justify = "left"))

      for (j in idx) {
        cat(format(tab[i, j], width = col_widths[j], justify = "right"))
      }
      cat("\n")
    }

    if (b < length(blocks)) cat("\n")
  }

  invisible(NULL)
}

#' Summarize a unitregTMB Model
#'
#' Computes Wald summaries for the fixed effects and collects ordinary random
#' effects and proper CAR spatial components for all active distributional
#' predictors.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param ... Additional arguments, currently unused.
#'
#' @return An object of class \code{summary.unitregTMB}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' summary(fit)
#' }
#' @export
summary.unitregTMB <- function(object, ...) {
  coef_tables <- setNames(
    lapply(.unitreg_components, function(x) .unitreg_coef_table(object, x)),
    .unitreg_components
  )
  ncoef <- vapply(coef_tables, function(x) if (is.null(x)) 0L else nrow(x), integer(1))
  coeftable <- do.call(rbind, coef_tables[!vapply(coef_tables, is.null, logical(1))])

  re_by_component <- .unitreg_re_list(object)
  spatial_info <- setNames(
    lapply(.unitreg_components, function(x) .unitreg_spatial_info(object, x)),
    .unitreg_components
  )
  spatial_parameters <- setNames(
    lapply(.unitreg_components, function(x) .unitreg_spatial_parameters(object, x)),
    .unitreg_components
  )

  out <- list(
    call = object$call,
    coeftable = coeftable,
    coefficients = coef_tables,
    links = .unitreg_link_names(object),
    family = object$family,
    model_for = .unitreg_model_for(object),
    tau = object$tau,
    logLik = object$logLik,
    nobs = object$nobs,
    npar = object$npar,
    ncoef = ncoef,
    random_effects_variance_by_component = re_by_component,
    random_effects_variance = object$random_effects_variance,
    spatial_info_by_component = spatial_info,
    spatial_parameters_by_component = spatial_parameters
  )

  for (component in .unitreg_components) {
    out[[paste0("has_random_effects_", component)]] <- .unitreg_has_re(object, component)
    out[[paste0("re_info_", component)]] <- .unitreg_re_info(object, component)
    out[[paste0("has_spatial_effects_", component)]] <- .unitreg_has_spatial(object, component)
    out[[paste0("spatial_info_", component)]] <- spatial_info[[component]]
    out[[paste0("spatial_parameters_", component)]] <- spatial_parameters[[component]]
  }

  out$n_mu_coefs <- ncoef["mu"]
  out$n_phi_coefs <- ncoef["phi"]
  out$n_p0_coefs <- ncoef["p0"]
  out$n_p1_coefs <- ncoef["p1"]

  out$link.mu <- out$links["mu"]
  out$link.phi <- out$links["phi"]
  out$link.p0 <- out$links["p0"]
  out$link.p1 <- out$links["p1"]

  class(out) <- "summary.unitregTMB"
  out
}
#' Print a unitregTMB Summary
#'
#' Prints fixed-effect summaries followed by ordinary random effects, proper CAR
#' spatial effects, and model-level information.
#'
#' @param x An object of class \code{summary.unitregTMB}.
#' @param digits Number of significant digits used in printed tables.
#' @param ... Additional arguments, currently unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' print(summary(fit))
#' }
#' @noRd
#' @export
print.summary.unitregTMB <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  if (!is.matrix(x$coeftable) || !is.numeric(x$coeftable))
    stop("The coefficient table is not numeric; model fitting may have failed.")

  cat("\nCall:  ", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")

  for (component in .unitreg_components) {
    tab <- x$coefficients[[component]]
    if (is.null(tab)) next

    cat(.unitreg_labels[component], " coefficients (",
        x$links[component], " link):\n", sep = "")
    printCoefmat(tab, digits = digits, has.Pvalue = TRUE)
    cat("\n")
  }

  re_by_component <- x$random_effects_variance_by_component
  for (component in .unitreg_components) {
    if (!isTRUE(x[[paste0("has_random_effects_", component)]])) next
    .unitreg_print_re(
      .unitreg_labels[component],
      re_by_component[[component]],
      x[[paste0("re_info_", component)]],
      digits
    )
  }

  for (component in .unitreg_components) {
    if (!isTRUE(x[[paste0("has_spatial_effects_", component)]])) next
    .unitreg_print_spatial(
      .unitreg_labels[component],
      x$spatial_info_by_component[[component]],
      x$spatial_parameters_by_component[[component]],
      digits
    )
  }

  family <- paste0("Family: ", x$family)
  if (!is.null(x$model_for) && length(x$model_for) && !is.na(x$model_for) &&
      !grepl(paste0("(", x$model_for, ")"), x$family, fixed = TRUE)) {
    family <- paste0(family, " [", x$model_for, "]")
  }
  if (!is.null(x$tau)) family <- paste0(family, " (tau = ", x$tau, ")")

  cat(family, "\n", sep = "")
  cat("Log-Likelihood: ",
      if (is.finite(x$logLik)) format(round(x$logLik, 2), nsmall = 2) else "NA",
      "\n", sep = "")
  cat("Number of observations: ", x$nobs, "\n", sep = "")

  invisible(x)
}
#' Print a unitregTMB Model
#'
#' Prints the model call, fixed-effect estimates, active spatial structures, and
#' basic model information. Use \code{summary()} for inferential tables and
#' variance-component details.
#'
#' @param x A fitted \code{unitregTMB} object.
#' @param digits Number of digits used to print estimates.
#' @param ... Additional arguments, currently unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' print(fit)
#' }
#' @export
print.unitregTMB <- function(x, digits = 4, ...) {
  cat("\nCall:  ", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")

  links <- .unitreg_link_names(x)
  for (component in .unitreg_components) {
    if (!.unitreg_active_component(x, component)) next
    b <- x$model_coef[[component]]
    if (is.null(b)) next

    cat(.unitreg_labels[component], " coefficients (",
        links[component], " link):\n", sep = "")
    print.default(
      formatC(b, digits = digits, format = "f"),
      print.gap = 2, quote = FALSE
    )
    cat("\n")
  }

  for (component in .unitreg_components) {
    if (!.unitreg_has_spatial(x, component)) next
    info <- .unitreg_spatial_info(x, component)
    par <- .unitreg_spatial_parameters(x, component)

    cat("Spatial effect in ", component, ": ", info$structure,
        " (", info$n_regions, " regions)", sep = "")
    if (!is.null(par)) {
      cat(", Scale = ", formatC(par["Scale", "Estimate"], digits = digits, format = "f"),
          ", Rho = ", formatC(par["Rho", "Estimate"], digits = digits, format = "f"))
    }
    cat("\n")
  }

  if (any(vapply(.unitreg_components, function(z) .unitreg_has_spatial(x, z), logical(1))))
    cat("\n")

  cat("Family: ", .unitreg_family_label(x), "\n", sep = "")
  cat("Log-Likelihood: ",
      if (is.finite(x$logLik)) format(round(x$logLik, 2), nsmall = 2) else "NA",
      "\n", sep = "")
  cat("Number of observations: ", x$nobs, "\n", sep = "")

  invisible(x)
}
#' Compare Goodness-of-Fit Measures from unitregTMB Models
#'
#' @description
#' Compares log-likelihood, Akaike information criterion (AIC), and Bayesian
#' information criterion (BIC) across one or more fitted
#' \code{unitregTMB} models.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param ... Additional fitted \code{unitregTMB} objects.
#' @param digits Non-negative integer giving the number of decimal places
#'   used when printing the goodness-of-fit measures.
#'
#' @details
#' Model labels are automatically constructed from the fitted family and
#' parameterization. Mean and modal models are displayed, for example, as
#' \code{"Beta (mean)"} and \code{"Beta (mode)"}.
#'
#' For quantile regression, the quantile level is displayed directly in the
#' model label. For example, a Kumaraswamy median regression is shown as
#' \code{"Kumaraswamy (tau = 0.5)"}.
#'
#' When many models are compared, the print method automatically divides the
#' table into blocks according to the current console width, as determined by
#' \code{getOption("width")}.
#'
#' @return An object of class \code{gof_tab.unitregTMB}. It is a list with
#'   one element per fitted model, containing the model label,
#'   log-likelihood, AIC, and BIC.
#'
#' @examples
#' \donttest{
#' data("wqi", package = "unitregTMB")
#'
#' fit_beta <- unitregTMB(y ~ location + quarter + (1 | id),
#'                        family = beta_fam(model_for = "mean"),
#'                        data = wqi)
#'
#' fit_ug <- unitregTMB(y ~ location + quarter + (1 | id),
#'                      family = unitgamma(model_for = "mean"),
#'                      data = wqi)
#'
#' gof_tab(fit_beta, fit_ug)
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom regcore gof_tab
#' @method gof_tab unitregTMB
#' @rdname gof_tab_unitregTMB
#' @export
gof_tab.unitregTMB <- function(object, ..., digits = 2) {

  models <- list(object, ...)

  if (!length(models)) {
    stop("At least one fitted 'unitregTMB' model must be supplied.", call. = FALSE)
  }

  ok <- vapply(
    models,
    inherits,
    logical(1),
    what = "unitregTMB"
  )

  if (!all(ok)) {
    stop(
      "All supplied objects must inherit from 'unitregTMB'.",
      call. = FALSE
    )
  }

  if (!is.numeric(digits) || length(digits) != 1L ||  !is.finite(digits) ||
      digits < 0 || digits != as.integer(digits)) {
    stop( "'digits' must be a non-negative integer.", call. = FALSE)
  }

  digits <- as.integer(digits)

  model_labels <- make.unique(
    vapply(models, .unitreg_compact_model_label,  character(1)), sep = " #")

  out <- lapply(seq_along(models), function(i) {
    model <- models[[i]]

    ll <- as.numeric(logLik(model))

    list(model_name = model_labels[i],
         logLik = ll, 
         AIC = AIC(model), 
         BIC = BIC(model))
  })

  class(out) <- c("gof_tab.unitregTMB", "list")

  attr(out, "digits") <- digits
  out
}


#' Print a gof_tab.unitregTMB Object
#'
#' @description
#' Prints goodness-of-fit measures returned by
#' \code{\link{gof_tab.unitregTMB}}.
#'
#' @param x An object of class \code{gof_tab.unitregTMB}.
#' @param digits Optional non-negative integer overriding the number of
#'   decimal places stored in \code{x}.
#' @param ... Additional arguments, currently unused.
#'
#' @details
#' Columns are sized independently. If the complete table is wider than the
#' available console width, the output is automatically divided into blocks.
#'
#' The available console width can be modified with, for example,
#' \code{options(width = 120)}.
#'
#' @return Invisibly returns \code{x}.
#'
#' @author Ricardo Rasmussen Petterle
#' @method print gof_tab.unitregTMB
#' @rdname print_gof_tab_unitregTMB
#' @export
print.gof_tab.unitregTMB <- function(x, digits = NULL, ...) {

  if (!inherits(x, "gof_tab.unitregTMB")) {
    stop("'x' must inherit from 'gof_tab.unitregTMB'.", call. = FALSE)
  }

  if (is.null(digits)) {
    digits <- attr(x, "digits")
  }

  if (is.null(digits)) {
    digits <- 2L
  }

  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 0 || digits != as.integer(digits)) {

    stop("'digits' must be a non-negative integer.", call. = FALSE)
  }

  digits <- as.integer(digits)

  model_names <- make.unique(vapply(x, function(z) {
        if (!is.null(z$model_name) && length(z$model_name)) {
          as.character(z$model_name)[1L]
        } else {
          "unitregTMB"
        }
      },
      character(1)
    ),
    sep = " #"
  )

  tab <- rbind(LogLik = vapply(x, function(z) z$logLik, numeric(1)),
               AIC = vapply(x, function(z) z$AIC, numeric(1)),
               BIC = vapply(x, function(z) z$BIC, numeric(1)))

  colnames(tab) <- model_names

  fmt <- paste0("%.", digits, "f")

  tab_chr <- matrix(
    vapply(
      as.vector(tab),
      function(v) {
        if (length(v) &&
            is.finite(v)) {
          sprintf(fmt, v)
        } else {
          "---"
        }
      },
      character(1)
    ),
    nrow = nrow(tab),
    ncol = ncol(tab),
    dimnames = dimnames(tab)
  )

  cat("Goodness-of-fit measures for unitregTMB models:\n")

  .unitreg_print_blocks(tab = tab_chr, row_labels = rownames(tab_chr),
                        first_header = "Criterion", min_col_width = 12L)
  invisible(x)
}

#' Extract Fixed-Effect Coefficients
#'
#' Extracts regression coefficients from one distributional component or from
#' all active components.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param type Component to extract: \code{"all"}, \code{"mu"}, \code{"phi"},
#'   \code{"p0"}, or \code{"p1"}.
#' @param ... Additional arguments, currently unused.
#'
#' @return A numeric vector of coefficients.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' coef(fit, type = "mu")
#' coef(fit, type = "all")
#' }
#' @export
coef.unitregTMB <- function(
    object, type = c("all", "mu", "phi", "p0", "p1"), ...) {
  type <- match.arg(type)

  if (type != "all") return(object$model_coef[[type]])

  unlist(lapply(.unitreg_components, function(component) {
    if (.unitreg_active_component(object, component))
      object$model_coef[[component]] else NULL
  }), use.names = TRUE)
}
#' Extract the Log-Likelihood
#'
#' Returns the maximized marginal log-likelihood of a fitted
#' \code{unitregTMB} model.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param ... Additional arguments, currently unused.
#'
#' @return An object of class \code{logLik}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' logLik(fit)
#' }
#' @importFrom stats logLik
#' @export
#' @noRd
logLik.unitregTMB <- function(object, ...) {
  out <- object$logLik
  attr(out, "df") <- object$npar
  attr(out, "nobs") <- object$nobs
  class(out) <- "logLik"
  out
}
#' Akaike Information Criterion
#'
#' Computes the AIC for a fitted \code{unitregTMB} model.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param ... Additional arguments passed to \code{\link[stats]{AIC}}.
#' @param k Penalty per parameter; the default \code{2} gives the usual AIC.
#'
#' @return A numeric AIC value.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' AIC(fit)
#' }
#' @importFrom stats AIC
#' @export
AIC.unitregTMB <- function(object, ..., k = 2)
  stats::AIC(logLik(object), k = k)
#' Bayesian Information Criterion
#'
#' Computes the BIC for a fitted \code{unitregTMB} model.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param ... Additional arguments, currently unused.
#'
#' @return A numeric BIC value.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' BIC(fit)
#' }
#' @importFrom stats BIC
#' @export
BIC.unitregTMB <- function(object, ...)
  stats::BIC(logLik(object))
#' Wald Confidence Intervals
#'
#' Computes Wald confidence intervals for fixed regression coefficients in the
#' active \code{mu}, \code{phi}, \code{p0}, and \code{p1} predictors. Random-
#' effect covariance parameters and CAR hyperparameters are excluded.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param parm Optional numeric or character vector selecting coefficients.
#'   Character names use the form \code{"mu.x"}, \code{"phi.x"},
#'   \code{"p0.x"}, or \code{"p1.x"}.
#' @param level Confidence level.
#' @param ... Additional arguments, currently unused.
#'
#' @return A matrix with lower and upper confidence limits.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' confint(fit)
#' confint(fit, parm = "mu.x")
#' }
#' @export
confint.unitregTMB <- function(object, parm, level = 0.95, ...) {
  if (is.null(object$sd_report))
    stop("Refit the model with control = unitregTMB.control(sd.report = TRUE).")
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("'level' must be a single number in (0, 1).")

  est <- .unitreg_coef_values(object)
  se <- .unitreg_coef_se(object)
  nms <- .unitreg_coef_names(object)

  alpha <- 1 - level
  z <- stats::qnorm(1 - alpha / 2)
  out <- cbind(est - z * se, est + z * se)
  rownames(out) <- nms
  colnames(out) <- paste0(
    format(100 * c(alpha / 2, 1 - alpha / 2), trim = TRUE, digits = 3),
    "%"
  )

  if (!missing(parm)) out <- out[parm, , drop = FALSE]
  out
}
#' Profile a Fixed Regression Coefficient
#'
#' Computes a TMB profile likelihood for one fixed regression coefficient.
#'
#' @param fitted A fitted \code{unitregTMB} object.
#' @param parm Name of one coefficient, using the form \code{"mu.x"},
#'   \code{"phi.x"}, \code{"p0.x"}, or \code{"p1.x"}.
#' @param level.max Maximum confidence level passed to the profiling procedure.
#' @param trace Logical; print profiling progress.
#' @param ... Additional arguments passed to \code{\link[TMB]{tmbprofile}}.
#'
#' @return An object inheriting from \code{profile.unitregTMB} and
#'   \code{tmbprofile}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(y ~ x, data = dat, family = unitgamma(model_for = "mean"))
#' profile(fit, parm = "mu.x")
#' }
#' @export
profile.unitregTMB <- function(
    fitted, parm, level.max = 0.99, trace = FALSE, ...) {

  friendly <- .unitreg_coef_names(fitted)
  if (missing(parm) || length(parm) != 1L || !parm %in% friendly)
    stop("Specify one valid coefficient: ", paste(friendly, collapse = ", "))

  opt_names <- names(fitted$opt$par)
  positions <- unlist(lapply(.unitreg_components, function(component) {
    if (!.unitreg_active_component(fitted, component)) return(integer())
    which(opt_names == paste0("beta_", component))
  }), use.names = FALSE)

  if (length(positions) != length(friendly))
    stop("Could not align regression coefficients with TMB fixed parameters.")

  lc <- numeric(length(fitted$opt$par))
  lc[positions[match(parm, friendly)]] <- 1

  if (trace) message("Profiling parameter '", parm, "'.")

  if (!is.numeric(level.max) || length(level.max) != 1L ||
      level.max <= 0 || level.max >= 1)
    stop("'level.max' must be a single number in (0, 1).")

  dots <- list(...)
  if (is.null(dots$ytol))
    dots$ytol <- stats::qchisq(level.max, df = 1) / 2

  prof <- do.call(
    TMB::tmbprofile,
    c(list(obj = fitted$obj, name = parm, lincomb = lc, trace = trace), dots)
  )
  class(prof) <- c("profile.unitregTMB", class(prof))
  prof
}
#' Extract the Number of Observations
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param ... Additional arguments, currently unused.
#'
#' @return The number of observations used in the fit.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' nobs(fit)
#' }
#' @method nobs unitregTMB
#' @export
nobs.unitregTMB <- function(object, ...) object$nobs
#' Extract Residual Degrees of Freedom
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param ... Additional arguments, currently unused.
#'
#' @return The number of observations minus the number of estimated fixed
#'   parameters and variance/hyperparameters counted in \code{npar}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' df.residual(fit)
#' }
#' @export
#' @noRd
df.residual.unitregTMB <- function(object, ...) object$nobs - object$npar
#' Extract Fixed Effects
#'
#' Generic function for extracting fixed effects.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' fixed_effects(fit)
#' }
#' @export
fixed_effects <- function(object, ...) UseMethod("fixed_effects")
#' Extract Fixed Effects from a unitregTMB Model
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param ... Additional arguments, currently unused.
#'
#' @return A named list with fixed effects for the active distributional
#'   predictors.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' fixed_effects(fit)
#' }
#' @rdname fixed_effects
#' @export
fixed_effects.unitregTMB <- function(object, ...) {
  out <- object$model_coef
  class(out) <- c("fixef.unitregTMB", "list")
  out
}
#' Print Fixed Effects
#'
#' @param x An object of class \code{fixef.unitregTMB}.
#' @param digits Number of significant digits used to print estimates.
#' @param ... Additional arguments, currently unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' print(fixed_effects(fit))
#' }
#' @export
#' @noRd
print.fixef.unitregTMB <- function(
    x, digits = max(3, getOption("digits") - 3), ...) {

  for (component in .unitreg_components) {
    b <- x[[component]]
    if (is.null(b)) next
    cat("\n", .unitreg_labels[component], " coefficients:\n", sep = "")
    print.default(format(b, digits = digits), print.gap = 2L, quote = FALSE)
  }

  invisible(x)
}
#' Extract Spatial Effects
#'
#' Generic function for extracting region-specific proper CAR effects.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments.
#'
#' @return A model-specific representation of spatial effects.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' spatial_effects(fit)
#' }
#' @export
spatial_effects <- function(object, ...) UseMethod("spatial_effects")
#' Extract Spatial Effects from a unitregTMB Model
#'
#' Extracts region-specific proper CAR effects from one distributional predictor
#' or from all predictors.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param component Component to extract: \code{"mu"}, \code{"phi"},
#'   \code{"p0"}, \code{"p1"}, or \code{"all"}.
#' @param ... Additional arguments, currently unused.
#'
#' @return For a single component, a data frame containing region-specific
#'   estimates and standard errors, or \code{NULL} if no spatial effect is
#'   present. For \code{component = "all"}, a named list is returned.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' spatial_effects(fit, component = "all")
#' }
#' @rdname spatial_effects
#' @export
spatial_effects.unitregTMB <- function(
    object, component = c("mu", "phi", "p0", "p1", "all"), ...) {
  component <- match.arg(component)

  if (component == "all") {
    return(setNames(
      lapply(.unitreg_components, function(x)
        object[[paste0("spatial_effects_", x)]]),
      .unitreg_components
    ))
  }

  if (!.unitreg_has_spatial(object, component)) return(NULL)
  object[[paste0("spatial_effects_", component)]]
}
#' Extract a Fixed-Effect Design Matrix
#'
#' Extracts the fixed-effect design matrix used for one distributional
#' predictor.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param component Component to extract: \code{"mu"}, \code{"phi"},
#'   \code{"p0"}, or \code{"p1"}.
#' @param ... Additional arguments, currently unused.
#'
#' @return The corresponding fixed-effect design matrix, or \code{NULL} for an
#'   inactive boundary component.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' head(model.matrix(fit, component = "mu"))
#' }
#' @export
model.matrix.unitregTMB <- function(
    object, component = c("mu", "phi", "p0", "p1"), ...) {
  component <- match.arg(component)
  if (!.unitreg_active_component(object, component)) return(NULL)

  X <- object$obj$env$data[[paste0("X_", component)]]
  if (is.null(X)) return(NULL)

  nm <- names(object$model_coef[[component]])
  if (length(nm) == ncol(X)) colnames(X) <- nm
  X
}
#' Extract a Model Formula
#'
#' Extracts the formula used for one distributional predictor.
#'
#' @param x A fitted \code{unitregTMB} object.
#' @param component Component to extract: \code{"mu"}, \code{"phi"},
#'   \code{"p0"}, or \code{"p1"}.
#' @param ... Additional arguments, currently unused.
#'
#' @return A formula, or \code{NULL} for an inactive boundary component.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' formula(fit)
#' formula(fit, component = "phi")
#' }
#' @export
#' @noRd
formula.unitregTMB <- function( x, component = c("mu", "phi", "p0", "p1"), ...) {
        component <- match.arg(component)
        switch(component,
               mu = x$formula,
               phi = x$phi.formula,
               p0 = x$p0.formula,
               p1 = x$p1.formula)
}
#' Extract Fitted Values
#'
#' Reconstructs fitted values from the optimized fixed effects, conditional
#' modes of ordinary random effects, and proper CAR effects in all active
#' distributional predictors.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param type Quantity to return: \code{"response"}, \code{"mu"},
#'   \code{"phi"}, \code{"p0"}, or \code{"p1"}.
#' @param ... Additional arguments, currently unused.
#'
#' @details
#' The interpretation of \code{type = "response"} depends on the
#' parameterization of the continuous component.
#'
#' For mean-parameterized families, the fitted response is
#' \deqn{
#'   p_1 + (1 - p_0 - p_1)\mu,
#' }
#' which is the unconditional conditional mean of the augmented response.
#'
#' For modal and quantile parameterizations, \code{mu} represents the modeled
#' mode or quantile of the continuous component rather than its mean.
#' Consequently, \code{type = "response"} returns the modeled response-scale
#' \code{mu} parameter instead of interpreting
#' \eqn{p_1 + (1-p_0-p_1)\mu} as an expectation.
#'
#' The options \code{"mu"}, \code{"phi"}, \code{"p0"}, and \code{"p1"}
#' always return the corresponding fitted distributional parameter.
#'
#' @return A numeric vector of fitted values.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' head(fitted(fit, type = "mu"))
#' head(fitted(fit, type = "phi"))
#' head(fitted(fit, type = "response"))
#' }
#' @export
fitted.unitregTMB <- function(object, type = c("response", "mu", "phi", "p0", "p1"), ...) {

  type <- match.arg(type)
  dat <- object$obj$env$data
  par_full <- .unitreg_full_par(object)
  n <- nrow(dat$X_mu)

  eta <- function(component, offset = NULL) {
    X <- dat[[paste0("X_", component)]]
    b <- object$model_coef[[component]]
    out <- if (!is.null(b) && length(b)) as.vector(X %*% b) else rep(0, n)

    if (!is.null(offset)) out <- out + as.numeric(offset)

    if (.unitreg_has_re(object, component)) {
      Z <- dat[[paste0("Z_", component)]]
      u <- par_full[[paste0("re_coef_", component)]]
      if (is.null(Z) || is.null(u)) stop("Random-effect information is missing for component '", component, "'.", call. = FALSE)
      if (length(u) != ncol(Z)) stop("Dimension mismatch for random effects in component '", component, "'.", call. = FALSE)
      out <- out + as.vector(Z %*% u)
    }

    if (.unitreg_has_spatial(object, component)) {
      Zs <- dat[[paste0("Z_spatial_", component)]]
      us <- par_full[[paste0("u_spatial_", component)]]
      if (is.null(Zs) || is.null(us)) stop("Spatial-effect information is missing for component '", component, "'.", call. = FALSE)
      if (length(us) != ncol(Zs)) stop("Dimension mismatch for spatial effects in component '", component, "'.", call. = FALSE)
      out <- out + as.vector(Zs %*% us)
    }

    out
  }

  eta_mu <- eta("mu", dat$offset_mu)
  mu <- object$link.mu$linkinv(eta_mu)
  if (type == "mu") return(as.numeric(mu))

  eta_phi <- eta("phi")
  phi_link <- object$family_object$phi_link_r_name %||% "log"
  phi <- stats::make.link(phi_link)$linkinv(eta_phi)
  if (type == "phi") return(as.numeric(phi))

  eta_p0 <- if (isTRUE(object$has_p0)) eta("p0") else rep(0, n)
  eta_p1 <- if (isTRUE(object$has_p1)) eta("p1") else rep(0, n)

  if (isTRUE(object$has_p0) && isTRUE(object$has_p1)) {
    m <- pmax(0, eta_p0, eta_p1)
    e0 <- exp(eta_p0 - m)
    e1 <- exp(eta_p1 - m)
    ec <- exp(-m)
    den <- ec + e0 + e1
    p0 <- e0 / den
    p1 <- e1 / den
  } else if (isTRUE(object$has_p0)) {
    p0 <- stats::plogis(eta_p0)
    p1 <- rep(0, n)
  } else if (isTRUE(object$has_p1)) {
    p0 <- rep(0, n)
    p1 <- stats::plogis(eta_p1)
  } else {
    p0 <- rep(0, n)
    p1 <- rep(0, n)
  }

  if (type == "p0") return(as.numeric(p0))
  if (type == "p1") return(as.numeric(p1))

  if (identical(.unitreg_model_for(object), "mean")) {
    return(as.numeric(p1 + (1 - p0 - p1) * mu))
  }

  as.numeric(mu)
}
#' Extract the Fixed-Effect Variance-Covariance Matrix
#'
#' Extracts the covariance matrix for fixed regression coefficients only.
#' Ordinary random-effect covariance parameters and CAR hyperparameters are
#' excluded.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param ... Additional arguments, currently unused.
#'
#' @return A covariance matrix with component-prefixed coefficient names.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' vcov(fit)
#' }
#' @export
vcov.unitregTMB <- function(object, ...) {
  if (!is.null(object$sd_report) && !is.null(object$sd_report$cov.fixed)) {
    V <- object$sd_report$cov.fixed
  } else if (!is.null(object$opt$hessian)) {
    V <- tryCatch(solve(object$opt$hessian), error = function(e) NULL)
    if (is.null(V)) stop("The Hessian is singular; covariance matrix unavailable.")
  } else {
    stop("Refit the model with control = unitregTMB.control(sd.report = TRUE).")
  }

  tmb_names <- rownames(V)
  if (is.null(tmb_names)) tmb_names <- names(object$opt$par)

  idx <- unlist(lapply(.unitreg_components, function(component) {
    if (!.unitreg_active_component(object, component)) return(integer())
    which(tmb_names == paste0("beta_", component))
  }), use.names = FALSE)

  V <- V[idx, idx, drop = FALSE]
  nm <- .unitreg_coef_names(object)

  if (length(nm) != nrow(V))
    stop("Could not align fixed coefficients with the covariance matrix.")

  dimnames(V) <- list(nm, nm)
  V
}
#' Update and Refit a unitregTMB Model
#'
#' Updates the original model call and, by default, refits the model. Fixed,
#' random, inflation, and spatial specifications can be changed through
#' \code{...}.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param formula. Optional update to the \code{mu} formula, following
#'   \code{\link[stats]{update.formula}} syntax.
#' @param ... Updated call arguments, such as \code{phi.formula},
#'   \code{p0.formula}, \code{p1.formula}, \code{spatial},
#'   \code{spatial.phi}, \code{spatial.p0}, \code{spatial.p1},
#'   \code{adjacency}, \code{family}, \code{data}, or \code{control}.
#' @param evaluate Logical; if \code{TRUE}, refit the model. Otherwise return
#'   the updated call.
#'
#' @return A fitted \code{unitregTMB} object or an unevaluated call.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x, data = dat, family = unitgamma(model_for = "mean")
#' )
#' update(fit, . ~ . + I(x^2), evaluate = FALSE)
#' }
#' @export
#' @noRd
update.unitregTMB <- function(object, formula., ..., evaluate = TRUE) {
  if (is.null(object$call)) stop("The model object does not contain its original call.")

  call <- object$call
  extras <- match.call(expand.dots = FALSE)$...

  if (!missing(formula.))
    call$formula <- stats::update.formula(object$formula, formula.)

  if (length(extras)) {
    for (nm in names(extras)) call[[nm]] <- extras[[nm]]
  }

  if (evaluate) eval(call, parent.frame()) else call
}


`%||%` <- function(x, y) if (is.null(x)) y else x
