#' Likelihood Ratio Tests for unitregTMB Models
#'
#' @description
#' Compares two or more fitted \code{unitregTMB} models by likelihood ratio tests
#' and returns an analysis-of-deviance table with model degrees of freedom,
#' information criteria, log-likelihoods, deviances, likelihood ratio statistics,
#' degrees-of-freedom differences, and p-values.
#'
#' @param object A fitted object of class \code{unitregTMB}.
#' @param ... Additional fitted \code{unitregTMB} objects to be compared with
#'   \code{object}.
#' @param test Character string specifying the reference test. Currently only
#'   \code{"Chisq"} is supported.
#'
#' @details
#' Models are ordered by their number of estimated parameters before consecutive
#' likelihood ratio comparisons are computed. The models must use the same
#' response values and observation order, the same prior weights and offset, and
#' the same distributional family, parameterization, link, and quantile level
#' when applicable.
#'
#' The usual chi-square reference distribution assumes nested models and regular
#' likelihood conditions. Comparisons that add or remove random effects,
#' spatial effects, or zero/one-inflation components can involve parameters on
#' the boundary of the parameter space. In such cases the standard chi-square
#' p-value can be non-standard and should be interpreted with caution.
#'
#' @return An object of class \code{anova.unitregTMB} inheriting from
#'   \code{data.frame}. The table contains columns \code{Df}, \code{AIC},
#'   \code{BIC}, \code{logLik}, \code{deviance}, \code{Chisq},
#'   \code{Chi Df}, and \code{Pr(>Chisq)}.
#'
#' @examples
#' \donttest{
#' data("wqi", package = "unitregTMB")
#'
#' fit0 <- unitregTMB(y ~ 1,
#'                    family = unitgamma(model_for = "mean"),
#'                    data = wqi)
#'
#' fit1 <- unitregTMB(y ~ location,
#'                    family = unitgamma(model_for = "mean"),
#'                    data = wqi)
#'
#' tab <- anova(fit0, fit1)
#' print(tab)
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @method anova unitregTMB
#' @export
anova.unitregTMB <- function(object, ..., test = c("Chisq")) {
  test <- match.arg(test)
  models <- list(object, ...)

  if (length(models) < 2L) {
    stop("At least two models must be supplied.", call. = FALSE)
  }

  ok <- vapply(models, inherits, logical(1), what = "unitregTMB")
  if (!all(ok)) {
    stop("All supplied objects must inherit from 'unitregTMB'.", call. = FALSE)
  }

  formula_text <- function(x) {
    if (is.null(x)) return(NULL)
    paste(deparse(x, width.cutoff = 500L), collapse = "")
  }

  family_signature <- function(mod) {
    fam <- mod$family_object
    code <- if (!is.null(fam$family_code)) as.integer(fam$family_code)[1L] else NA_integer_
    link <- if (!is.null(fam$link_code)) as.integer(fam$link_code)[1L] else {
      if (!is.null(mod$link.mu$name)) as.character(mod$link.mu$name)[1L] else NA_character_
    }
    tau <- if (is.null(mod$tau)) NA_real_ else as.numeric(mod$tau)[1L]
    list(
      family = as.character(mod$family)[1L],
      code = code,
      link = link,
      tau = tau
    )
  }

  same_signature <- function(a, b) {
    identical(a$family, b$family) &&
      identical(a$code, b$code) &&
      identical(a$link, b$link) &&
      isTRUE(all.equal(a$tau, b$tau, tolerance = 0, check.attributes = FALSE))
  }

  model_spec <- function(mod) {
    forms <- list(
      mu = mod$formula,
      phi = mod$phi.formula,
      p0 = mod$p0.formula,
      p1 = mod$p1.formula,
      spatial.mu = mod$spatial.formula,
      spatial.phi = mod$spatial.phi.formula,
      spatial.p0 = mod$spatial.p0.formula,
      spatial.p1 = mod$spatial.p1.formula
    )
    txt <- lapply(forms, formula_text)
    keep <- !vapply(txt, is.null, logical(1))
    txt <- txt[keep]
    if (!length(txt)) return("<formula unavailable>")
    paste(paste0(names(txt), ": ", unlist(txt, use.names = FALSE)), collapse = "; ")
  }

  random_signature <- function(mod) {
    forms <- list(mu = mod$formula, phi = mod$phi.formula, p0 = mod$p0.formula, p1 = mod$p1.formula)
    out <- vapply(names(forms), function(nm) {
      f <- forms[[nm]]
      if (is.null(f)) return("")
      bars <- reformulas::findbars(f)
      if (!length(bars)) return("")
      paste(vapply(bars, function(b) paste(deparse(b, width.cutoff = 500L), collapse = ""), character(1)), collapse = " + ")
    }, character(1))
    paste(out, collapse = "|")
  }

  boundary_signature <- function(mod) {
    spatial <- c(
      formula_text(mod$spatial.formula),
      formula_text(mod$spatial.phi.formula),
      formula_text(mod$spatial.p0.formula),
      formula_text(mod$spatial.p1.formula)
    )
    spatial <- spatial[nzchar(spatial)]
    paste(
      random_signature(mod),
      paste(spatial, collapse = "|"),
      isTRUE(mod$has_p0),
      isTRUE(mod$has_p1),
      sep = "||"
    )
  }

  ref <- models[[1L]]
  ref_family <- family_signature(ref)
  ref_rows <- if (!is.null(ref$model)) rownames(ref$model) else NULL

  for (i in seq_along(models)[-1L]) {
    mod <- models[[i]]

    if (!same_signature(ref_family, family_signature(mod))) {
      stop(
        "Likelihood-ratio comparisons require the same distribution, parameterization, link, and quantile level.",
        call. = FALSE
      )
    }

    if (!identical(as.integer(ref$nobs), as.integer(mod$nobs))) {
      stop("Models were fitted using different numbers of observations.", call. = FALSE)
    }

    if (!isTRUE(all.equal(ref$Y, mod$Y, tolerance = 0, check.attributes = FALSE))) {
      stop("Models were fitted to different response values or observation orders.", call. = FALSE)
    }

    if (!is.null(ref_rows) && !is.null(mod$model) && !identical(ref_rows, rownames(mod$model))) {
      stop("Models were fitted to different sets or orders of observations.", call. = FALSE)
    }

    if (!isTRUE(all.equal(ref$weights, mod$weights, tolerance = 0, check.attributes = FALSE))) {
      stop("Models use different prior weights and cannot be compared by an LRT.", call. = FALSE)
    }

    if (!isTRUE(all.equal(ref$offset, mod$offset, tolerance = 0, check.attributes = FALSE))) {
      stop("Models use different offsets and cannot be compared by an LRT.", call. = FALSE)
    }
  }

  info <- lapply(models, function(mod) {
    ll <- logLik(mod)
    df <- attr(ll, "df")
    if (is.null(df) || !length(df) || !is.finite(df)) df <- mod$npar

    nobs_ll <- attr(ll, "nobs")
    if (is.null(nobs_ll) || !length(nobs_ll) || !is.finite(nobs_ll)) nobs_ll <- mod$nobs

    aic <- tryCatch(AIC(mod), error = function(e) mod$AIC)
    bic <- tryCatch(BIC(mod), error = function(e) mod$BIC)

    data.frame(
      Df = as.integer(df),
      AIC = as.numeric(aic),
      BIC = as.numeric(bic),
      logLik = as.numeric(ll),
      deviance = -2 * as.numeric(ll),
      nobs = as.integer(nobs_ll),
      model = model_spec(mod),
      stringsAsFactors = FALSE
    )
  })

  tab <- do.call(rbind, info)

  if (length(unique(tab$nobs)) > 1L) {
    stop("Models were fitted using different numbers of observations.", call. = FALSE)
  }

  model_names <- as.character(match.call(expand.dots = TRUE))[-1L]
  model_names <- make.unique(model_names[seq_along(models)])

  ord <- order(tab$Df, seq_len(nrow(tab)))
  tab <- tab[ord, , drop = FALSE]
  models <- models[ord]
  model_names <- model_names[ord]

  chisq <- c(NA_real_, 2 * diff(tab$logLik))
  chi_df <- c(NA_integer_, diff(tab$Df))
  pval <- rep(NA_real_, nrow(tab))

  valid <- seq_len(nrow(tab)) > 1L &
    is.finite(chisq) &
    is.finite(chi_df) &
    chi_df > 0L &
    chisq >= 0

  pval[valid] <- pchisq(chisq[valid], df = chi_df[valid], lower.tail = FALSE)

  if (any(chi_df[-1L] <= 0L, na.rm = TRUE)) {
    warning(
      "At least one consecutive comparison has no increase in model degrees of freedom; its LRT p-value was not computed.",
      call. = FALSE
    )
  }

  if (any(chisq[-1L] < 0, na.rm = TRUE)) {
    warning(
      "At least one model with more parameters has a lower log-likelihood; check nesting and convergence. Its LRT p-value was not computed.",
      call. = FALSE
    )
  }

  specs <- tab$model
  tab$model <- NULL
  tab$nobs <- NULL
  tab$Chisq <- chisq
  tab$`Chi Df` <- chi_df
  tab$`Pr(>Chisq)` <- pval
  rownames(tab) <- model_names

  attr(tab, "heading") <- "Likelihood ratio tests of unitregTMB models:"
  attr(tab, "models") <- specs
  attr(tab, "model_names") <- model_names
  class(tab) <- c("anova.unitregTMB", "data.frame")
  tab
}


#' Print an anova.unitregTMB Object
#'
#' @description
#' Prints the analysis-of-deviance table produced by
#' \code{\link{anova.unitregTMB}}, including the model specifications and
#' optional significance codes.
#'
#' @param x An object of class \code{anova.unitregTMB}.
#' @param ... Additional arguments passed to or from methods.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#' data("wqi", package = "unitregTMB")
#'
#' fit0 <- unitregTMB(y ~ 1,
#'                    family = unitgamma(model_for = "mean"),
#'                    data = wqi)
#'
#' fit1 <- unitregTMB(y ~ location,
#'                    family = unitgamma(model_for = "mean"),
#'                    data = wqi)
#'
#' x <- anova(fit0, fit1)
#' print(x)
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @method print anova.unitregTMB
#' @export
#' @noRd
print.anova.unitregTMB <- function(x, ...) {
  heading <- attr(x, "heading")
  model_names <- attr(x, "model_names")
  models <- attr(x, "models")

  if (!is.null(heading)) cat(heading, "\n")

  if (!is.null(model_names) && !is.null(models)) {
    cat("Models:\n")
    for (i in seq_along(model_names)) {
      cat(model_names[i], ": ", models[i], "\n", sep = "")
    }
    cat("\n")
  }

  df_to_print <- x
  attr(df_to_print, "heading") <- NULL
  attr(df_to_print, "models") <- NULL
  attr(df_to_print, "model_names") <- NULL
  class(df_to_print) <- "data.frame"

  pvals <- x$`Pr(>Chisq)`
  chi_df <- x$`Chi Df`

  df_to_print$Df <- as.integer(df_to_print$Df)
  df_to_print$AIC <- format(round(df_to_print$AIC, 2), nsmall = 2)
  df_to_print$BIC <- format(round(df_to_print$BIC, 2), nsmall = 2)
  df_to_print$logLik <- format(round(df_to_print$logLik, 2), nsmall = 2)
  df_to_print$deviance <- format(round(df_to_print$deviance, 2), nsmall = 2)
  df_to_print$Chisq <- ifelse(is.na(x$Chisq), "", format(round(x$Chisq, 4), nsmall = 4))
  df_to_print$`Chi Df` <- ifelse(is.na(chi_df), "", format(as.integer(chi_df), trim = TRUE))
  df_to_print$`Pr(>Chisq)` <- ifelse(
    is.na(pvals),
    "",
    format(format.pval(pvals, digits = 4, eps = 1e-4), justify = "right")
  )

  has_pvals <- any(!is.na(pvals))
  show_stars <- getOption("show.signif.stars", TRUE)

  if (has_pvals && show_stars) {
    signif <- symnum(
      pvals,
      corr = FALSE,
      na = FALSE,
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols = c("***", "**", "*", ".", " ")
    )
    df_to_print$Signif <- as.character(signif)
    df_to_print$Signif[is.na(pvals)] <- ""
    names(df_to_print)[ncol(df_to_print)] <- ""
  }

  print(df_to_print, right = TRUE, quote = FALSE)

  if (has_pvals && show_stars) {
    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }

  invisible(x)
}
