#' Compare Fixed-Effect Coefficients from unitregTMB Models
#'
#' @description
#' Extracts fixed-effect estimates and standard errors from one or more fitted
#' \code{unitregTMB} models and displays them in a compact comparison table.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param ... Additional fitted \code{unitregTMB} objects.
#' @param component Character vector selecting distributional components:
#'   \code{"mu"}, \code{"phi"}, \code{"p0"}, \code{"p1"}, or \code{"all"}.
#'   The default is \code{"mu"}.
#' @param digits Non-negative integer giving the number of decimal places.
#' @param print_tex Logical. If \code{TRUE}, prints a LaTeX table.
#'   The default is \code{FALSE}.
#' @param signif_stars Logical. If \code{TRUE}, adds Wald-test significance
#'   symbols to the coefficient estimates. The default is \code{FALSE}.
#'
#' @details
#' The method supports the distributional predictors \code{mu}, \code{phi},
#' \code{p0}, and \code{p1}. Coefficients are distinguished internally by
#' component; the \code{mu.} prefix is omitted when printed.
#'
#' Two-sided asymptotic Wald p-values are computed as
#' \code{2 * pnorm(-abs(estimate / se))} using full-precision estimates and
#' standard errors stored in the fitted models. Invalid or nonpositive
#' standard errors produce missing p-values and no significance symbols.
#' These are unadjusted p-values; no multiple-testing correction is applied.
#' Significance codes: \code{***} (< 0.001), \code{**} (< 0.01),
#' \code{*} (< 0.05), \code{.} (< 0.1).
#'
#' In console output, significance codes use a fixed three-character
#' field so that estimates and standard errors remain vertically aligned.
#' The console printer automatically splits wide tables into blocks according
#' to \code{getOption("width")}. LaTeX output is a simple \code{table} and
#' \code{tabular} environment.
#'
#' @return An object of class \code{summary_coef.unitregTMB}, containing one
#'   element per fitted model with \code{model_name}, \code{coefficients},
#'   \code{std_errors}, \code{p_values}, and \code{coef_names}. The original
#'   structure and default printed output are retained; the \code{p_values}
#'   member is an additive extension.
#'
#' @examples
#' \donttest{
#' data("wqi", package = "unitregTMB")
#' fit_beta <- unitregTMB(y ~ location, data = wqi,
#'                        family = beta_fam(model_for = "mean"))
#' fit_ug <- unitregTMB(y ~ location, data = wqi,
#'                      family = unitgamma(model_for = "mean"))
#' summary_coef(fit_beta, fit_ug)
#' summary_coef(fit_beta, fit_ug, component = "all", signif_stars = TRUE)
#' summary_coef(fit_beta, fit_ug, signif_stars = TRUE, print_tex = TRUE)
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom regcore summary_coef
#' @method summary_coef unitregTMB
#' @rdname summary_coef_unitregTMB
#' @export
summary_coef.unitregTMB <- function(object, ..., component = "mu", digits = 3,
                                    print_tex = FALSE, signif_stars = FALSE) {
  models <- list(object, ...)
  
  if (!length(models)) stop("At least one fitted 'unitregTMB' model must be supplied.", call. = FALSE)
  if (!all(vapply(models, inherits, logical(1), what = "unitregTMB")))
    stop("All supplied objects must inherit from 'unitregTMB'.", call. = FALSE)
  if (!is.character(component) || !length(component) || anyNA(component))
    stop("'component' must be a non-empty character vector.", call. = FALSE)
  
  component <- unique(tolower(component))
  components <- c("mu", "phi", "p0", "p1")
  valid <- c("all", components)
  if (!all(component %in% valid))
    stop("'component' must contain only: ", paste(sprintf("'%s'", valid), collapse = ", "), ".", call. = FALSE)
  
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 0 || digits != as.integer(digits))
    stop("'digits' must be a non-negative integer.", call. = FALSE)
  
  for (arg in c("print_tex", "signif_stars")) {
    value <- get(arg)
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop("'", arg, "' must be TRUE or FALSE.", call. = FALSE)
  }
  
  digits <- as.integer(digits)
  selected <- if ("all" %in% component) components else component
  
  model_labels <- make.unique(vapply(models, .unitreg_compact_model_label, character(1)), sep = " #")
  
  results <- lapply(seq_along(models), function(i) {
    model <- models[[i]]
    b_all <- numeric()
    se_all <- numeric()
    p_all <- numeric()
    
    for (comp in selected) {
      b <- model$model_coef[[comp]]
      if (is.null(b) || !length(b)) next
      
      se <- model$model_std[[paste0("std.", comp)]]
      if (is.null(se)) se <- rep(NA_real_, length(b))
      if (length(se) != length(b))
        stop("Coefficient and standard-error lengths differ for component '", comp,
             "' in model '", model_labels[i], "'.", call. = FALSE)
      
      nm0 <- names(b)
      if (is.null(nm0) || length(nm0) != length(b)) nm0 <- paste0("coef", seq_along(b))
      nm <- paste0(comp, ".", nm0)
      names(b) <- nm
      names(se) <- nm
      
      p <- rep(NA_real_, length(b))
      valid_se <- is.finite(b) & is.finite(se) & se > 0
      p[valid_se] <- 2 * stats::pnorm(-abs(b[valid_se] / se[valid_se]))
      names(p) <- nm
      
      b_all <- c(b_all, b)
      se_all <- c(se_all, se)
      p_all <- c(p_all, p)
    }
    
    list(model_name = model_labels[i], coefficients = b_all,
         std_errors = se_all, p_values = p_all, coef_names = names(b_all))
  })
  
  class(results) <- c("summary_coef.unitregTMB", "list")
  attr(results, "digits") <- digits
  attr(results, "print_tex") <- print_tex
  attr(results, "component") <- selected
  attr(results, "signif_stars") <- signif_stars
  results
}

#' Print a summary_coef.unitregTMB Object
#'
#' @description
#' Prints coefficient estimates and standard errors returned by
#' \code{\link{summary_coef.unitregTMB}}, optionally with Wald-test stars.
#'
#' @param x An object of class \code{summary_coef.unitregTMB}.
#' @param digits Optional non-negative integer overriding the value in \code{x}.
#' @param print_tex Optional logical overriding the output format in \code{x}.
#' @param signif_stars Optional logical overriding the significance display in
#'   \code{x}. When \code{TRUE}, adds significance codes after (SE).
#' @param ... Additional arguments, currently unused.
#'
#' @return Invisibly returns \code{x}.
#' @author Ricardo Rasmussen Petterle
#' @method print summary_coef.unitregTMB
#' @rdname print_summary_coef_unitregTMB
#' @export
print.summary_coef.unitregTMB <- function(x, digits = NULL, print_tex = NULL,
                                          signif_stars = NULL, ...) {
  if (!inherits(x, "summary_coef.unitregTMB"))
    stop("'x' must inherit from 'summary_coef.unitregTMB'.", call. = FALSE)
  
  if (is.null(digits)) digits <- attr(x, "digits")
  if (is.null(digits)) digits <- 3L
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 0 || digits != as.integer(digits))
    stop("'digits' must be a non-negative integer.", call. = FALSE)
  
  if (is.null(print_tex)) print_tex <- isTRUE(attr(x, "print_tex"))
  if (is.null(signif_stars)) signif_stars <- isTRUE(attr(x, "signif_stars"))
  for (arg in c("print_tex", "signif_stars")) {
    value <- get(arg)
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop("'", arg, "' must be TRUE or FALSE.", call. = FALSE)
  }
  
  digits <- as.integer(digits)
  model_names <- make.unique(vapply(x, function(z) {
    if (!is.null(z$model_name) && length(z$model_name)) as.character(z$model_name)[1L] else "unitregTMB"
  }, character(1)), sep = " #")
  
  coef_names <- Reduce(union, lapply(x, function(z) {
    if (is.null(z$coef_names)) character() else z$coef_names
  }))
  if (!length(coef_names)) {
    cat("No fixed-effect coefficients to display.\n")
    return(invisible(x))
  }
  
  tab <- matrix("---", nrow = length(coef_names), ncol = length(x),
                dimnames = list(coef_names, model_names))
  fmt <- paste0("%.", digits, "f")
  
  .wald_stars <- function(p) {
    if (length(p) != 1L || !is.finite(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    if (p < 0.1) return(".")
    ""
  }
  
  for (j in seq_along(x)) {
    z <- x[[j]]
    if (!length(z$coefficients)) next
    
    for (nm in names(z$coefficients)) {
      est <- unname(z$coefficients[nm])
      se <- unname(z$std_errors[nm])
      cell <- if (length(se) && is.finite(se)) {
        paste0(sprintf(fmt, est), " (", sprintf(fmt, se), ")")
      } else {
        paste0(sprintf(fmt, est), " (---)")
      }
      
      if (signif_stars) {
        p <- if (!is.null(z$p_values)) unname(z$p_values[nm]) else NA_real_
        if (length(p) != 1L || is.na(p)) {
          p <- if (length(se) && is.finite(est) && is.finite(se) && se > 0)
            2 * stats::pnorm(-abs(est / se)) else NA_real_
        }
        stars <- .wald_stars(p)
        if (print_tex) {
          if (nzchar(stars)) cell <- paste0(cell, "\\textsuperscript{", stars, "}")
        } else {
          cell <- paste0(cell, sprintf("%-3s", stars))
        }
      }
      tab[nm, j] <- cell
    }
  }
  
  display_names <- sub("^mu\\.", "", rownames(tab))
  
  if (print_tex) {
    escape_tex <- function(z) {
      z <- gsub("_", "\\\\_", z, fixed = TRUE)
      z <- gsub("%", "\\\\%", z, fixed = TRUE)
      z <- gsub("&", "\\\\&", z, fixed = TRUE)
      z <- gsub("#", "\\\\#", z, fixed = TRUE)
      z <- gsub("$", "\\\\$", z, fixed = TRUE)
      z
    }
    
    cat("\\begin{table}[ht]\n")
    cat("\\centering\n")
    cat("\\begin{tabular}{", paste0("l", paste(rep("c", ncol(tab)), collapse = "")), "}\n", sep = "")
    cat("\\hline\n")
    cat("Predictor & ", paste(escape_tex(colnames(tab)), collapse = " & "), " \\\\\n", sep = "")
    cat("\\hline\n")
    for (i in seq_len(nrow(tab)))
      cat(escape_tex(display_names[i]), " & ", paste(tab[i, ], collapse = " & "), " \\\\\n", sep = "")
    cat("\\hline\n")
    cat("\\end{tabular}\n")
    if (signif_stars) {
      cat("\\par\\smallskip\n{\\footnotesize\\emph{Note:} ")
      notes <- character()
      if (signif_stars)
        notes <- c(notes, "Two-sided Wald tests: \\textsuperscript{***} $p<0.001$, \\textsuperscript{**} $p<0.01$, \\textsuperscript{*} $p<0.05$, \\textsuperscript{.} $p<0.10$ (unadjusted).")
      cat(paste(notes, collapse = " "), "}\n", sep = "")
    }
    cat("\\end{table}\n")
  } else {
    cat("Regression coefficients (standard errors) ", "for unitregTMB models:\n", sep = "")
    .unitreg_print_blocks(tab = tab, row_labels = display_names,
                          first_header = "Predictor", min_col_width = 15L)
    if (signif_stars) {
      cat("---\n")
      cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
    }
  }
  
  invisible(x)
}