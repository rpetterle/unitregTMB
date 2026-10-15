#' Compare Fixed-Effect Coefficients from unitregTMB Models
#'
#' @description
#' Extracts fixed-effect estimates and standard errors from one or more fitted
#' \code{unitregTMB} models and displays them in a compact comparison table.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param ... Additional fitted \code{unitregTMB} objects.
#' @param component Character vector selecting the distributional components
#'   to be displayed. Available options are \code{"mu"}, \code{"phi"},
#'   \code{"p0"}, \code{"p1"}, and \code{"all"}. The default is \code{"mu"}.
#' @param digits Non-negative integer giving the number of decimal places used
#'   when printing estimates and standard errors.
#' @param print_tex Logical. If \code{TRUE}, a simple LaTeX table is printed.
#'   The default is \code{FALSE}.
#'
#' @details
#' The method supports the four distributional predictors available in
#' \code{unitregTMB}: the location predictor \code{mu}, the precision or
#' dispersion predictor \code{phi}, and the boundary-mass predictors
#' \code{p0} and \code{p1}.
#'
#' Coefficients are internally identified by both their distributional
#' component and coefficient name. For display purposes, the \code{mu.}
#' prefix is omitted, whereas the \code{phi.}, \code{p0.}, and \code{p1.}
#' prefixes are retained.
#'
#' Model labels are automatically constructed from the fitted family and
#' parameterization. Mean and modal models are displayed, for example, as
#' \code{"Beta (mean)"} and \code{"Beta (mode)"}. For quantile regression,
#' the quantile level is used directly in the label, for example,
#' \code{"Kumaraswamy (tau = 0.5)"}, without an additional
#' \code{"quantile"} tag.
#'
#' When several models are compared and the resulting table exceeds the
#' current console width, the print method automatically splits the output
#' into blocks. The available width is determined by \code{getOption("width")}.
#'
#' Standard errors are taken from the fitted \code{unitregTMB} object.
#' If a standard error is unavailable, the estimate is retained and its
#' standard error is displayed as \code{---}.
#'
#' @return An object of class \code{summary_coef.unitregTMB}. It is a list
#'   with one element per fitted model, containing the model label,
#'   coefficient estimates, standard errors, and coefficient names.
#'
#' @examples
#' \donttest{
#' data("wqi", package = "unitregTMB")
#'
#' fit_beta <- unitregTMB(
#'   y ~ location,
#'   data = wqi,
#'   family = beta_fam(model_for = "mean")
#' )
#'
#' fit_ug <- unitregTMB(
#'   y ~ location,
#'   data = wqi,
#'   family = unitgamma(model_for = "mean")
#' )
#'
#' summary_coef(
#'   fit_beta,
#'   fit_ug,
#'   component = "mu",
#'   digits = 2
#' )
#'
#' ## Display both location and dispersion predictors
#' summary_coef(
#'   fit_beta,
#'   fit_ug,
#'   component = c("mu", "phi"),
#'   digits = 3
#' )
#'
#' ## Display all available distributional predictors
#' summary_coef(
#'   fit_beta,
#'   fit_ug,
#'   component = "all",
#'   digits = 2
#' )
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom regcore summary_coef
#' @method summary_coef unitregTMB
#' @rdname summary_coef_unitregTMB
#' @export
summary_coef.unitregTMB <- function(object, ..., component = "mu", digits = 3, print_tex = FALSE) {

  models <- list(object, ...)

  if (!length(models)) {
    stop(
      "At least one fitted 'unitregTMB' model must be supplied.",
      call. = FALSE
    )
  }

  ok <- vapply(models, inherits, logical(1), what = "unitregTMB")

  if (!all(ok)) {
    stop(
      "All supplied objects must inherit from 'unitregTMB'.",
      call. = FALSE
    )
  }

  if (!is.character(component) ||
      !length(component) ||
      anyNA(component)) {
    stop(
      "'component' must be a non-empty character vector.",
      call. = FALSE
    )
  }

  component <- unique(tolower(component))

  components <- c("mu", "phi", "p0", "p1")
  valid <- c("all", components)

  if (!all(component %in% valid)) {
    stop("'component' must contain only: ", paste(sprintf("'%s'", valid), collapse = ", "), ".", call. = FALSE)
  }

  if (!is.numeric(digits) ||
      length(digits) != 1L ||
      !is.finite(digits) ||
      digits < 0 ||
      digits != as.integer(digits)) {
    stop(
      "'digits' must be a non-negative integer.",
      call. = FALSE
    )
  }

  if (length(print_tex) != 1L ||
      !is.logical(print_tex) ||
      is.na(print_tex)) {
    stop(
      "'print_tex' must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  digits <- as.integer(digits)

  selected <- if ("all" %in% component) {
    components
  } else {
    component
  }

  model_labels <- make.unique(
    vapply(
      models,
      .unitreg_compact_model_label,
      character(1)
    ),
    sep = " #"
  )

  results <- lapply(seq_along(models), function(i) {
    model <- models[[i]]

    b_all <- numeric()
    se_all <- numeric()

    for (comp in selected) {
      b <- model$model_coef[[comp]]

      if (is.null(b) || !length(b)) next

      se <- model$model_std[[paste0("std.", comp)]]

      if (is.null(se)) {
        se <- rep(NA_real_, length(b))
      }

      if (length(se) != length(b)) {
        stop(
          "Coefficient and standard-error lengths differ for component '",
          comp,
          "' in model '",
          model_labels[i],
          "'.",
          call. = FALSE
        )
      }

      nm0 <- names(b)

      if (is.null(nm0) || length(nm0) != length(b)) {
        nm0 <- paste0("coef", seq_along(b))
      }

      nm <- paste0(comp, ".", nm0)

      names(b) <- nm
      names(se) <- nm

      b_all <- c(b_all, b)
      se_all <- c(se_all, se)
    }

    list(
      model_name = model_labels[i],
      coefficients = b_all,
      std_errors = se_all,
      coef_names = names(b_all)
    )
  })

  class(results) <- c("summary_coef.unitregTMB", "list")

  attr(results, "digits") <- digits
  attr(results, "print_tex") <- print_tex
  attr(results, "component") <- selected

  results
}


#' Print a summary_coef.unitregTMB Object
#'
#' @description
#' Prints coefficient estimates and standard errors returned by
#' \code{\link{summary_coef.unitregTMB}}.
#'
#' @param x An object of class \code{summary_coef.unitregTMB}.
#' @param digits Optional non-negative integer overriding the number of
#'   decimal places stored in \code{x}.
#' @param print_tex Optional logical value overriding the output format
#'   stored in \code{x}.
#' @param ... Additional arguments, currently unused.
#'
#' @details
#' For console output, columns are sized independently and the table is
#' automatically divided into blocks whenever the complete table would
#' exceed the current console width.
#'
#' The console width can be modified with, for example,
#' \code{options(width = 120)}.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#' data("wqi", package = "unitregTMB")
#'
#' fit <- unitregTMB(
#'   y ~ location,
#'   data = wqi,
#'   family = unitgamma(model_for = "mean")
#' )
#'
#' x <- summary_coef(
#'   fit,
#'   component = c("mu", "phi")
#' )
#'
#' print(x, digits = 2)
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @method print summary_coef.unitregTMB
#' @rdname print_summary_coef_unitregTMB
#' @export
print.summary_coef.unitregTMB <- function(x, digits = NULL, print_tex = NULL, ...) {

  if (!inherits(x, "summary_coef.unitregTMB")) {
    stop("'x' must inherit from 'summary_coef.unitregTMB'.", call. = FALSE)
  }

  if (is.null(digits)) {
    digits <- attr(x, "digits")
  }

  if (is.null(digits)) {
    digits <- 3L
  }

  if (!is.numeric(digits) ||
      length(digits) != 1L ||
      !is.finite(digits) ||
      digits < 0 ||
      digits != as.integer(digits)) {
    stop(
      "'digits' must be a non-negative integer.",
      call. = FALSE
    )
  }

  if (is.null(print_tex)) {
    print_tex <- isTRUE(attr(x, "print_tex"))
  }

  if (length(print_tex) != 1L ||
      !is.logical(print_tex) ||
      is.na(print_tex)) {
    stop(
      "'print_tex' must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  digits <- as.integer(digits)

  model_names <- make.unique(
    vapply(
      x,
      function(z) {
        if (!is.null(z$model_name) &&
            length(z$model_name)) {
          as.character(z$model_name)[1L]
        } else {
          "unitregTMB"
        }
      },
      character(1)
    ),
    sep = " #"
  )

  coef_names <- Reduce(
    union,
    lapply(
      x,
      function(z) {
        if (is.null(z$coef_names)) {
          character()
        } else {
          z$coef_names
        }
      }
    )
  )

  if (!length(coef_names)) {
    cat("No fixed-effect coefficients to display.\n")
    return(invisible(x))
  }

  tab <- matrix(
    "---",
    nrow = length(coef_names),
    ncol = length(x),
    dimnames = list(
      coef_names,
      model_names
    )
  )

  fmt <- paste0("%.", digits, "f")

  for (j in seq_along(x)) {
    z <- x[[j]]

    if (!length(z$coefficients)) next

    for (nm in names(z$coefficients)) {
      est <- unname(z$coefficients[nm])
      se <- unname(z$std_errors[nm])

      tab[nm, j] <- if (length(se) &&
                        is.finite(se)) {
        paste0(sprintf(fmt, est), " (", sprintf(fmt, se), ")")
      } else {
        paste0(sprintf(fmt, est), " (---)")
      }
    }
  }

  display_names <- sub(
    "^mu\\.",
    "",
    rownames(tab)
  )

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

    for (i in seq_len(nrow(tab))) {
      cat(escape_tex(display_names[i]), " & ", paste(tab[i, ], collapse = " & "), " \\\\\n", sep = "")
    }

    cat("\\hline\n")
    cat("\\end{tabular}\n")
    cat("\\end{table}\n")

  } else {

    cat("Regression coefficients (standard errors) ", "for unitregTMB models:\n", sep = "")

    .unitreg_print_blocks(tab = tab, row_labels = display_names, first_header = "Predictor", min_col_width = 15L)
  }

  invisible(x)
}
