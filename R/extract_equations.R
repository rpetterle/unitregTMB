#' Extract Model Equations in LaTeX Format for unitregTMB
#'
#' @description
#' Generates LaTeX equations for fitted \code{unitregTMB} models in symbolic or
#' numeric form. Equations can be extracted for the location, dispersion/shape,
#' zero-inflation, and one-inflation predictors, including ordinary random
#' effects and proper CAR spatial effects when present.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param mode Character string specifying \code{"symbolic"} or
#'   \code{"numeric"} output.
#' @param component Character string specifying \code{"mu"}, \code{"phi"},
#'   \code{"p0"}, \code{"p1"}, or \code{"all"}.
#' @param link_style Character string. \code{"name"} prints the link name;
#'   \code{"inverse"} prints its mathematical expression.
#' @param grouping_var Optional character string identifying the grouping
#'   variable used to display compound \eqn{ij} indexing for predictors with
#'   ordinary random effects.
#' @param digits Non-negative integer giving the number of decimal places used
#'   in numeric mode.
#' @param custom_labels Optional named character vector mapping model-matrix
#'   column names to custom LaTeX labels.
#' @param index_labels Logical; if \code{TRUE}, append observation or cluster
#'   indexing information.
#' @param ... Additional arguments, currently unused.
#'
#' @details
#' For models containing both zero and one augmentation, the boundary
#' predictors are represented using the multinomial-logit parameterization
#' employed by the TMB likelihood, with the continuous component as the
#' baseline category. Numeric equations use the fitted fixed-effect
#' coefficients stored in the fitted model object. Random and spatial effects
#' are represented symbolically through their conditional modes.
#'
#' @return A character string containing the LaTeX \code{align} environment,
#'   returned invisibly after being printed to the console.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data.frame(y = plogis(rnorm(40)), x = rnorm(40))
#' fit <- unitregTMB(
#'   y ~ x,
#'   phi.formula = ~ x,
#'   data = dat,
#'   family = unitgamma(model_for = "mean"),
#' )
#' regcore::extract_equations(fit, mode = "symbolic", component = "all")
#' regcore::extract_equations(fit, mode = "numeric", component = "mu")
#' }
#'
#' @importFrom regcore extract_equations
#' @rdname extract_equations_unitregTMB
#' @export
extract_equations.unitregTMB <- function(
    object,
    mode = c("symbolic", "numeric"),
    component = c("mu", "phi", "p0", "p1", "all"),
    link_style = c("name", "inverse"),
    grouping_var = NULL,
    digits = 2,
    custom_labels = NULL,
    index_labels = TRUE,
    ...) {

  if (!inherits(object, "unitregTMB")) {
    stop("'object' must inherit from 'unitregTMB'.", call. = FALSE)
  }

  mode <- match.arg(mode)
  component <- match.arg(component)
  link_style <- match.arg(link_style)

  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) || digits < 0 || digits != as.integer(digits)) {
    stop("'digits' must be a single non-negative integer.", call. = FALSE)
  }
  digits <- as.integer(digits)

  if (!is.logical(index_labels) || length(index_labels) != 1L || is.na(index_labels)) {
    stop("'index_labels' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(grouping_var) && (!is.character(grouping_var) || length(grouping_var) != 1L || is.na(grouping_var) || !nzchar(grouping_var))) {
    stop("'grouping_var' must be NULL or a single non-empty character string.", call. = FALSE)
  }

  if (!is.null(custom_labels) && (!is.character(custom_labels) || is.null(names(custom_labels)) || any(!nzchar(names(custom_labels))))) {
    stop("'custom_labels' must be a named character vector.", call. = FALSE)
  }

  components <- c("mu", "phi", "p0", "p1")
  dat <- tryCatch(object$obj$env$data, error = function(e) NULL)
  if (is.null(dat)) dat <- list()

  `%||%` <- function(x, y) {
    if (!is.null(x) && length(x) > 0L && !all(is.na(x))) x else y
  }

  object_flag <- function(name) {
    x <- object[[name]]
    !is.null(x) && length(x) > 0L && isTRUE(as.logical(x[1L]))
  }

  data_flag <- function(name) {
    x <- dat[[name]]
    !is.null(x) && length(x) > 0L && isTRUE(as.logical(x[1L]))
  }

  has_p0 <- object_flag("has_p0") || data_flag("has_p0_inflation")
  has_p1 <- object_flag("has_p1") || data_flag("has_p1_inflation")

  active_component <- function(comp) {
    switch(comp, mu = TRUE, phi = TRUE, p0 = has_p0, p1 = has_p1, FALSE)
  }

  has_re <- function(comp) {
    object_flag(paste0("has_random_effects_", comp))
  }

  has_spatial <- function(comp) {
    object_flag(paste0("has_spatial_effects_", comp)) || data_flag(paste0("has_spatial_effects_", comp))
  }

  formula_for <- function(comp) {
    switch(comp, mu = object$formula, phi = object$phi.formula, p0 = object$p0.formula, p1 = object$p1.formula)
  }

  family_object <- if (is.list(object$family_object)) object$family_object else list()

  mu_link <- NULL
  if (!is.null(object$link.mu) && is.list(object$link.mu)) mu_link <- object$link.mu$name
  mu_link <- mu_link %||% family_object$link_r_name %||% "logit"
  mu_link <- as.character(mu_link[1L])

  phi_link <- family_object$phi_link_r_name %||% "log"
  if (is.list(phi_link)) phi_link <- phi_link$name
  phi_link <- as.character(phi_link[1L])

  link_for <- function(comp) {
    switch(
      comp,
      mu = mu_link,
      phi = phi_link,
      p0 = if (has_p0 && has_p1) "bivariate_logit" else "logit",
      p1 = if (has_p0 && has_p1) "bivariate_logit" else "logit"
    )
  }

  component_index <- function(comp) {
    if (has_re(comp) && !is.null(grouping_var)) "ij" else "i"
  }

  parameter_symbol <- function(comp) {
    idx <- component_index(comp)
    switch(
      comp,
      mu = sprintf("\\mu_{%s}", idx),
      phi = sprintf("\\phi_{%s}", idx),
      p0 = sprintf("p_{0,%s}", idx),
      p1 = sprintf("p_{1,%s}", idx)
    )
  }

  get_link_latex <- function(comp) {
    lnk <- link_for(comp)
    par_name <- parameter_symbol(comp)

    if (lnk == "bivariate_logit") {
      idx <- component_index(comp)
      denom <- sprintf("1-p_{0,%s}-p_{1,%s}", idx, idx)
      num <- if (comp == "p0") sprintf("p_{0,%s}", idx) else sprintf("p_{1,%s}", idx)
      return(sprintf("\\log\\left(\\frac{%s}{%s}\\right)", num, denom))
    }

    if (link_style == "name") {
      if (lnk == "log") return(sprintf("\\log(%s)", par_name))
      return(sprintf("\\mathrm{%s}(%s)", lnk, par_name))
    }

    switch(
      lnk,
      logit = sprintf("\\log\\left(\\frac{%s}{1-%s}\\right)", par_name, par_name),
      probit = sprintf("\\Phi^{-1}(%s)", par_name),
      cloglog = sprintf("\\log[-\\log(1-%s)]", par_name),
      cauchit = sprintf("\\tan\\{\\pi(%s-0.5)\\}", par_name),
      log = sprintf("\\log(%s)", par_name),
      sprintf("\\mathrm{%s}(%s)", lnk, par_name)
    )
  }

  coef_symbol <- c(mu = "\\beta", phi = "\\psi", p0 = "\\gamma", p1 = "\\delta")

  format_var_math <- function(v) {
    if (!is.null(custom_labels) && v %in% names(custom_labels)) return(paste0("{", custom_labels[[v]], "}"))
    parts <- strsplit(v, "_", fixed = TRUE)[[1L]]
    if (length(parts) == 1L) return(sprintf("\\mathrm{%s}", parts[1L]))
    out <- sprintf("\\mathrm{%s}", parts[1L])
    if (length(parts) > 1L) {
      for (j in 2:length(parts)) out <- paste0(out, "_{\\mathrm{", parts[j], "}}")
    }
    paste0("{", out, "}")
  }

  variable_index <- function(comp) component_index(comp)

  chunk_symbolic <- function(terms, chunk_size = 3L) {
    if (!length(terms)) return("0")
    chunks <- split(terms, ceiling(seq_along(terms) / chunk_size))
    out <- paste(chunks[[1L]], collapse = " + ")
    if (length(chunks) > 1L) {
      for (j in 2:length(chunks)) {
        out <- paste0(out, " + \\nonumber \\\\ \n  & \\quad + ", paste(chunks[[j]], collapse = " + "))
      }
    }
    out
  }

  chunk_numeric <- function(terms, chunk_size = 3L) {
    if (!length(terms)) return("0")
    chunks <- split(terms, ceiling(seq_along(terms) / chunk_size))
    out <- paste(chunks[[1L]], collapse = " ")
    if (length(chunks) > 1L) {
      for (j in 2:length(chunks)) {
        out <- paste0(out, " \\nonumber \\\\ \n  & \\quad ", paste(chunks[[j]], collapse = " "))
      }
    }
    out
  }

  fixed_columns <- function(comp) {
    b <- object$model_coef[[comp]]
    if (!is.null(b) && length(b) && !is.null(names(b))) return(names(b))
    X <- dat[[paste0("X_", comp)]]
    if (!is.null(X) && !is.null(colnames(X))) return(colnames(X))
    form <- formula_for(comp)
    if (is.null(form)) return(character())
    fixed_form <- reformulas::nobars(form)
    colnames(stats::model.matrix(fixed_form, object$model))
  }

  symbolic_rhs <- function(comp) {
    cols <- fixed_columns(comp)
    sym <- coef_symbol[[comp]]
    rhs <- character()

    if (length(cols)) {
      for (k in seq_along(cols)) {
        term <- cols[k]
        if (identical(term, "(Intercept)")) {
          rhs <- c(rhs, sprintf("%s_{%d}", sym, k - 1L))
        } else {
          rhs <- c(rhs, sprintf("%s_{%d}~%s_{%s}", sym, k - 1L, format_var_math(term), variable_index(comp)))
        }
      }
    }

    if (has_re(comp)) {
      rhs <- c(rhs, sprintf("\\mathbf{z}_{\\mathrm{%s},%s}^{\\top}\\mathbf{u}_{\\mathrm{%s}}", comp, component_index(comp), comp))
    }

    if (has_spatial(comp)) {
      rhs <- c(rhs, sprintf("s_{\\mathrm{%s},r(%s)}", comp, component_index(comp)))
    }

    rhs
  }

  numeric_rhs <- function(comp) {
    b <- object$model_coef[[comp]]
    if (is.null(b)) stop("Fixed-effect estimates are not available for component '", comp, "'.", call. = FALSE)

    term_names <- names(b)
    if (is.null(term_names) || length(term_names) != length(b)) {
      X <- dat[[paste0("X_", comp)]]
      if (!is.null(X) && !is.null(colnames(X)) && ncol(X) == length(b)) {
        term_names <- colnames(X)
      } else {
        term_names <- paste0("x", seq_along(b) - 1L)
      }
    }

    rhs <- character()
    for (k in seq_along(b)) {
      value <- unname(b[k])
      term <- term_names[k]
      if (identical(term, "(Intercept)")) {
        piece <- sprintf(paste0("%.", digits, "f"), value)
      } else {
        abs_value <- sprintf(paste0("%.", digits, "f"), abs(value))
        sign <- if (value < 0) "-" else "+"
        piece <- sprintf("%s %s~%s_{%s}", sign, abs_value, format_var_math(term), variable_index(comp))
      }
      rhs <- c(rhs, piece)
    }

    if (length(rhs) && !identical(term_names[1L], "(Intercept)")) rhs[1L] <- sub("^\\+ ", "", rhs[1L])

    if (has_re(comp)) {
      rhs <- c(rhs, sprintf("+ \\mathbf{z}_{\\mathrm{%s},%s}^{\\top}\\widehat{\\mathbf{u}}_{\\mathrm{%s}}", comp, component_index(comp), comp))
    }

    if (has_spatial(comp)) {
      rhs <- c(rhs, sprintf("+ \\widehat{s}_{\\mathrm{%s},r(%s)}", comp, component_index(comp)))
    }

    rhs
  }

  comps_to_print <- if (component == "all") components else component
  comps_to_print <- comps_to_print[vapply(comps_to_print, active_component, logical(1))]
  if (!length(comps_to_print)) stop("The requested component is not active in this model.", call. = FALSE)

  lines_out <- character()

  for (comp in comps_to_print) {
    lhs <- get_link_latex(comp)
    rhs <- if (mode == "symbolic") chunk_symbolic(symbolic_rhs(comp)) else chunk_numeric(numeric_rhs(comp))
    lines_out <- c(lines_out, sprintf("  %s &= %s \\nonumber", lhs, rhs))
  }

  for (comp in comps_to_print) {
    if (has_re(comp)) {
      lines_out <- c(lines_out, sprintf("  & \\quad \\mathbf{u}_{\\mathrm{%s}} \\sim N(\\mathbf{0}, \\mathbf{\\Sigma}_{\\mathrm{%s}}) \\nonumber", comp, comp))
    }
  }

  if (index_labels) {
    re_components <- comps_to_print[vapply(comps_to_print, has_re, logical(1))]
    if (length(re_components) && !is.null(grouping_var)) {
      first_comp <- re_components[1L]
      re_info <- object[[paste0("re_info_", first_comp)]]
      n_clusters <- NULL
      if (!is.null(re_info) && !is.null(re_info$n_re_levels_list) && length(re_info$n_re_levels_list)) n_clusters <- re_info$n_re_levels_list[1L]
      if (!is.null(n_clusters) && is.finite(n_clusters)) {
        lbl <- sprintf("  & \\quad \\text{for } i = 1, \\dots, %s, \\quad j = 1, \\dots, n_i \\nonumber", n_clusters)
      } else {
        lbl <- "  & \\quad \\text{for clusters } i \\text{ and within-cluster observations } j \\nonumber"
      }
    } else {
      n_obs <- object$nobs %||% if (!is.null(dat$X_mu)) nrow(dat$X_mu) else NA_integer_
      if (is.finite(n_obs)) {
        lbl <- sprintf("  & \\quad \\text{for } i = 1, \\dots, %s \\nonumber", n_obs)
      } else {
        lbl <- "  & \\quad \\text{for observations } i \\nonumber"
      }
    }
    lines_out <- c(lines_out, lbl)
  }

  final_tex <- paste0("\\begin{align}\n", paste(lines_out, collapse = " \\\\\n"), "\n\\end{align}\n")
  cat(final_tex)
  invisible(final_tex)
}
