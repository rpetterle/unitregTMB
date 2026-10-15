#' Extract Ordinary Random Effects from a unitregTMB Model
#'
#' @description
#' Extracts conditional modes of ordinary Gaussian random effects from one
#' distributional predictor or from all active predictors in a fitted
#' \code{unitregTMB} model.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param component Component to extract: \code{"mu"}, \code{"phi"},
#'   \code{"p0"}, \code{"p1"}, or \code{"all"}. The default is \code{"mu"}.
#' @param ... Additional arguments, currently unused.
#'
#' @details
#' Random effects are returned term by term using the same grouping-level and
#' random-coefficient ordering used to construct the TMB random-effects design
#' matrix. For \code{component = "all"}, inactive components or components
#' without ordinary random effects are returned as \code{NULL}.
#'
#' @return For a single component, a named list with one data frame per
#'   random-effects term. Rows correspond to grouping-factor levels and columns
#'   to random coefficients. For \code{component = "all"}, a named list with
#'   entries \code{mu}, \code{phi}, \code{p0}, and \code{p1} is returned.
#'
#' @examples
#' \donttest{
#' data("wqi", package = "unitregTMB")
#'
#' fit <- unitregTMB(
#'   y ~ location + quarter + (1 | id),
#'   data = wqi,
#'   family = unitgamma(model_for = "mean")
#' )
#'
#' re <- random_effects(fit, component = "mu")
#' head(re[[1]])
#' random_effects(fit, component = "all")
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom regcore random_effects
#' @rdname random_effects_unitregTMB
#' @export
random_effects.unitregTMB <- function(object, component = c("mu", "phi", "p0", "p1", "all"), ...) {

  if (!inherits(object, "unitregTMB")) {
    stop("'object' must inherit from 'unitregTMB'.", call. = FALSE)
  }

  component <- match.arg(component)

  extract_one <- function(comp) {
    if (!.unitreg_active_component(object, comp) ||
        !.unitreg_has_re(object, comp)) {
      return(NULL)
    }

    info <- .unitreg_re_info(object, comp)

    if (is.null(info) || !isTRUE(info$has_random_effects)) {
      stop(
        "Random-effect metadata are missing or malformed for component '",
        comp, "'.",
        call. = FALSE
      )
    }

    sizes <- .unitreg_re_term_sizes(info)
    n_terms <- as.integer(info$n_re_terms)

    if (length(n_terms) != 1L || !is.finite(n_terms) || n_terms < 1L ||
        length(sizes) != n_terms ||
        length(info$n_re_levels_list) != n_terms ||
        length(info$group_var_names) != n_terms ||
        length(info$cnms) != n_terms ||
        length(info$group_levels) != n_terms) {
      stop(
        "Random-effect metadata are inconsistent for component '",
        comp, "'.",
        call. = FALSE
      )
    }

    modes <- NULL

    if (!is.null(object$random_effects_modes)) {
      modes <- object$random_effects_modes[[comp]]
    }

    if (is.null(modes) || !length(modes)) {
      par_full <- .unitreg_full_par(object)
      if (!is.null(par_full)) {
        modes <- par_full[[paste0("re_coef_", comp)]]
      }
    }

    if (is.null(modes)) {
      stop(
        "Conditional modes could not be recovered for component '",
        comp, "'.",
        call. = FALSE
      )
    }

    modes <- as.numeric(modes)
    expected_total <- sum(as.integer(info$n_re_levels_list) * sizes)

    if (length(modes) != expected_total) {
      stop(
        "Dimension mismatch for random effects in component '", comp,
        "': expected ", expected_total, " values but found ", length(modes),
        ".",
        call. = FALSE
      )
    }

    group_names <- as.character(info$group_var_names)
    list_names <- make.unique(group_names)
    out <- vector("list", n_terms)
    names(out) <- list_names
    pos <- 1L

    for (i in seq_len(n_terms)) {
      nlev <- as.integer(info$n_re_levels_list[i])
      ncomp <- as.integer(sizes[i])
      len <- nlev * ncomp
      idx <- pos:(pos + len - 1L)
      vals <- modes[idx]

      mat <- matrix(
        vals,
        nrow = nlev,
        ncol = ncomp,
        byrow = TRUE
      )

      cn <- info$cnms[[i]]
      lev <- info$group_levels[[i]]

      if (is.null(cn) || length(cn) != ncomp) {
        cn <- paste0("effect", seq_len(ncomp))
      }

      if (is.null(lev) || length(lev) != nlev) {
        lev <- as.character(seq_len(nlev))
      }

      colnames(mat) <- as.character(cn)
      rownames(mat) <- as.character(lev)

      out[[i]] <- as.data.frame(
        mat,
        check.names = FALSE
      )

      attr(out[[i]], "component") <- comp
      attr(out[[i]], "group") <- group_names[i]
      attr(out[[i]], "term") <- i

      pos <- pos + len
    }

    class(out) <- c("unitregTMB_random_effects", "list")
    attr(out, "component") <- comp
    out
  }

  if (component == "all") {
    out <- setNames(
      lapply(.unitreg_components, extract_one),
      .unitreg_components
    )
    class(out) <- c("unitregTMB_random_effects_all", "list")
    return(out)
  }

  extract_one(component)
}

#' @keywords internal
#' @noRd
.unitreg_re_term_sizes <- function(info) {
  sizes <- info$re_term_size

  if (is.null(sizes)) {
    sizes <- info$re_term_n_components
  }

  if (is.null(sizes)) {
    return(integer())
  }

  as.integer(sizes)
}

#' @keywords internal
#' @noRd
.unitreg_random_effects_long <- function(object, component) {
  if (!.unitreg_active_component(object, component) ||
      !.unitreg_has_re(object, component)) {
    return(NULL)
  }

  info <- .unitreg_re_info(object, component)
  re <- random_effects(object, component = component)

  if (is.null(re)) return(NULL)

  sizes <- .unitreg_re_term_sizes(info)
  total <- sum(as.integer(info$n_re_levels_list) * sizes)
  se_vec <- rep(NA_real_, total)

  if (!is.null(object$sd_report)) {
    sdr <- tryCatch(
      summary(object$sd_report, "random"),
      error = function(e) NULL
    )

    if (is.matrix(sdr) && nrow(sdr)) {
      rn <- rownames(sdr)
      target <- paste0("re_coef_", component)
      idx <- rn == target |
        grepl(paste0("^", target, "\\["), rn)

      if (sum(idx) == length(se_vec)) {
        se_col <- if ("Std. Error" %in% colnames(sdr)) {
          "Std. Error"
        } else {
          2L
        }
        se_vec <- as.numeric(sdr[idx, se_col])
      }
    }
  }

  group_names <- as.character(info$group_var_names)
  list_names <- names(re)
  ans <- vector("list", info$n_re_terms)
  pos <- 1L

  for (i in seq_len(info$n_re_terms)) {
    nlev <- as.integer(info$n_re_levels_list[i])
    ncomp <- as.integer(sizes[i])
    len <- nlev * ncomp
    est_mat <- as.matrix(re[[i]])

    se_mat <- matrix(
      se_vec[pos:(pos + len - 1L)],
      nrow = nlev,
      ncol = ncomp,
      byrow = TRUE
    )

    cn <- colnames(est_mat)
    lev <- rownames(est_mat)

    ans[[i]] <- data.frame(
      Component = component,
      Term = i,
      Group = group_names[i],
      TermLabel = list_names[i],
      Level = rep(lev, each = ncomp),
      Effect = rep(cn, times = nlev),
      Estimate = as.vector(t(est_mat)),
      Std.Error = as.vector(t(se_mat)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    pos <- pos + len
  }

  do.call(rbind, ans)
}

#' Print Extracted Random Effects
#'
#' @param x An object of class \code{unitregTMB_random_effects}.
#' @param digits Number of significant digits used when printing estimates.
#' @param header Logical; if \code{TRUE}, prints the main heading.
#' @param ... Additional arguments passed to \code{print.data.frame}.
#'
#' @return Invisibly returns \code{x}.
#' @export
#' @noRd
print.unitregTMB_random_effects <- function(x, digits = max(3L, getOption("digits") - 2L), header = TRUE, ...) {

  component <- attr(x, "component")
  if (is.null(component) || !length(component)) component <- "unknown"

  component_label <- switch(
    component,
    mu = "mu",
    phi = "phi",
    p0 = "p0",
    p1 = "p1",
    component
  )

  if (isTRUE(header)) cat("\nRandom effects (conditional modes)\n")
  cat("Component:", component_label, "\n")

  if (!length(x)) {
    cat("No ordinary random effects.\n")
    return(invisible(x))
  }

  for (i in seq_along(x)) {
    tab <- x[[i]]
    group <- attr(tab, "group")
    term <- attr(tab, "term")

    if (is.null(group) || !length(group)) group <- names(x)[i]
    if (is.null(term) || !length(term)) term <- i

    effect_names <- colnames(tab)
    rhs <- effect_names
    rhs[rhs == "(Intercept)"] <- "1"
    term_label <- paste0("(", paste(rhs, collapse = " + "), " | ", group, ")")

    cat("\nTerm ", term, ": ", term_label, "\n", sep = "")
    cat("Groups: ", group, " (", nrow(tab), " levels)\n", sep = "")

    out <- data.frame(
      Level = rownames(tab),
      tab,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    rownames(out) <- NULL

    print.data.frame(
      out,
      row.names = FALSE,
      digits = digits,
      ...
    )
  }

  invisible(x)
}

#' Print Extracted Random Effects from All Components
#'
#' @param x An object of class \code{unitregTMB_random_effects_all}.
#' @param digits Number of significant digits used when printing estimates.
#' @param ... Additional arguments passed to component print methods.
#'
#' @return Invisibly returns \code{x}.
#' @export
#' @noRd
print.unitregTMB_random_effects_all <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {

  cat("\nRandom effects (conditional modes)\n")
  cat("All distributional components\n")

  active <- names(x)[!vapply(x, is.null, logical(1))]

  if (!length(active)) {
    cat("No ordinary random effects.\n")
    return(invisible(x))
  }

  for (comp in active) {
    cat("\n", paste(rep("-", 60L), collapse = ""), "\n", sep = "")
    print.unitregTMB_random_effects(
      x[[comp]],
      digits = digits,
      header = FALSE,
      ...
    )
  }

  invisible(x)
}

#' Convert Extracted Random Effects to a Long Data Frame
#'
#' @param x An object of class \code{unitregTMB_random_effects} or
#'   \code{unitregTMB_random_effects_all}.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Additional arguments, currently unused.
#'
#' @return A data frame with columns \code{Component}, \code{Term},
#'   \code{Group}, \code{Level}, \code{Effect}, and \code{Estimate}.
#' @noRd
as.data.frame.unitregTMB_random_effects <- function(x, row.names = NULL, optional = FALSE, ...) {
  component <- attr(x, "component")
  if (is.null(component) || !length(component)) component <- NA_character_

  if (!length(x)) {
    return(data.frame(
      Component = character(), Term = integer(), Group = character(),
      Level = character(), Effect = character(), Estimate = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  ans <- lapply(seq_along(x), function(i) {
    tab <- x[[i]]
    group <- attr(tab, "group")
    term <- attr(tab, "term")
    if (is.null(group) || !length(group)) group <- names(x)[i]
    if (is.null(term) || !length(term)) term <- i

    nr <- nrow(tab)
    nc <- ncol(tab)
    data.frame(
      Component = rep(component, nr * nc),
      Term = rep(as.integer(term), nr * nc),
      Group = rep(as.character(group), nr * nc),
      Level = rep(rownames(tab), each = nc),
      Effect = rep(colnames(tab), times = nr),
      Estimate = as.vector(t(as.matrix(tab))),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  ans <- do.call(rbind, ans)
  rownames(ans) <- NULL
  ans
}

#' @noRd
as.data.frame.unitregTMB_random_effects_all <- function(x, row.names = NULL, optional = FALSE, ...) {
  ans <- lapply(x, function(z) {
    if (is.null(z)) return(NULL)
    as.data.frame(z)
  })
  ans <- ans[!vapply(ans, is.null, logical(1))]

  if (!length(ans)) {
    return(data.frame(
      Component = character(), Term = integer(), Group = character(),
      Level = character(), Effect = character(), Estimate = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, ans)
  rownames(out) <- NULL
  out
}

#' Summarize Extracted Random Effects
#'
#' @param object An object of class \code{unitregTMB_random_effects} or
#'   \code{unitregTMB_random_effects_all}.
#' @param ... Additional arguments, currently unused.
#'
#' @return A summary data frame with one row per component, term, group, and
#'   random coefficient.
#' @export
#' @noRd
summary.unitregTMB_random_effects <- function(object, ...) {
  d <- as.data.frame(object)
  if (!nrow(d)) return(d)

  key <- interaction(d$Component, d$Term, d$Group, d$Effect, drop = TRUE, lex.order = TRUE)
  spl <- split(d, key)

  ans <- lapply(spl, function(z) {
    x <- z$Estimate
    data.frame(
      Component = z$Component[1L],
      Term = z$Term[1L],
      Group = z$Group[1L],
      Effect = z$Effect[1L],
      Levels = length(x),
      Mean = mean(x),
      SD = if (length(x) > 1L) stats::sd(x) else NA_real_,
      Min = min(x),
      Q1 = unname(stats::quantile(x, 0.25, names = FALSE)),
      Median = stats::median(x),
      Q3 = unname(stats::quantile(x, 0.75, names = FALSE)),
      Max = max(x),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  out <- do.call(rbind, ans)
  rownames(out) <- NULL
  class(out) <- c("summary.unitregTMB_random_effects", "data.frame")
  out
}

#' @noRd
summary.unitregTMB_random_effects_all <- function(object, ...) {
  d <- as.data.frame(object)
  if (!nrow(d)) return(d)

  key <- interaction(d$Component, d$Term, d$Group, d$Effect, drop = TRUE, lex.order = TRUE)
  spl <- split(d, key)

  ans <- lapply(spl, function(z) {
    x <- z$Estimate
    data.frame(
      Component = z$Component[1L],
      Term = z$Term[1L],
      Group = z$Group[1L],
      Effect = z$Effect[1L],
      Levels = length(x),
      Mean = mean(x),
      SD = if (length(x) > 1L) stats::sd(x) else NA_real_,
      Min = min(x),
      Q1 = unname(stats::quantile(x, 0.25, names = FALSE)),
      Median = stats::median(x),
      Q3 = unname(stats::quantile(x, 0.75, names = FALSE)),
      Max = max(x),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  out <- do.call(rbind, ans)
  rownames(out) <- NULL
  class(out) <- c("summary.unitregTMB_random_effects", "data.frame")
  out
}

#' @noRd
print.summary.unitregTMB_random_effects <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("\nSummary of random effects (conditional modes)\n")
  print.data.frame(x, row.names = FALSE, digits = digits, ...)
  invisible(x)
}
