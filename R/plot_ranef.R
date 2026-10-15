#' Caterpillar Plot of Ordinary Random Effects
#'
#' @description
#' Produces a caterpillar plot of conditional modes for one ordinary
#' random-effects term in any active \code{unitregTMB} distributional
#' predictor.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param component Component containing the random effect: \code{"mu"},
#'   \code{"phi"}, \code{"p0"}, or \code{"p1"}.
#' @param term Random-effects term to plot. Supply either a positive integer,
#'   an unambiguous grouping-factor name, or the term label returned by
#'   \code{random_effects()}.
#' @param effect Random coefficient to plot, such as \code{"(Intercept)"} or a
#'   slope name. If \code{NULL}, the intercept is used when available;
#'   otherwise the first random coefficient is selected.
#' @param level Numeric scalar in \eqn{(0,1)} giving the confidence level for
#'   intervals based on conditional-mode standard errors.
#' @param ... Additional graphical arguments passed to
#'   \code{\link[graphics]{plot.default}}.
#'
#' @details
#' The plotted intervals are Wald-type intervals based on the conditional-mode
#' standard errors returned by TMB. They summarize uncertainty in the
#' conditional modes and should not be interpreted as confidence intervals for
#' population-level fixed effects.
#'
#' If the same grouping factor occurs in more than one random-effects term, its
#' grouping-factor name is ambiguous; in that case use an integer term index or
#' the explicit term label.
#'
#' @return Invisibly returns a data frame containing the grouping-factor level,
#'   conditional mode, standard error, and lower and upper interval limits used
#'   in the plot.
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
#' plot_ranef(fit, component = "mu")
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom regcore plot_ranef
#' @importFrom stats qnorm
#' @importFrom graphics plot abline segments points grid axis
#' @rdname plot_ranef_unitregTMB
#' @export
plot_ranef.unitregTMB <- function(object, component = c("mu", "phi", "p0", "p1"),
                                  term = 1L, effect = NULL, level = 0.95, ...) {

  if (!inherits(object, "unitregTMB")) {
    stop("'object' must inherit from 'unitregTMB'.", call. = FALSE)
  }

  component <- match.arg(component)

  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("'level' must be a single number in (0, 1).", call. = FALSE)
  }

  if (!.unitreg_active_component(object, component)) {
    stop("Component '", component, "' is not active in this model.", call. = FALSE)
  }

  if (!.unitreg_has_re(object, component)) {
    stop(
      "Component '", component,
      "' does not contain ordinary random effects.",
      call. = FALSE
    )
  }

  dat <- .unitreg_random_effects_long(object, component)

  if (is.null(dat) || !nrow(dat)) {
    stop("Random effects could not be extracted.", call. = FALSE)
  }

  term_map <- unique(
    dat[c("Term", "TermLabel", "Group")]
  )
  term_map <- term_map[order(term_map$Term), , drop = FALSE]

  if (is.character(term)) {
    if (length(term) != 1L || is.na(term) || !nzchar(term)) {
      stop("'term' must be one non-empty character string.", call. = FALSE)
    }

    exact_label <- which(term_map$TermLabel == term)

    if (length(exact_label) == 1L) {
      term_id <- term_map$Term[exact_label]
    } else {
      group_hit <- which(term_map$Group == term)

      if (!length(group_hit)) {
        stop(
          "Unknown random-effects term '", term,
          "'. Available term labels are: ",
          paste(term_map$TermLabel, collapse = ", "),
          ".",
          call. = FALSE
        )
      }

      if (length(group_hit) > 1L) {
        stop(
          "Grouping factor '", term,
          "' occurs in more than one random-effects term; use a term index or ",
          "an explicit term label.",
          call. = FALSE
        )
      }

      term_id <- term_map$Term[group_hit]
    }
  } else {
    if (!is.numeric(term) || length(term) != 1L || !is.finite(term) ||
        term < 1 || term != as.integer(term) ||
        !as.integer(term) %in% term_map$Term) {
      stop(
        "'term' must identify one available random-effects term.",
        call. = FALSE
      )
    }

    term_id <- as.integer(term)
  }

  d <- dat[dat$Term == term_id, , drop = FALSE]
  effects <- unique(d$Effect)

  if (is.null(effect)) {
    effect <- if ("(Intercept)" %in% effects) "(Intercept)" else effects[1L]
  }

  if (!is.character(effect) || length(effect) != 1L ||
      is.na(effect) || !nzchar(effect)) {
    stop("'effect' must be one non-empty character string.", call. = FALSE)
  }

  if (!effect %in% effects) {
    stop(
      "Unknown random coefficient '", effect,
      "'. Available effects are: ",
      paste(effects, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  d <- d[d$Effect == effect, , drop = FALSE]

  if (all(!is.finite(d$Std.Error))) {
    stop(
      "Conditional-mode standard errors are unavailable. Refit with ",
      "control = unitregTMB.control(sd.report = TRUE).",
      call. = FALSE
    )
  }

  if (any(!is.finite(d$Std.Error))) {
    warning(
      "Some conditional-mode standard errors are unavailable; their intervals ",
      "will not be displayed.",
      call. = FALSE
    )
  }

  z <- stats::qnorm(1 - (1 - level) / 2)
  d$Lower <- d$Estimate - z * d$Std.Error
  d$Upper <- d$Estimate + z * d$Std.Error

  ord <- order(d$Estimate)
  d <- d[ord, , drop = FALSE]
  rownames(d) <- NULL

  ypos <- seq_len(nrow(d))
  xr <- range(c(d$Lower, d$Upper, d$Estimate), finite = TRUE)

  main_txt <- paste0(
    "Random effects - ", component,
    " - ", unique(d$Group),
    " - ", effect
  )

  graphics::plot(
    d$Estimate,
    ypos,
    xlim = xr,
    type = "n",
    xlab = "Conditional mode",
    ylab = "Grouping-factor level",
    main = main_txt,
    yaxt = "n",
    bty = "n",
    ...
  )

  graphics::grid(nx = NULL, ny = NA)
  graphics::abline(v = 0, lty = 2)

  ok_interval <- is.finite(d$Lower) & is.finite(d$Upper)
  if (any(ok_interval)) {
    graphics::segments(
      x0 = d$Lower[ok_interval],
      y0 = ypos[ok_interval],
      x1 = d$Upper[ok_interval],
      y1 = ypos[ok_interval]
    )
  }

  graphics::points(d$Estimate, ypos, pch = 16)

  if (nrow(d) <= 40L) {
    graphics::axis(
      side = 2,
      at = ypos,
      labels = d$Level,
      las = 2,
      cex.axis = 0.7
    )
  }

  invisible(d)
}
