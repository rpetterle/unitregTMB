#' Apply Machine-Level Guards to Open-Unit Values
#'
#' @description
#' Internal utility that validates numeric values intended to lie in the open
#' unit interval and replaces only exact floating-point boundary values by
#' machine-representable interior values.
#'
#' @param x Numeric vector.
#' @param name Character label used in error messages.
#'
#' @return A numeric vector strictly inside \eqn{(0,1)}.
#'
#' @examples
#' unitregTMB:::.unitreg_open_unit_guard(c(0, 0.25, 1), "x")
#'
#' @noRd
.unitreg_open_unit_guard <- function(x, name = "value") {
  x <- as.numeric(x)

  if (any(!is.finite(x))) {
    stop("'", name, "' contains non-finite values.", call. = FALSE)
  }

  if (any(x < 0 | x > 1)) {
    stop("'", name, "' contains values outside [0, 1].", call. = FALSE)
  }

  x[x == 0] <- .Machine$double.xmin
  x[x == 1] <- 1 - .Machine$double.eps
  x
}

#' Validate Zero/One-Augmentation Probabilities
#'
#' @description
#' Internal utility that validates zero-, one-, and continuous-component
#' probabilities without substantive clipping or renormalization.
#'
#' @param p0 Numeric zero-mass probabilities.
#' @param p1 Numeric one-mass probabilities.
#'
#' @return A list containing \code{p0}, \code{p1}, and \code{pc}.
#'
#' @examples
#' unitregTMB:::.unitreg_validate_aug_probs(c(0.1, 0.2), c(0.05, 0.1))
#'
#' @noRd
.unitreg_validate_aug_probs <- function(p0, p1) {
  p0 <- as.numeric(p0)
  p1 <- as.numeric(p1)

  if (length(p0) != length(p1)) {
    stop("'p0' and 'p1' must have the same length.", call. = FALSE)
  }

  if (any(!is.finite(p0)) || any(!is.finite(p1))) {
    stop("Augmentation probabilities must be finite.", call. = FALSE)
  }

  tol <- 100 * .Machine$double.eps

  if (any(p0 < -tol | p0 > 1 + tol) ||
      any(p1 < -tol | p1 > 1 + tol)) {
    stop("Augmentation probabilities must lie in [0, 1].", call. = FALSE)
  }

  p0[p0 < 0] <- 0
  p0[p0 > 1] <- 1
  p1[p1 < 0] <- 0
  p1[p1 > 1] <- 1

  pc <- 1 - p0 - p1

  if (any(pc < -tol)) {
    stop(
      "The sum of zero- and one-augmentation probabilities exceeds one.",
      call. = FALSE
    )
  }

  pc[pc < 0 & pc >= -tol] <- 0
  list(p0 = p0, p1 = p1, pc = pc)
}

#' Simulate Responses from a Fitted unitregTMB Model
#'
#' @description
#' Generates responses from the fitted continuous distribution and optional
#' zero/one augmentation. By default, simulations are conditional on the fitted
#' ordinary random effects and proper CAR spatial effects in all active
#' distributional predictors.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param nsim Positive integer giving the number of response vectors to
#'   simulate.
#' @param seed Optional finite numeric random-number seed.
#' @param newdata Optional data frame at which to simulate.
#' @param re.form \code{NULL} uses fitted ordinary and spatial latent effects.
#'   \code{NA} or \code{~0} excludes them and gives population-level
#'   simulations.
#' @param allow.new.levels Logical; passed to \code{predict.unitregTMB} when
#'   \code{newdata} is supplied.
#' @param offset Optional numeric offset for the \code{mu} predictor in
#'   \code{newdata}.
#' @param ... Additional arguments passed to \code{predict.unitregTMB}.
#'
#' @details
#' This method conditions on fitted latent-effect modes when
#' \code{re.form = NULL}; it does not draw new Gaussian random effects or a new
#' CAR field. Use \code{\link{sim_unitreg}} when a complete data-generating
#' process with newly generated latent effects is required.
#'
#' The continuous response is not clipped at an arbitrary tolerance. Exact
#' floating-point values of zero or one returned by a continuous random
#' generator are moved only to machine-representable interior values because
#' the supported continuous components have open-unit support. Positive values
#' arbitrarily close to either boundary are otherwise left unchanged.
#'
#' @return A data frame with \code{nsim} simulated response vectors.
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
#' sim1 <- simulate(fit, nsim = 2, seed = 123)
#' head(sim1)
#' sim2 <- simulate(fit, nsim = 2, seed = 123, re.form = NA)
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom stats predict runif
#' @method simulate unitregTMB
#' @export
simulate.unitregTMB <- function(object, nsim = 1, seed = NULL, newdata = NULL, re.form = NULL,
                                allow.new.levels = FALSE, offset = NULL, ...) {

  if (!inherits(object, "unitregTMB")) {
    stop("'object' must inherit from 'unitregTMB'.", call. = FALSE)
  }

  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim < 1 || nsim != as.integer(nsim)) {
    stop("'nsim' must be a positive integer.", call. = FALSE)
  }

  nsim <- as.integer(nsim)

  if (length(allow.new.levels) != 1L ||
      !is.logical(allow.new.levels) ||
      is.na(allow.new.levels)) {
    stop("'allow.new.levels' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("'seed' must be NULL or one finite number.", call. = FALSE)
    }
    set.seed(seed)
  }

  pred <- stats::predict(object, newdata = newdata, type = "parameter", re.form = re.form,
                         allow.new.levels = allow.new.levels, offset = offset, ...)

  if (!is.data.frame(pred) ||
      !all(c("mu", "phi", "p0", "p1") %in% names(pred))) {
    stop("Parameter predictions are malformed.", call. = FALSE)
  }

  n <- nrow(pred)

  if (n < 1L) {
    stop("No observations are available for simulation.", call. = FALSE)
  }

  family_code <- object$family_object$family_code
  tau <- object$tau
  mu <- .unitreg_open_unit_guard(pred$mu, "predicted mu")
  phi <- as.numeric(pred$phi)

  if (length(phi) != n || any(!is.finite(phi)) || any(phi <= 0)) {
    stop("Predicted phi values must be positive and finite.", call. = FALSE)
  }

  probs <- .unitreg_validate_aug_probs(pred$p0, pred$p1)
  p0 <- probs$p0
  p1 <- probs$p1

  out <- vector("list", nsim)

  for (s in seq_len(nsim)) {
    u <- stats::runif(n)
    is0 <- u < p0
    is1 <- !is0 & u < (p0 + p1)
    isc <- !(is0 | is1)

    y <- numeric(n)
    y[is0] <- 0
    y[is1] <- 1

    if (any(isc)) {
      y_cont <- get_random_continuous(
        family_code = family_code,
        n = sum(isc),
        mu = mu[isc],
        phi = phi[isc],
        tau = tau
      )

      if (length(y_cont) != sum(isc)) {
        stop(
          "The continuous random generator returned an incompatible length.",
          call. = FALSE
        )
      }

      y[isc] <- .unitreg_open_unit_guard(
        y_cont,
        "simulated continuous response"
      )
    }

    out[[s]] <- y
  }

  ans <- as.data.frame(out)
  names(ans) <- paste0("sim_", seq_len(nsim))
  attr(ans, "seed") <- seed
  attr(ans, "re.form") <- re.form
  ans
}
