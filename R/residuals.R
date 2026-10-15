#' Residuals for unitregTMB Models
#'
#' @description
#' Computes randomized quantile, Cox-Snell, or response-scale residuals for
#' fitted \code{unitregTMB} models. Zero/one augmentation, ordinary random
#' effects, and proper CAR spatial effects are handled through the current
#' prediction layer.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param type Residual type: \code{"quantile"}, \code{"cox-snell"}, or
#'   \code{"response"}.
#' @param re.form \code{NULL} computes conditional residuals using fitted latent
#'   effects. \code{NA} or \code{~0} computes population-level residuals.
#' @param seed Optional finite numeric seed used for Dunn-Smyth randomization at
#'   exact zero and one observations.
#' @param ... Additional arguments passed to \code{predict.unitregTMB}.
#'
#' @details
#' Randomized quantile residuals use the probability integral transform of the
#' complete augmented distribution. For an interior observation,
#' \deqn{
#' F_Y(y)=p_0+(1-p_0-p_1)F_c(y),
#' }
#' where \eqn{F_c} is the CDF of the continuous component. Exact zeros are
#' randomized uniformly over \eqn{[0,p_0]}, and exact ones over
#' \eqn{[1-p_1,1]}.
#'
#' The continuous CDF is evaluated only for observations strictly inside
#' \eqn{(0,1)}. This avoids evaluating continuous-family transformations at
#' structural boundary observations.
#'
#' For mean-parameterized families, \code{type = "response"} is the observed
#' response minus the fitted augmented conditional mean. For modal and quantile
#' parameterizations it is the observed response minus the modeled
#' response-scale mode or quantile, respectively, and is therefore not a
#' mean residual.
#'
#' Machine-level guards are applied only when transforming CDF values to
#' Gaussian or Cox-Snell residuals, preventing infinite transformed residuals
#' without altering the observed response.
#'
#' @return A numeric vector of residuals with the same length as the response
#'   used to fit the model.
#'
#' @examples
#' \donttest{
#' data("wqi", package = "unitregTMB")
#'
#' fit <- unitregTMB(
#'   y ~ location + quarter,
#'   data = wqi,
#'   family = unitgamma(model_for = "mean")
#' )
#'
#' head(residuals(fit, type = "quantile", seed = 123))
#' head(residuals(fit, type = "cox-snell", seed = 123))
#' head(residuals(fit, type = "response"))
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom stats predict qnorm runif
#' @export
residuals.unitregTMB <- function(object, type = c("quantile", "cox-snell", "response"),  re.form = NULL, seed = NULL, ...) {

  if (!inherits(object, "unitregTMB")) {
    stop("'object' must inherit from 'unitregTMB'.", call. = FALSE)
  }

  type <- match.arg(type)
  Y <- as.numeric(object$Y)

  if (!length(Y) || any(!is.finite(Y)) || any(Y < 0 | Y > 1)) {
    stop("The stored response must contain finite values in [0, 1].", call. = FALSE)
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("'seed' must be NULL or one finite number.", call. = FALSE)
    }
    set.seed(seed)
  }

  if (type == "response") {
    pred <- stats::predict(
      object,
      type = "response",
      re.form = re.form,
      ...
    )

    if (length(pred) != length(Y)) {
      stop("Prediction length does not match the response length.", call. = FALSE)
    }

    return(Y - as.numeric(pred))
  }

  pred <- stats::predict(object, type = "parameter", re.form = re.form, ...)

  if (!is.data.frame(pred) ||
      !all(c("mu", "phi", "p0", "p1") %in% names(pred)) ||
      nrow(pred) != length(Y)) {
    stop("Distributional-parameter predictions are malformed.", call. = FALSE)
  }

  idx_cont <- Y > 0 & Y < 1
  cdf_cont <- rep(NA_real_, length(Y))

  if (any(idx_cont)) {
    cdf_cont[idx_cont] <- get_cdf_continuous(
      family_code = object$family_object$family_code,
      Y = Y[idx_cont],
      mu = pred$mu[idx_cont],
      phi = pred$phi[idx_cont],
      tau = object$tau
    )

    cdf_cont[idx_cont] <- .unitreg_validate_cdf(
      cdf_cont[idx_cont],
      context = "continuous CDF"
    )
  }

  if (isTRUE(object$has_p0) || isTRUE(object$has_p1)) {
    cdf <- correct_for_inflation(
      Y = Y,
      cdf_cont = cdf_cont,
      p0 = pred$p0,
      p1 = pred$p1
    )
  } else {
    if (any(!idx_cont)) {
      stop(
        "Boundary observations require the corresponding zero/one augmentation.",
        call. = FALSE
      )
    }
    cdf <- cdf_cont
  }

  compute_residuals(cdf, type)
}

#' @keywords internal
#' @noRd
get_cdf_continuous <- function(family_code, Y, mu, phi, tau) {
  code <- as.integer(family_code)[1L]

  if (is.na(code) || !code %in% 0:15) {
    stop("Unknown family code: ", family_code, ".", call. = FALSE)
  }

  switch(
    as.character(code),
    "0" = pbeta_mean(Y, mu, phi),
    "1" = psimplex(Y, mu, phi),
    "2" = pvasicek_mean(Y, mu, phi),
    "3" = pugamma_mean(Y, mu, phi),
    "4" = pbessel(Y, mu, phi),
    "5" = pbeta_mode(Y, mu, phi),
    "6" = pkum_mode(Y, mu, phi),
    "7" = pugamma_mode(Y, mu, phi),
    "8" = pugompertz_mode(Y, mu, phi),
    "9" = pkum_quantile(Y, mu, phi, tau),
    "10" = pvasicek_quantile(Y, mu, phi, tau),
    "11" = puweibull(Y, mu, phi, tau),
    "12" = pugompertz_quantile(Y, mu, phi, tau),
    "13" = pjohnsonsb(Y, mu, phi, tau),
    "14" = pashw(Y, mu, phi, tau),
    "15" = pubs(Y, mu, phi, tau)
  )
}

#' @keywords internal
#' @noRd
get_random_continuous <- function(family_code, n, mu, phi, tau) {
  code <- as.integer(family_code)[1L]

  if (is.na(code) || !code %in% 0:15) {
    stop("Unknown family code: ", family_code, ".", call. = FALSE)
  }

  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 0 || n != as.integer(n)) {
    stop("'n' must be a non-negative integer.", call. = FALSE)
  }

  n <- as.integer(n)

  switch(
    as.character(code),
    "0" = rbeta_mean(n, mu, phi),
    "1" = rsimplex(n, mu, phi),
    "2" = rvasicek_mean(n, mu, phi),
    "3" = rugamma_mean(n, mu, phi),
    "4" = rbessel(n, mu, phi),
    "5" = rbeta_mode(n, mu, phi),
    "6" = rkum_mode(n, mu, phi),
    "7" = rugamma_mode(n, mu, phi),
    "8" = rugompertz_mode(n, mu, phi),
    "9" = rkum_quantile(n, mu, phi, tau),
    "10" = rvasicek_quantile(n, mu, phi, tau),
    "11" = ruweibull(n, mu, phi, tau),
    "12" = rugompertz_quantile(n, mu, phi, tau),
    "13" = rjohnsonsb(n, mu, phi, tau),
    "14" = rashw(n, mu, phi, tau),
    "15" = rubs(n, mu, phi, tau)
  )
}

#' @keywords internal
#' @noRd
.unitreg_validate_cdf <- function(x, context = "CDF") {
  x <- as.numeric(x)

  if (any(!is.finite(x))) {
    stop(context, " contains non-finite values.", call. = FALSE)
  }

  tol <- 100 * .Machine$double.eps

  if (any(x < -tol | x > 1 + tol)) {
    stop(context, " contains values outside [0, 1].", call. = FALSE)
  }

  x[x < 0] <- 0
  x[x > 1] <- 1
  x
}

#' @keywords internal
#' @noRd
correct_for_inflation <- function(Y, cdf_cont, p0, p1) {
  n <- length(Y)

  if (length(cdf_cont) != n) {
    stop("Continuous CDF values have an incompatible length.", call. = FALSE)
  }

  if (length(p0) == 1L) p0 <- rep(p0, n)
  if (length(p1) == 1L) p1 <- rep(p1, n)

  if (length(p0) != n || length(p1) != n) {
    stop("Inflation probabilities have incompatible lengths.", call. = FALSE)
  }

  p0 <- as.numeric(p0)
  p1 <- as.numeric(p1)

  if (any(!is.finite(p0)) || any(!is.finite(p1))) {
    stop("Inflation probabilities must be finite.", call. = FALSE)
  }

  tol <- 100 * .Machine$double.eps

  if (any(p0 < -tol | p0 > 1 + tol) ||
      any(p1 < -tol | p1 > 1 + tol)) {
    stop("Inflation probabilities must lie in [0, 1].", call. = FALSE)
  }

  p0[p0 < 0] <- 0
  p0[p0 > 1] <- 1
  p1[p1 < 0] <- 0
  p1[p1 > 1] <- 1

  sum_p <- p0 + p1

  if (any(sum_p > 1 + tol)) {
    stop("The sum of zero- and one-inflation probabilities exceeds one.", call. = FALSE)
  }

  pc <- 1 - sum_p
  pc[pc < 0 & pc >= -tol] <- 0

  out <- numeric(n)
  idx0 <- Y == 0
  idx1 <- Y == 1
  idxc <- Y > 0 & Y < 1

  if (any(idxc)) {
    if (any(!is.finite(cdf_cont[idxc]))) {
      stop("Continuous CDF values are unavailable for interior observations.", call. = FALSE)
    }

    out[idxc] <- p0[idxc] + pc[idxc] * cdf_cont[idxc]
  }

  if (any(idx0)) {
    out[idx0] <- stats::runif(
      sum(idx0),
      min = 0,
      max = p0[idx0]
    )
  }

  if (any(idx1)) {
    out[idx1] <- stats::runif(
      sum(idx1),
      min = 1 - p1[idx1],
      max = 1
    )
  }

  .unitreg_validate_cdf(out, context = "augmented CDF")
}

#' @keywords internal
#' @noRd
compute_residuals <- function(cdf_val, type) {
  cdf_val <- .unitreg_validate_cdf(cdf_val, context = "CDF")
  eps <- .Machine$double.eps
  cdf_safe <- pmin(pmax(cdf_val, eps), 1 - eps)

  switch(
    type,
    quantile = stats::qnorm(cdf_safe),
    `cox-snell` = -log1p(-cdf_safe),
    stop("Unknown residual type.", call. = FALSE)
  )
}
