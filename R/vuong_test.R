#' Vuong Test for Non-Nested unitregTMB Models
#'
#' @description
#' Performs the Vuong test for two non-nested \code{unitregTMB} models using
#' observation-level log-likelihood contributions from the response model.
#'
#' @param object1,object2 Fitted \code{unitregTMB} objects.
#' @param alternative Alternative hypothesis: \code{"two.sided"},
#'   \code{"less"}, or \code{"greater"}. The \code{"greater"} alternative tests
#'   whether \code{object1} has the larger expected pointwise log-likelihood.
#' @param allow_conditional Logical. Models containing ordinary random effects
#'   or spatial effects do not have an observation-wise additive marginal
#'   likelihood after latent effects are integrated out. If \code{FALSE}, the
#'   comparison stops. If \code{TRUE}, a conditional diagnostic is computed.
#' @param allow_weighted Logical. Standard Vuong asymptotics do not directly
#'   cover arbitrary non-unit prior weights. If \code{FALSE}, weighted models
#'   are rejected. If \code{TRUE}, a weighted pointwise diagnostic is computed.
#' @param ... Additional arguments, currently unused.
#'
#' @details
#' Let \eqn{m_i = \ell_{1i} - \ell_{2i}} denote the pointwise log-likelihood
#' difference. The statistic is
#' \deqn{
#' V = \frac{\sqrt{n}\,\bar m}{s_m},
#' }
#' where \eqn{s_m^2} is the sample variance of the \eqn{m_i}. Under the usual
#' regularity conditions for independent observations and genuinely non-nested
#' models, the null distribution is asymptotically standard normal.
#'
#' The function does not attempt to determine whether two arbitrary
#' multi-parameter \code{unitregTMB} specifications are nested. The user should
#' use this test only for non-nested comparisons.
#'
#' When latent effects are present, the marginal likelihood contains
#' cluster- or field-level integrals and is not additive over observations.
#' Setting \code{allow_conditional = TRUE} therefore produces only a diagnostic
#' comparison based on the conditional response log-likelihood evaluated at the
#' fitted latent-effect modes.
#'
#' @return An object of class \code{htest}.
#'
#' @references
#' Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
#' non-nested hypotheses. \emph{Econometrica}, \bold{57}(2), 307--333.
#'
#' @examples
#' \donttest{
#' data("wqi", package = "unitregTMB")
#'
#' fit_beta <- unitregTMB(y ~ location + quarter,
#'                        family = beta_fam(model_for = "mean"),
#'                        data = wqi)
#'
#' fit_ug <- unitregTMB(y ~ location + quarter,
#'                      family = unitgamma(model_for = "mean"),
#'                      data = wqi)
#'
#' vuong_test(fit_beta, fit_ug)
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom regcore vuong_test
#' @importFrom stats pnorm var
#' @method vuong_test unitregTMB
#' @rdname vuong_test_unitregTMB
#' @export
vuong_test.unitregTMB <- function(object1, object2, alternative = c("two.sided", "less", "greater"),
                                  allow_conditional = FALSE, allow_weighted = FALSE, ...) {

  if (!inherits(object1, "unitregTMB") ||
      !inherits(object2, "unitregTMB")) {
    stop("Both objects must inherit from 'unitregTMB'.", call. = FALSE)
  }

  alternative <- match.arg(alternative)

  for (x in c("allow_conditional", "allow_weighted")) {
    value <- get(x, inherits = FALSE)

    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      stop("'", x, "' must be TRUE or FALSE.", call. = FALSE)
    }
  }

  if (!identical(as.integer(object1$nobs), as.integer(object2$nobs))) {
    stop("Models must use the same number of observations.", call. = FALSE)
  }

  if (is.null(object1$Y) || is.null(object2$Y) ||
      !isTRUE(all.equal(
        as.numeric(object1$Y),
        as.numeric(object2$Y),
        tolerance = 0,
        check.attributes = FALSE
      ))) {
    stop(
      "Models must be fitted to the same response values and observation order.",
      call. = FALSE
    )
  }

  if (!is.null(object1$model) && !is.null(object2$model)) {
    rn1 <- rownames(object1$model)
    rn2 <- rownames(object2$model)

    if (!is.null(rn1) && !is.null(rn2) && !identical(rn1, rn2)) {
      stop(
        "Models must use the same analysis rows and observation order.",
        call. = FALSE
      )
    }
  }

  w1 <- object1$weights
  w2 <- object2$weights

  if (is.null(w1)) w1 <- rep(1, object1$nobs)
  if (is.null(w2)) w2 <- rep(1, object2$nobs)

  if (!isTRUE(all.equal(
    as.numeric(w1),
    as.numeric(w2),
    tolerance = 0,
    check.attributes = FALSE
  ))) {
    stop(
      "Models must use the same observation weights for a Vuong comparison.",
      call. = FALSE
    )
  }

  use_weighted_ll <- any(as.numeric(w1) != 1)

  if (use_weighted_ll && !allow_weighted) {
    stop(
      "Non-unit observation weights are present. Standard Vuong asymptotics ",
      "do not directly cover arbitrary prior weights. Set allow_weighted = TRUE ",
      "only to request a weighted diagnostic comparison.",
      call. = FALSE
    )
  }

  if (use_weighted_ll && allow_weighted) {
    warning(
      "Using weighted pointwise log-likelihood contributions. Interpret this ",
      "as a diagnostic comparison unless the weighting scheme justifies the ",
      "corresponding asymptotics.",
      call. = FALSE
    )
  }

  has_latent <- function(model) {
    any(vapply(
      .unitreg_components,
      function(comp) {
        .unitreg_has_re(model, comp) ||
          .unitreg_has_spatial(model, comp)
      },
      logical(1)
    ))
  }

  latent <- has_latent(object1) || has_latent(object2)

  if (latent && !allow_conditional) {
    stop(
      "At least one model contains ordinary random effects or spatial effects. ",
      "The marginal likelihood is not observation-wise additive after latent ",
      "effects are integrated out, so the standard Vuong test is not valid. ",
      "Set allow_conditional = TRUE only for a conditional diagnostic.",
      call. = FALSE
    )
  }

  if (latent && allow_conditional) {
    warning(
      "Using conditional response log-likelihoods evaluated at fitted latent-",
      "effect modes. Interpret this as a diagnostic comparison, not as the ",
      "standard marginal Vuong test.",
      call. = FALSE
    )
  }

  get_pointwise_ll <- function(model) {
    rep <- tryCatch(
      model$obj$report(),
      error = function(e) NULL
    )

    if (is.null(rep)) {
      stop(
        "The TMB report could not be evaluated for one of the models.",
        call. = FALSE
      )
    }

    if (use_weighted_ll && !is.null(rep$conditional_log_lik_weighted)) {
      ll <- rep$conditional_log_lik_weighted
    } else if (!is.null(rep$conditional_log_lik_unweighted)) {
      ll <- rep$conditional_log_lik_unweighted
    } else if (!is.null(rep$ll_i)) {
      ll <- rep$ll_i
    } else if (!is.null(rep$nll_i)) {
      ll <- -rep$nll_i
    } else {
      ll <- NULL
    }

    if (is.null(ll)) {
      stop(
        "Pointwise log-likelihood contributions are unavailable. The C++ ",
        "template must REPORT conditional_log_lik_unweighted and, when ",
        "weighted diagnostics are needed, conditional_log_lik_weighted.",
        call. = FALSE
      )
    }

    ll <- as.numeric(ll)

    if (length(ll) != model$nobs || any(!is.finite(ll))) {
      stop(
        "Invalid pointwise log-likelihood contributions were returned.",
        call. = FALSE
      )
    }

    ll
  }

  ll1 <- get_pointwise_ll(object1)
  ll2 <- get_pointwise_ll(object2)
  m <- ll1 - ll2
  mean_m <- mean(m)
  sm2 <- stats::var(m)

  if (!is.finite(sm2) || sm2 <= .Machine$double.eps) {
    warning(
      "The variance of the pointwise log-likelihood difference is effectively ",
      "zero; the models are numerically indistinguishable by this statistic.",
      call. = FALSE
    )
    stat <- 0
    sd_m <- 0
  } else {
    sd_m <- sqrt(sm2)
    stat <- sqrt(length(m)) * mean_m / sd_m
  }

  names(stat) <- "Vuong z"

  pvalue <- switch(
    alternative,
    two.sided = 2 * stats::pnorm(abs(stat), lower.tail = FALSE),
    greater = stats::pnorm(stat, lower.tail = FALSE),
    less = stats::pnorm(stat, lower.tail = TRUE)
  )

  family_label <- function(model) {
    fam <- if (!is.null(model$family) && length(model$family)) {
      as.character(model$family)[1L]
    } else {
      "unitregTMB"
    }

    role <- model$family_object$model_for

    if (!is.null(role) && length(role)) {
      role <- as.character(role)[1L]
      fam <- paste0(fam, " [", role, "]")
    }

    if (identical(role, "quantile") &&
        !is.null(model$tau) &&
        length(model$tau) == 1L &&
        is.finite(model$tau)) {
      fam <- paste0(fam, " (tau = ", format(model$tau, trim = TRUE), ")")
    }

    fam
  }

  fam1 <- family_label(object1)
  fam2 <- family_label(object2)

  method <- if (latent) {
    "Conditional pointwise Vuong diagnostic for unitregTMB models"
  } else if (use_weighted_ll) {
    "Weighted pointwise Vuong diagnostic for unitregTMB models"
  } else {
    "Vuong likelihood-ratio test for non-nested unitregTMB models"
  }

  out <- list(statistic = stat, p.value = unname(pvalue),
              estimate = c(`mean log-likelihood difference` = mean_m,`sd log-likelihood difference` = sd_m),
              null.value = c(`mean log-likelihood difference` = 0),
              alternative = alternative,
              method = method,
              data.name = paste0(fam1, " versus ", fam2))

  class(out) <- "htest"
  out
}

#' Pairwise Vuong Comparisons for unitregTMB Models
#'
#' @description
#' Computes pairwise Vuong comparisons among two or more fitted
#' \code{unitregTMB} models and adjusts the resulting p-values for multiple
#' testing.
#'
#' @param ... Fitted \code{unitregTMB} objects.
#' @param lt Optional list of fitted \code{unitregTMB} objects. If supplied,
#'   objects in \code{...} are ignored.
#' @param p.adjust.method Multiple-testing adjustment method passed to
#'   \code{\link[stats]{p.adjust}} through \code{pairwise.table}. The default is
#'   \code{"holm"}.
#' @param alternative Alternative hypothesis passed to
#'   \code{\link{vuong_test.unitregTMB}}.
#' @param allow_conditional Logical; passed to
#'   \code{\link{vuong_test.unitregTMB}}.
#' @param allow_weighted Logical; passed to
#'   \code{\link{vuong_test.unitregTMB}}.
#'
#' @details
#' The same limitations as for \code{\link{vuong_test.unitregTMB}} apply to
#' every pair. In particular, the standard test is intended for non-nested
#' models fitted to independent observations.
#'
#' @return An object of class \code{pairwise.htest}.
#'
#' @examples
#' \donttest{
#' data("wqi", package = "unitregTMB")
#'
#' fit_beta <- unitregTMB(y ~ location,
#'                        family = beta_fam(model_for = "mean"),
#'                        data = wqi)
#'
#' fit_ug <- unitregTMB(y ~ location,
#'                      family = unitgamma(model_for = "mean"),
#'                      data = wqi)
#'
#' pairwise_vuong_test(fit_beta, fit_ug)
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom regcore pairwise_vuong_test vuong_test
#' @importFrom stats pairwise.table p.adjust.methods
#' @method pairwise_vuong_test unitregTMB
#' @rdname pairwise_vuong_test_unitregTMB
#' @export
pairwise_vuong_test.unitregTMB <- function(lt = NULL, p.adjust.method = "holm", alternative = c("two.sided", "less", "greater"),
                                           allow_conditional = FALSE, allow_weighted = FALSE, ...) {

  if (is.null(lt)) lt <- list(...)

  if (!is.list(lt) || length(lt) < 2L) {
    stop("At least two 'unitregTMB' models are required.", call. = FALSE)
  }

  ok <- vapply(lt, inherits, logical(1), what = "unitregTMB")

  if (!all(ok)) {
    stop("All supplied objects must inherit from 'unitregTMB'.", call. = FALSE)
  }

  for (x in c("allow_conditional", "allow_weighted")) {
    value <- get(x, inherits = FALSE)

    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      stop("'", x, "' must be TRUE or FALSE.", call. = FALSE)
    }
  }

  p.adjust.method <- match.arg(p.adjust.method, choices = stats::p.adjust.methods)

  alternative <- match.arg(alternative)

  family_label <- function(model) {
    fam <- if (!is.null(model$family) && length(model$family)) {
      as.character(model$family)[1L]
    } else {
      "unitregTMB"
    }

    role <- model$family_object$model_for

    if (!is.null(role) && length(role)) {
      role <- as.character(role)[1L]
      fam <- paste0(fam, " [", role, "]")
    }

    if (identical(role, "quantile") &&
        !is.null(model$tau) &&
        length(model$tau) == 1L &&
        is.finite(model$tau)) {
      fam <- paste0(fam, " (tau = ", format(model$tau, trim = TRUE), ")")
    }

    fam
  }

  labels <- make.unique(vapply(lt, family_label, character(1)))
  nmod <- length(lt)

  compare_levels <- function(i, j) {
    vuong_test(
      lt[[as.integer(i)]],
      lt[[as.integer(j)]],
      alternative = alternative,
      allow_conditional = allow_conditional,
      allow_weighted = allow_weighted
    )$p.value
  }

  pval <- stats::pairwise.table(
    compare.levels = compare_levels,
    level.names = seq_len(nmod),
    p.adjust.method = p.adjust.method
  )

  colnames(pval) <- labels[seq_len(nmod - 1L)]
  rownames(pval) <- labels[2L:nmod]

  out <- list(method = "Pairwise Vuong comparisons for unitregTMB models",
              data.name = paste(labels, collapse = ", "),
              p.value = pval,
              p.adjust.method = p.adjust.method,
              alternative = alternative,
              conditional = allow_conditional,
              weighted = allow_weighted)

  class(out) <- "pairwise.htest"
  out
}
