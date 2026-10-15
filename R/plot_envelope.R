# ==============================================================================
# Simulated residual envelopes for unitregTMB
# ==============================================================================

.unitreg_env_or <- function(x, y) if (is.null(x)) y else x

.unitreg_env_components <- c("mu", "phi", "p0", "p1")

.unitreg_env_matvec <- function(M, v, label = "matrix product") {
  if (is.null(M)) stop("Missing matrix in ", label, ".", call. = FALSE)

  ans <- tryCatch(as.numeric(M %*% v), error = function(e) e)

  if (inherits(ans, "error")) {
    stop(label, " failed: ", conditionMessage(ans), call. = FALSE)
  }

  ans
}

.unitreg_env_has_re <- function(object, component) {
  isTRUE(object[[paste0("has_random_effects_", component)]])
}

.unitreg_env_has_spatial <- function(object, component) {
  isTRUE(object[[paste0("has_spatial_effects_", component)]])
}

.unitreg_env_has_any_re <- function(object) {
  any(vapply(
    .unitreg_env_components,
    function(component) .unitreg_env_has_re(object, component),
    logical(1)
  ))
}

.unitreg_env_has_any_spatial <- function(object) {
  any(vapply(
    .unitreg_env_components,
    function(component) .unitreg_env_has_spatial(object, component),
    logical(1)
  ))
}

.unitreg_env_re_info <- function(object, component) {
  object[[paste0("re_info_", component)]]
}

.unitreg_env_open_unit_guard <- function(x, name = "value") {
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

.unitreg_env_mvn_rows <- function(n, Sigma) {
  n <- as.integer(n)
  Sigma <- as.matrix(Sigma)
  q <- nrow(Sigma)

  if (q < 1L || ncol(Sigma) != q || any(!is.finite(Sigma))) {
    stop("Invalid random-effect covariance matrix.", call. = FALSE)
  }

  Sigma <- (Sigma + t(Sigma)) / 2

  ee <- eigen(Sigma, symmetric = TRUE)
  tol <- max(1, max(abs(ee$values))) * 1e-10

  if (min(ee$values) < -tol) {
    stop(
      "Estimated random-effect covariance matrix is not positive semidefinite.",
      call. = FALSE
    )
  }

  vals <- pmax(ee$values, 0)
  A <- ee$vectors %*% diag(sqrt(vals), nrow = q)
  Z <- matrix(stats::rnorm(n * q), nrow = n, ncol = q)

  Z %*% t(A)
}

.unitreg_env_variance_tables <- function(object, component) {
  by_component <- object$random_effects_variance_by_component

  if (!is.null(by_component) && !is.null(by_component[[component]])) {
    return(by_component[[component]])
  }

  if (identical(component, "mu") && !is.null(object$random_effects_variance)) {
    return(object$random_effects_variance)
  }

  NULL
}

.unitreg_env_term_covariance <- function(object, component, term) {
  info <- .unitreg_env_re_info(object, component)
  tabs <- .unitreg_env_variance_tables(object, component)

  if (is.null(info) || is.null(tabs)) {
    stop(
      "Random-effect variance metadata are missing for component '",
      component, "'.",
      call. = FALSE
    )
  }

  if (length(tabs) < term || is.null(tabs[[term]])) {
    stop(
      "Random-effect variance table is missing for component '",
      component, "', term ", term, ".",
      call. = FALSE
    )
  }

  tab <- tabs[[term]]
  qvec <- .unitreg_env_or(info$re_term_size, info$re_term_n_components)
  q <- as.integer(qvec[term])

  if (is.na(q) || q < 1L) {
    stop(
      "Invalid random-effect term size for component '",
      component, "', term ", term, ".",
      call. = FALSE
    )
  }

  sd_rows <- startsWith(rownames(tab), "Std.Dev.")
  sd <- as.numeric(tab[sd_rows, "Estimate"])

  if (length(sd) != q || any(!is.finite(sd)) || any(sd < 0)) {
    stop(
      "Invalid random-effect standard deviations for component '",
      component, "', term ", term, ".",
      call. = FALSE
    )
  }

  if (q == 1L) return(matrix(sd^2, 1L, 1L))

  corr_rows <- startsWith(rownames(tab), "Corr(")
  corr_vals <- as.numeric(tab[corr_rows, "Estimate"])
  need <- q * (q - 1L) / 2L

  if (length(corr_vals) != need || any(!is.finite(corr_vals))) {
    stop(
      "Invalid random-effect correlations for component '",
      component, "', term ", term, ".",
      call. = FALSE
    )
  }

  R <- diag(q)
  idx <- which(lower.tri(R), arr.ind = TRUE)

  for (k in seq_len(nrow(idx))) {
    i <- idx[k, 1L]
    j <- idx[k, 2L]
    R[i, j] <- corr_vals[k]
    R[j, i] <- corr_vals[k]
  }

  D <- diag(sd, nrow = q)
  D %*% R %*% D
}

.unitreg_env_draw_re_vector <- function(object, component) {
  if (!.unitreg_env_has_re(object, component)) return(numeric(0))

  info <- .unitreg_env_re_info(object, component)

  if (is.null(info) || !isTRUE(info$has_random_effects)) {
    stop(
      "Malformed random-effect metadata for component '",
      component, "'.",
      call. = FALSE
    )
  }

  qvec <- .unitreg_env_or(info$re_term_size, info$re_term_n_components)
  n_terms <- as.integer(info$n_re_terms)

  if (length(n_terms) != 1L || is.na(n_terms) || n_terms < 1L) {
    stop(
      "Invalid random-effect term count for component '",
      component, "'.",
      call. = FALSE
    )
  }

  out <- vector("list", n_terms)

  for (term in seq_len(n_terms)) {
    nlev <- as.integer(info$n_re_levels_list[term])
    q <- as.integer(qvec[term])
    Sigma <- .unitreg_env_term_covariance(object, component, term)

    if (is.na(nlev) || nlev < 1L || is.na(q) || q < 1L) {
      stop(
        "Invalid random-effect dimensions for component '",
        component, "', term ", term, ".",
        call. = FALSE
      )
    }

    B <- .unitreg_env_mvn_rows(nlev, Sigma)

    out[[term]] <- as.numeric(t(B))
  }

  unlist(out, use.names = FALSE)
}

.unitreg_env_fixed_eta <- function(object, component) {
  n <- object$nobs

  active <- component %in% c("mu", "phi") ||
    isTRUE(object[[paste0("has_", component)]])

  if (!active) return(rep(0, n))

  dat <- object$obj$env$data
  X <- dat[[paste0("X_", component)]]
  beta <- object$model_coef[[component]]

  if (is.null(X) || is.null(beta) ||
      nrow(X) != n || ncol(X) != length(beta)) {
    stop(
      "Fixed-effect information is inconsistent for component '",
      component, "'.",
      call. = FALSE
    )
  }

  eta <- .unitreg_env_matvec(
    X, beta,
    paste0("fixed predictor for ", component)
  )

  if (component == "mu" && !is.null(dat$offset_mu)) {
    eta <- eta + as.numeric(dat$offset_mu)
  }

  eta
}

.unitreg_env_probabilities <- function(object, eta_p0, eta_p1) {
  n <- length(eta_p0)
  has0 <- isTRUE(object$has_p0)
  has1 <- isTRUE(object$has_p1)

  if (has0 && has1) {
    m <- pmax(0, eta_p0, eta_p1)
    ec <- exp(-m)
    e0 <- exp(eta_p0 - m)
    e1 <- exp(eta_p1 - m)
    den <- ec + e0 + e1

    return(list(
      p0 = e0 / den,
      p1 = e1 / den,
      pc = ec / den
    ))
  }

  if (has0) {
    p0 <- stats::plogis(eta_p0)
    return(list(p0 = p0, p1 = rep(0, n), pc = 1 - p0))
  }

  if (has1) {
    p1 <- stats::plogis(eta_p1)
    return(list(p0 = rep(0, n), p1 = p1, pc = 1 - p1))
  }

  list(
    p0 = rep(0, n),
    p1 = rep(0, n),
    pc = rep(1, n)
  )
}

.unitreg_env_parameters_from_eta <- function(object, eta) {
  mu <- object$link.mu$linkinv(eta$mu)

  phi_link_name <- .unitreg_env_or(
    object$family_object$phi_link_r_name,
    "log"
  )

  phi <- stats::make.link(as.character(phi_link_name)[1L])$linkinv(eta$phi)
  probs <- .unitreg_env_probabilities(object, eta$p0, eta$p1)

  mu <- .unitreg_env_open_unit_guard(mu, "simulated mu")
  phi <- as.numeric(phi)

  if (any(!is.finite(phi)) || any(phi <= 0)) {
    stop("Simulated phi values must be positive and finite.", call. = FALSE)
  }

  list(
    mu = mu,
    phi = phi,
    p0 = probs$p0,
    p1 = probs$p1,
    pc = probs$pc
  )
}

.unitreg_env_draw_marginal_parameters <- function(object) {
  if (.unitreg_env_has_any_spatial(object)) {
    stop(
      "Marginal envelope generation for proper CAR spatial effects is not ",
      "implemented. Use bootstrap = 'conditional' for spatial models.",
      call. = FALSE
    )
  }

  eta <- setNames(
    lapply(
      .unitreg_env_components,
      function(component) .unitreg_env_fixed_eta(object, component)
    ),
    .unitreg_env_components
  )

  for (component in .unitreg_env_components) {
    if (!.unitreg_env_has_re(object, component)) next

    Z <- object$obj$env$data[[paste0("Z_", component)]]
    u <- .unitreg_env_draw_re_vector(object, component)

    if (is.null(Z) || ncol(Z) != length(u)) {
      stop(
        "Random-effect design information is inconsistent for component '",
        component, "'.",
        call. = FALSE
      )
    }

    eta[[component]] <- eta[[component]] +
      .unitreg_env_matvec(
        Z, u,
        paste0("random predictor for ", component)
      )
  }

  .unitreg_env_parameters_from_eta(object, eta)
}

.unitreg_env_draw_response <- function(object, par) {
  n <- length(par$mu)
  u <- stats::runif(n)

  is0 <- u < par$p0
  is1 <- !is0 & u < (par$p0 + par$p1)
  isc <- !(is0 | is1)

  y <- numeric(n)
  y[is0] <- 0
  y[is1] <- 1

  if (any(isc)) {
    y_cont <- get_random_continuous(
      family_code = object$family_object$family_code,
      n = sum(isc),
      mu = par$mu[isc],
      phi = par$phi[isc],
      tau = object$tau
    )

    if (length(y_cont) != sum(isc) || any(!is.finite(y_cont))) {
      stop(
        "Continuous random generator returned invalid values.",
        call. = FALSE
      )
    }

    y[isc] <- .unitreg_env_open_unit_guard(
      y_cont,
      "simulated response"
    )
  }

  y
}

.unitreg_env_simulate_response <- function(object, bootstrap) {
  if (identical(bootstrap, "conditional")) {
    sim <- stats::simulate(
      object,
      nsim = 1L,
      re.form = NULL
    )

    y <- as.numeric(sim[[1L]])

    if (length(y) != object$nobs || any(!is.finite(y))) {
      stop(
        "Conditional simulate() returned an invalid response vector.",
        call. = FALSE
      )
    }

    return(y)
  }

  par <- .unitreg_env_draw_marginal_parameters(object)
  .unitreg_env_draw_response(object, par)
}

.unitreg_env_analysis_data <- function(object) {
  dat <- object$data
  n <- object$nobs
  used_rows <- rownames(object$model)

  if (!is.null(used_rows) &&
      length(used_rows) == n &&
      !is.null(rownames(dat)) &&
      all(used_rows %in% rownames(dat))) {
    dat <- dat[used_rows, , drop = FALSE]
  } else if (nrow(dat) != n) {
    stop(
      "Could not reconstruct the exact analysis data used by the fitted model.",
      call. = FALSE
    )
  }

  dat
}

.unitreg_env_refit_args <- function(object, data, control, start = NULL) {
  args <- list(
    formula = object$formula,
    data = data,
    family = object$family_object,
    phi.formula = object$phi.formula,
    control = control
  )

  if (isTRUE(object$has_p0)) args$p0.formula <- object$p0.formula
  if (isTRUE(object$has_p1)) args$p1.formula <- object$p1.formula

  if (!is.null(object$spatial.formula)) {
    args$spatial <- object$spatial.formula
  }

  if (!is.null(object$spatial.phi.formula)) {
    args$spatial.phi <- object$spatial.phi.formula
  }

  if (!is.null(object$spatial.p0.formula)) {
    args$spatial.p0 <- object$spatial.p0.formula
  }

  if (!is.null(object$spatial.p1.formula)) {
    args$spatial.p1 <- object$spatial.p1.formula
  }

  if (!is.null(object$adjacency)) {
    args$adjacency <- object$adjacency
  }

  if (!is.null(object$weights) &&
      length(object$weights) == object$nobs &&
      any(abs(as.numeric(object$weights) - 1) > 1e-15)) {
    args$weights <- as.numeric(object$weights)
  }

  if (!is.null(object$offset) &&
      length(object$offset) == object$nobs &&
      any(abs(as.numeric(object$offset)) > 1e-15)) {
    args$offset <- as.numeric(object$offset)
  }

  if (!is.null(start)) args$start <- start

  args
}

.unitreg_env_full_start <- function(object) {
  if (!exists(".unitreg_full_par", mode = "function", inherits = TRUE)) {
    return(NULL)
  }

  tryCatch(
    .unitreg_full_par(object),
    error = function(e) NULL
  )
}

.unitreg_env_preflight <- function(object) {
  if (.unitreg_env_has_any_spatial(object)) {
    return(list(
      ok = FALSE,
      stage = "spatial",
      error = paste0(
        "marginal bootstrap does not currently generate new proper CAR ",
        "spatial effects"
      )
    ))
  }

  par <- tryCatch(
    .unitreg_env_draw_marginal_parameters(object),
    error = function(e) e
  )

  if (inherits(par, "error")) {
    return(list(
      ok = FALSE,
      stage = "parameters",
      error = conditionMessage(par)
    ))
  }

  y <- tryCatch(
    .unitreg_env_draw_response(object, par),
    error = function(e) e
  )

  if (inherits(y, "error")) {
    return(list(
      ok = FALSE,
      stage = "response",
      error = conditionMessage(y)
    ))
  }

  list(ok = TRUE, stage = "ok", error = NULL)
}

.unitreg_env_plot <- function(x, main = NULL, ...) {
  d <- x$data

  yr <- range(
    c(d$Lower, d$Upper, d$Observed),
    finite = TRUE
  )

  xr <- range(d$Theoretical, finite = TRUE)

  graphics::plot(
    d$Theoretical,
    d$Observed,
    type = "n",
    xlim = xr,
    ylim = yr,
    xlab = paste(
      "Theoretical",
      x$distribution_label,
      "quantiles"
    ),
    ylab = x$ylab,
    main = .unitreg_env_or(
      main,
      paste0(
        as.character(x$family)[1L],
        " [",
        x$bootstrap,
        " bootstrap]"
      )
    ),
    bty = "o",
    ...
  )

  graphics::grid(
    col = "gray92",
    lty = 1,
    lwd = 1.1
  )

  graphics::lines(
    d$Theoretical,
    d$Lower,
    lty = 1,
    col = "gray35",
    lwd = 1.4
  )

  graphics::lines(
    d$Theoretical,
    d$Upper,
    lty = 1,
    col = "gray35",
    lwd = 1.4
  )

  graphics::lines(
    d$Theoretical,
    d$Median,
    lty = 2,
    col = "black",
    lwd = 1.2
  )

  graphics::points(
    d$Theoretical,
    d$Observed,
    pch = 3,
    col = "black",
    cex = 0.8
  )

  invisible(x)
}

#' Simulated Residual Envelope for unitregTMB Models
#'
#' @description
#' Constructs a pointwise simulated envelope for residual diagnostics from a
#' fitted \code{unitregTMB} model. Each bootstrap data set is generated from the
#' fitted model, the model is refitted, residuals are recomputed, and their
#' ordered values are used to obtain the lower, median, and upper pointwise
#' envelope curves.
#'
#' @param object A fitted object of class \code{unitregTMB}.
#' @param nsim Integer giving the target number of valid bootstrap refits.
#'   The default is \code{50}. For final diagnostic figures, larger values such
#'   as 300--500 are recommended when computationally feasible.
#' @param bootstrap Character string specifying how latent ordinary random
#'   effects are handled when bootstrap responses are generated. With
#'   \code{"marginal"} (the default), new Gaussian random effects are drawn from
#'   the fitted covariance matrices before each response is generated. With
#'   \code{"conditional"}, responses are generated conditionally on the fitted
#'   latent modes. For models without ordinary random effects the two options
#'   coincide. Marginal generation of new proper CAR spatial effects is not
#'   currently implemented; spatial models should use
#'   \code{bootstrap = "conditional"}.
#' @param halfnormal Logical. If \code{TRUE} and
#'   \code{resid.type = "quantile"}, absolute randomized quantile residuals are
#'   plotted against theoretical half-normal quantiles. If \code{FALSE}, signed
#'   residuals are plotted against standard normal quantiles. This argument is
#'   ignored for Cox--Snell residuals.
#' @param plot Logical. If \code{TRUE}, draw the envelope.
#' @param level Numeric scalar in \code{(0, 1)} giving the pointwise envelope
#'   level. The default is \code{0.95}.
#' @param resid.type Character string specifying the residual type. Available
#'   choices are \code{"quantile"} and \code{"cox-snell"}.
#' @param ncpus Positive integer giving the number of workers used for the outer
#'   bootstrap. Each bootstrap refit is forced to use a single internal
#'   \code{unitregTMB} worker when \code{ncpus > 1}, avoiding nested
#'   parallelism.
#' @param show.progress Logical. If \code{TRUE}, show progress toward the target
#'   number of valid bootstrap refits.
#' @param max_attempts Positive integer giving the maximum number of bootstrap
#'   data sets that may be attempted. If \code{NULL}, twice \code{nsim} is used.
#' @param seed Optional numeric scalar used to initialize random-number
#'   generation. Supplying a seed makes a serial run reproducible. Parallel
#'   runs use independent reproducible worker streams for a fixed
#'   \code{ncpus}.
#' @param grad.tol Positive numeric scalar giving the maximum accepted absolute
#'   final gradient for a bootstrap refit. The default is \code{1e-3}. Set to
#'   \code{NULL} to disable the gradient filter and rely only on the optimizer
#'   convergence code and finite-solution checks.
#' @param strict Logical. If \code{TRUE}, stop when fewer than \code{nsim}
#'   valid refits are obtained after \code{max_attempts}. If \code{FALSE}, use
#'   the available valid refits (at least two) and issue a warning.
#' @param retry.start Logical. If \code{TRUE}, a failed bootstrap refit that was
#'   initialized at the original fitted parameter values is retried once using
#'   the package default starting values.
#' @param warm.start Logical. If \code{TRUE}, use the fitted parameter values
#'   from \code{object} as starting values for the first refit attempt.
#' @param quantile.type Integer from 1 to 9 passed to
#'   \code{stats::quantile()} when the pointwise envelope limits are computed.
#'   The default is \code{7}, matching the base R default.
#' @param keep.sim Logical. If \code{TRUE}, retain the matrix of ordered
#'   simulated residuals in the returned object.
#' @param ... Additional graphical arguments passed to \code{graphics::plot()}.
#'
#' @details
#' The marginal algorithm for ordinary mixed-effects models is
#'
#' \deqn{
#'   b^{*} \sim N(0,\widehat{\Sigma}), \qquad
#'   Y^{*}\mid b^{*} \sim F(\widehat{\theta},b^{*}),
#' }
#'
#' followed by a complete refit of the model to each simulated response. Thus,
#' the residual envelope reflects both response variation and the estimation
#' step. In contrast, \code{bootstrap = "conditional"} generates responses
#' conditional on the fitted latent modes and then refits the model.
#'
#' Failed refits do not enter the envelope. A refit is accepted only when the
#' optimizer reports convergence, the optimized objective and parameter vector
#' are finite, and, when \code{grad.tol} is not \code{NULL}, the final maximum
#' absolute gradient does not exceed the requested tolerance. The returned
#' \code{bootstrap_diagnostics} component records accepted and rejected
#' attempts.
#'
#' The envelope is pointwise rather than simultaneous. With a 95 percent
#' envelope, for example, each tail is estimated from 2.5 percent of the valid
#' bootstrap replicates. Consequently, a small value of \code{nsim} can make
#' the extreme envelope limits unstable even when the fitted model is correct.
#'
#' For zero- and/or one-augmented models, randomized quantile residuals involve
#' additional random uniforms at boundary observations. The \code{seed}
#' argument controls this randomization together with the bootstrap simulation.
#'
#' @return
#' Invisibly returns an object of class \code{envelope.unitregTMB}, a list with
#' the following main components:
#' \describe{
#'   \item{\code{data}}{Data frame containing theoretical quantiles, observed
#'   ordered residuals, and the lower, median, and upper envelope curves.}
#'   \item{\code{bootstrap}}{Bootstrap type used.}
#'   \item{\code{generates_new_random_effects}}{Logical indicating whether new
#'   ordinary random effects were generated.}
#'   \item{\code{nsim.requested}, \code{nsim.valid}}{Requested and retained
#'   numbers of bootstrap refits.}
#'   \item{\code{attempts}, \code{failures}, \code{retries.used}}{Bootstrap
#'   fitting diagnostics.}
#'   \item{\code{bootstrap_diagnostics}}{Data frame with the status,
#'   convergence code, maximum gradient, retry indicator, and failure reason
#'   for each attempted bootstrap data set.}
#'   \item{\code{simulated_residuals}}{Matrix of ordered simulated residuals,
#'   included only when \code{keep.sim = TRUE}.}
#' }
#'
#' @examples
#' \dontrun{
#' library(unitregTMB)
#'
#' dat <- sim_unitreg(
#'   n_id = 40,
#'   n_rep = 4,
#'   family = beta_fam(model_for = "mode"),
#'   formula_mu = ~ x + (1 | id),
#'   formula_phi = ~ x,
#'   beta_mu = c(-0.5, 0.8),
#'   beta_phi = c(log(10), 0.30),
#'   Sigma_mu = matrix(0.30^2),
#'   cov_config = list(
#'     x = list(dist = "normal", mean = 0, sd = 1)
#'   ),
#'   seed = 123
#' )
#'
#' fit <- unitregTMB(Y ~ x + (1 | id),
#'                   phi.formula = ~ x,
#'                   family = beta_fam(model_for = "mode"),
#'                   data = dat)
#'
#' env <- plot_envelope(fit,
#'                      nsim = 100,
#'                      bootstrap = "marginal",
#'                      seed = 123)
#'
#' env$nsim.valid
#' head(env$bootstrap_diagnostics)
#'
#' # Conditional bootstrap as a complementary diagnostic:
#' plot_envelope(fit, 
#'               nsim = 100,
#'               bootstrap = "conditional",
#'               seed = 123)
#' }
#'
#' @seealso \code{residuals.unitregTMB()}, \code{simulate.unitregTMB()}
#'
#' @importFrom regcore plot_envelope
#' @method plot_envelope unitregTMB
#' @export
plot_envelope.unitregTMB <- function(object, 
                                     nsim = 50L, 
                                     bootstrap = c("marginal", "conditional"),
                                     halfnormal = TRUE, 
                                     plot = TRUE, 
                                     level = 0.95, 
                                     resid.type = c("quantile", "cox-snell"),
                                     ncpus = 1L, 
                                     show.progress = TRUE, 
                                     max_attempts = NULL, 
                                     seed = NULL, 
                                     grad.tol = 1e-3, 
                                     strict = TRUE, 
                                     retry.start = TRUE, 
                                     warm.start = TRUE, 
                                     quantile.type = 7L, 
                                     keep.sim = FALSE, ...) {

  if (!inherits(object, "unitregTMB")) {
    stop("'object' must inherit from 'unitregTMB'.", call. = FALSE)
  }

  bootstrap <- match.arg(bootstrap)
  resid.type <- match.arg(resid.type)

  check_integer <- function(x, name, lower = 1L) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower || x != as.integer(x)) {
      stop("'", name, "' must be an integer >= ", lower, ".", call. = FALSE)
    }

    as.integer(x)
  }

  nsim <- check_integer(nsim, "nsim", 2L)
  ncpus <- check_integer(ncpus, "ncpus", 1L)
  quantile.type <- check_integer(quantile.type, "quantile.type", 1L)

  if (quantile.type > 9L) {
    stop("'quantile.type' must be between 1 and 9.", call. = FALSE)
  }

  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("'level' must be one finite number in (0, 1).", call. = FALSE)
  }

  if (!is.null(grad.tol)) {
    if (!is.numeric(grad.tol) || length(grad.tol) != 1L || !is.finite(grad.tol) || grad.tol <= 0) {
      stop("'grad.tol' must be NULL or one positive finite number.", call. = FALSE)
    }
  }

  if (is.null(max_attempts)) {
    max_attempts <- max(nsim, as.integer(ceiling(2 * nsim)))
  }

  max_attempts <- check_integer(max_attempts, "max_attempts", nsim)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("'seed' must be NULL or one finite number.", call. = FALSE)
    }

    seed <- as.integer(seed)

    if (is.na(seed)) {
      stop("'seed' is outside the supported integer range.", call. = FALSE)
    }
  }

  if (identical(resid.type, "cox-snell")) {
    halfnormal <- FALSE
  }

  has_re <- .unitreg_env_has_any_re(object)

  if (identical(bootstrap, "marginal") && .unitreg_env_has_any_spatial(object)) {
    stop("bootstrap = 'marginal' currently supports ordinary Gaussian random ",
         "effects but not proper CAR spatial effects. Use ",
         "bootstrap = 'conditional' for this model.", call. = FALSE)
  }

  if (identical(bootstrap, "marginal") && !has_re && show.progress) {
    message("No ordinary random effects detected: marginal and conditional ",
            "generation coincide.")
  }

  if (identical(bootstrap, "marginal")) {
    preflight <- .unitreg_env_preflight(object)

    if (!isTRUE(preflight$ok)) {
      stop("Marginal bootstrap preflight failed at stage '", preflight$stage, "': ",
           preflight$error, ".", call. = FALSE)
    }
  }

  alpha <- (1 - level) / 2
  
  dat0 <- .unitreg_env_analysis_data(object)
  lhs <- object$formula[[2L]]

  if (!is.symbol(lhs)) {
    stop("The envelope currently requires a simple response variable on the ", 
         "left-hand side.", call. = FALSE)
  }

  response_name <- as.character(lhs)

  if (!response_name %in% names(dat0)) {
    stop("Response variable '", response_name, "' is not present in the stored analysis data.", call. = FALSE)
  }

  control <- object$control
  control$sd.report <- FALSE
  control$get.joint.precision <- FALSE
  control$verbose <- FALSE

  if (ncpus > 1L) {
    control$parallel <- list(n = 1L, autopar = FALSE)
  }

  start0 <- if (isTRUE(warm.start)) {
    .unitreg_env_full_start(object)
  } else {
    NULL
  }

  if (!is.null(seed)) set.seed(seed)

  observed_raw <- tryCatch(as.numeric(stats::residuals(object, type = resid.type, re.form = NULL)),
                           error = function(e) e)

  if (inherits(observed_raw, "error")) {
    stop("Observed residual calculation failed: ", conditionMessage(observed_raw), call. = FALSE)
  }

  n <- object$nobs

  if (length(observed_raw) != n || any(!is.finite(observed_raw))) {
    stop("Observed residuals are invalid.", call. = FALSE)
  }

  bootstrap_seed <- NULL

  if (!is.null(seed)) {
    bootstrap_seed <- as.integer((as.double(seed) %% 2147483646) + 1)
    set.seed(bootstrap_seed)
  }

  refit_one <- function(attempt) {
    fail <- function(stage, reason, convergence = NA_integer_, max_gradient = NA_real_, retried = FALSE) {

      list(ok = FALSE, residuals = NULL, stage = stage, reason = reason, convergence = convergence,
           max_gradient = max_gradient, retried = retried)
    }

    y_sim <- tryCatch(
      .unitreg_env_simulate_response(object, bootstrap),
      error = function(e) e)

    if (inherits(y_sim, "error")) {
      return(fail("simulate", conditionMessage(y_sim)))
    }

    if (length(y_sim) != n || any(!is.finite(y_sim))) {
      return(fail("simulate", "generated response vector is invalid"))
    }

    dat <- dat0
    dat[[response_name]] <- y_sim

    fit_once <- function(use_start, retried) {
      args <- .unitreg_env_refit_args(object = object, data = dat, control = control,
               start = if (use_start) start0 else NULL)

      fit_sim <- tryCatch(do.call(unitregTMB, args), error = function(e) e)

      if (inherits(fit_sim, "error")) {
        return(fail("refit", conditionMessage(fit_sim), retried = retried))
      }

      convergence <- as.integer(.unitreg_env_or(fit_sim$convergence, NA_integer_))[1L]

      max_gradient <- as.numeric(.unitreg_env_or(fit_sim$max_gradient, NA_real_))[1L]

      objective <- tryCatch(as.numeric(fit_sim$opt$objective)[1L], error = function(e) NA_real_)

      opt_par <- tryCatch(as.numeric(fit_sim$opt$par), error = function(e) numeric(0))

      if (!isTRUE(convergence == 0L)) {
        return(fail("convergence", paste0("optimizer convergence code = ", convergence),
                    convergence, max_gradient, retried))
      }

      if (!is.finite(objective) || !length(opt_par) || any(!is.finite(opt_par))) {
        return(fail("solution", "non-finite optimized objective or parameter vector",
                    convergence, max_gradient, retried))
      }

      if (!is.null(grad.tol)) {
        if (!is.finite(max_gradient)) {
          return(fail("gradient", "final maximum absolute gradient is non-finite",
                      convergence, max_gradient, retried))
        }

        if (max_gradient > grad.tol) {
          return(fail("gradient", paste0("max_gradient = ", format(max_gradient, scientific = TRUE),
                      " > grad.tol = ", format(grad.tol, scientific = TRUE)), convergence,
                      max_gradient, retried))
        }
      }

      rr <- tryCatch(as.numeric(stats::residuals(fit_sim, type = resid.type, re.form = NULL)),
                     error = function(e) e)

      if (inherits(rr, "error")) {
        return(fail("residuals", conditionMessage(rr), convergence, 
                    max_gradient, retried))
      }

      if (length(rr) != n || any(!is.finite(rr))) {
        return(fail("residuals", "invalid residual vector", convergence, max_gradient, retried))
      }

      list(ok = TRUE, residuals = rr, stage = "ok", reason = "", 
           convergence = convergence, max_gradient = max_gradient, retried = retried)
    }

    first <- fit_once(use_start = !is.null(start0), retried = FALSE)

    if (isTRUE(first$ok) || !isTRUE(retry.start) || is.null(start0)) {
      return(first)
    }

    second <- fit_once(use_start = FALSE, retried = TRUE)

    if (!isTRUE(second$ok)) {
      second$reason <- paste0("warm-start failure [", first$stage, ": ", first$reason, "]; retry failure [",
                              second$stage, ": ", second$reason, "]")
    }
    
    second
  }

  attempts <- list()
  n_valid <- 0L

  if (ncpus == 1L) {
    if (show.progress) {
      cat("Simulating ", bootstrap, " residual envelope (target ", nsim, 
          " valid fits; up to ", max_attempts, " attempts)...\n", sep = "")
    }

    pb <- NULL

    if (show.progress) {
      pb <- utils::txtProgressBar(min = 0, max = nsim, style = 3)

      on.exit(try(close(pb), silent = TRUE), add = TRUE)
    }

    for (attempt in seq_len(max_attempts)) {
      z <- refit_one(attempt)
      attempts[[length(attempts) + 1L]] <- z

      if (isTRUE(z$ok)) {
        n_valid <- n_valid + 1L

        if (show.progress) {
          utils::setTxtProgressBar(pb, min(n_valid, nsim))
        }
      }

      if (n_valid >= nsim) break
    }

    if (show.progress) {
      close(pb)
      pb <- NULL
      cat("\n")
    }

  } else {
    if (show.progress) {
      cat("Simulating ", bootstrap, " residual envelope with ", 
          ncpus, " workers (target ", nsim, " valid fits; up to ",
          max_attempts, " attempts)...\n", sep = "")
    }

    cl <- parallel::makeCluster(ncpus)

    tryCatch(
      {
        parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(unitregTMB)))

        if (!is.null(bootstrap_seed)) {
          parallel::clusterSetRNGStream(
            cl,
            iseed = bootstrap_seed
          )
        }

        pb <- NULL

        if (show.progress) {
          pb <- utils::txtProgressBar(min = 0, max = nsim, style = 3)
        }

        next_attempt <- 1L

        while (next_attempt <= max_attempts &&
               n_valid < nsim) {

          batch_end <- min(max_attempts, next_attempt + ncpus - 1L)

          idx <- seq.int(next_attempt, batch_end)

          batch <- parallel::parLapply(cl, idx, refit_one)

          attempts <- c(attempts, batch)

          n_valid <- sum(vapply(attempts, function(z) isTRUE(z$ok), logical(1)))

          if (show.progress) {
            utils::setTxtProgressBar(pb, min(n_valid, nsim))
          }

          next_attempt <- batch_end + 1L
        }

        if (show.progress) {
          close(pb)
          pb <- NULL
          cat("\n")
        }
      },
      finally = {
        parallel::stopCluster(cl)
      }
    )
  }

  ok <- vapply(attempts, function(z) isTRUE(z$ok), logical(1))

  diagnostics <- data.frame(
    attempt = seq_along(attempts),
    ok = ok,
    stage = vapply(attempts, function(z) .unitreg_env_or(z$stage, ""), character(1)),
    reason = vapply(attempts, function(z) .unitreg_env_or(z$reason, ""), character(1)),
    convergence = vapply(attempts, function(z) {
        value <- z$convergence
        if (!length(value) || is.na(value[1L])) {
          NA_integer_
        } else {
          as.integer(value[1L])
        }
      },
      integer(1)
    ),
    max_gradient = vapply(attempts,
      function(z) {
        value <- z$max_gradient

        if (!length(value) || is.na(value[1L])) {
          NA_real_
        } else {
          as.numeric(value[1L])
        }
      },
      numeric(1)
    ),
    retried = vapply(attempts, function(z) isTRUE(z$retried), logical(1)),
    stringsAsFactors = FALSE)

  valid <- attempts[ok]

  if (length(valid) < nsim) {
    tab_fail <- sort(table(diagnostics$stage[!diagnostics$ok]), decreasing = TRUE)

    failure_text <- if (length(tab_fail)) {
      paste(paste0(names(tab_fail), " [n=", as.integer(tab_fail), "]"), collapse = "; ")
    } else {
      "none"
    }

    msg <- paste0("Only ", length(valid), " valid bootstrap fits were obtained after ", length(attempts),
                  " attempts (target = ", nsim, "). Failures: ", failure_text, ".")

    if (isTRUE(strict)) {
      stop(msg, call. = FALSE)
    }

    warning(msg, call. = FALSE)
  }

  if (length(valid) < 2L) {
    stop("Fewer than two valid bootstrap fits remain.", call. = FALSE)
  }

  valid <- valid[
    seq_len(min(nsim, length(valid)))
  ]

  res_sim <- do.call(cbind, lapply(valid, `[[`, "residuals"))

  if (identical(resid.type, "cox-snell")) {
    ylab <- "Ordered Cox-Snell residuals"
    distribution_label <- "Exponential"

    observed <- sort(observed_raw)

    simulated_ordered <- apply(res_sim, 2L, sort, na.last = NA)

    p <- (seq_len(n) - 3 / 8) / (n + 1 / 4)

    theoretical <- stats::qexp(p)

  } else if (isTRUE(halfnormal)) {
    ylab <- "Ordered absolute randomized quantile residuals"
    distribution_label <- "Half-Normal"

    observed <- sort(abs(observed_raw))

    simulated_ordered <- apply(res_sim, 2L, function(x) sort(abs(x), na.last = NA))

    p <- (seq_len(n) + n - 1 / 8) / (2 * n + 0.5)

    theoretical <- stats::qnorm(p)

  } else {
    ylab <- "Ordered randomized quantile residuals"
    distribution_label <- "Normal"

    observed <- sort(observed_raw)

    simulated_ordered <- apply(res_sim, 2L, sort, na.last = NA)

    p <- (seq_len(n) - 3 / 8) / (n + 1 / 4)

    theoretical <- stats::qnorm(p)
  }

  lower <- apply(simulated_ordered, 1L, stats::quantile, probs = alpha,
                 type = quantile.type, names = FALSE, na.rm = TRUE)

  median_curve <- apply(simulated_ordered, 1L, stats::median, na.rm = TRUE)

  upper <- apply(simulated_ordered, 1L, stats::quantile, probs = 1 - alpha, 
                 type = quantile.type, names = FALSE, na.rm = TRUE)

  envelope_data <- data.frame(Theoretical = theoretical,
                              Observed = observed,
                              Lower = lower,
                              Median = median_curve,
                              Upper = upper,
                              row.names = NULL)

  ans <- list(data = envelope_data,
              family = object$family,
              bootstrap = bootstrap,
              generates_new_random_effects = identical(bootstrap, "marginal") && has_re,
              resid.type = resid.type,
              halfnormal = isTRUE(halfnormal),
              level = level,
              quantile.type = quantile.type,
              nsim.requested = nsim,
              nsim.valid = ncol(simulated_ordered),
              attempts = length(attempts),
              failures = sum(!ok),
              retries.used = sum(vapply(valid, function(z) isTRUE(z$retried), logical(1))),
              grad.tol = grad.tol,
              bootstrap_diagnostics = diagnostics,
              ylab = ylab,
              distribution_label = distribution_label)

  if (isTRUE(keep.sim)) {
    ans$simulated_residuals <- simulated_ordered
  }

  class(ans) <- "envelope.unitregTMB"

  if (isTRUE(plot)) {
    .unitreg_env_plot(ans, ...)
  }

  invisible(ans)
}
