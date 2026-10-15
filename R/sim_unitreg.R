#' Test Whether a Simulation Formula Is Inactive
#'
#' @description
#' Internal helper used by \code{sim_unitreg()} to identify \code{NULL},
#' \code{~0}, or \code{~-1} component formulas.
#'
#' @param formula A formula or \code{NULL}.
#'
#' @return Logical scalar.
#'
#' @examples
#' unitregTMB:::.sim_unitreg_is_zero_formula(~ 0)
#'
#' @noRd
.sim_unitreg_is_zero_formula <- function(formula) {
  if (is.null(formula)) return(TRUE)
  tt <- stats::terms(formula)
  identical(attr(tt, "intercept"), 0L) &&
    length(attr(tt, "term.labels")) == 0L
}

#' Generate a Covariate for sim_unitreg
#'
#' @description
#' Internal helper that generates one covariate from a configuration supplied
#' to \code{sim_unitreg()}.
#'
#' @param n Number of values.
#' @param conf Configuration list or \code{NULL}.
#' @param name Covariate name used in diagnostics.
#'
#' @return A numeric vector or factor.
#'
#' @examples
#' unitregTMB:::.sim_unitreg_generate_covariate(
#'   5, list(dist = "normal", mean = 0, sd = 1), "x"
#' )
#'
#' @noRd
.sim_unitreg_generate_covariate <- function(n, conf, name) {
  if (is.null(conf)) {
    return(stats::rnorm(n))
  }

  if (!is.list(conf) || is.null(conf$dist) ||
      length(conf$dist) != 1L || is.na(conf$dist)) {
    stop(
      "cov_config[['", name,
      "']] must be a list with one 'dist' element.",
      call. = FALSE
    )
  }

  dist <- match.arg(
    tolower(as.character(conf$dist)),
    c("normal", "binomial", "uniform", "poisson", "factor", "categorical")
  )

  if (dist == "normal") {
    mean <- if (is.null(conf$mean)) 0 else conf$mean
    sd <- if (is.null(conf$sd)) 1 else conf$sd

    if (!is.numeric(mean) || length(mean) != 1L || !is.finite(mean) ||
        !is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd < 0) {
      stop(
        "Normal covariate '", name,
        "' requires finite 'mean' and non-negative finite 'sd'.",
        call. = FALSE
      )
    }

    return(stats::rnorm(n, mean = mean, sd = sd))
  }

  if (dist == "binomial") {
    size <- if (is.null(conf$size)) 1 else conf$size
    prob <- if (is.null(conf$prob)) 0.5 else conf$prob

    if (!is.numeric(size) || length(size) != 1L || !is.finite(size) ||
        size < 0 || size != as.integer(size) ||
        !is.numeric(prob) || length(prob) != 1L || !is.finite(prob) ||
        prob < 0 || prob > 1) {
      stop(
        "Binomial covariate '", name,
        "' requires non-negative integer 'size' and 'prob' in [0, 1].",
        call. = FALSE
      )
    }

    return(stats::rbinom(n, size = as.integer(size), prob = prob))
  }

  if (dist == "uniform") {
    min <- if (is.null(conf$min)) 0 else conf$min
    max <- if (is.null(conf$max)) 1 else conf$max

    if (!is.numeric(min) || length(min) != 1L || !is.finite(min) ||
        !is.numeric(max) || length(max) != 1L || !is.finite(max) ||
        min > max) {
      stop(
        "Uniform covariate '", name,
        "' requires finite 'min' <= 'max'.",
        call. = FALSE
      )
    }

    return(stats::runif(n, min = min, max = max))
  }

  if (dist == "poisson") {
    lambda <- if (is.null(conf$lambda)) 1 else conf$lambda

    if (!is.numeric(lambda) || length(lambda) != 1L ||
        !is.finite(lambda) || lambda < 0) {
      stop(
        "Poisson covariate '", name,
        "' requires non-negative finite 'lambda'.",
        call. = FALSE
      )
    }

    return(stats::rpois(n, lambda = lambda))
  }

  lev <- conf$levels

  if (is.null(lev) || !length(lev) || anyNA(lev)) {
    stop(
      "Factor/categorical covariate '", name,
      "' requires non-missing 'levels'.",
      call. = FALSE
    )
  }

  prob <- conf$prob

  if (!is.null(prob)) {
    if (!is.numeric(prob) || length(prob) != length(lev) ||
        any(!is.finite(prob)) || any(prob < 0) || sum(prob) <= 0) {
      stop(
        "Probabilities for factor/categorical covariate '", name,
        "' must be non-negative, finite, and match 'levels'.",
        call. = FALSE
      )
    }
  }

  factor(
    sample(lev, size = n, replace = TRUE, prob = prob),
    levels = lev
  )
}

#' Validate a Random-Effect Covariance Matrix
#'
#' @description
#' Internal helper that selects and validates the covariance matrix associated
#' with one random-effects term in \code{sim_unitreg()}.
#'
#' @param Sigma A covariance matrix or list of covariance matrices.
#' @param i Term index.
#' @param n_terms Number of random-effects terms.
#' @param q Number of random coefficients in the term.
#' @param component Distributional component name.
#'
#' @return A symmetric positive-definite covariance matrix.
#'
#' @examples
#' unitregTMB:::.sim_unitreg_sigma_for_term(
#'   matrix(0.25, 1, 1), 1, 1, 1, "mu"
#' )
#'
#' @noRd
.sim_unitreg_sigma_for_term <- function(Sigma, i, n_terms, q, component) {
  if (is.null(Sigma)) {
    stop(
      "Random effects are present in component '", component,
      "', but its covariance matrix was not supplied.",
      call. = FALSE
    )
  }

  S <- if (is.list(Sigma) && !is.matrix(Sigma)) {
    if (length(Sigma) != n_terms) {
      stop(
        "The covariance list for component '", component,
        "' must contain one matrix per random-effects term.",
        call. = FALSE
      )
    }
    Sigma[[i]]
  } else {
    Sigma
  }

  S <- as.matrix(S)

  if (nrow(S) != q || ncol(S) != q) {
    stop(
      "Random-effects covariance matrix for component '", component,
      "', term ", i, " must be ", q, " x ", q, ".",
      call. = FALSE
    )
  }

  if (any(!is.finite(S))) {
    stop(
      "Random-effects covariance matrix for component '", component,
      "', term ", i, " must contain only finite values.",
      call. = FALSE
    )
  }

  if (!isTRUE(all.equal(S, t(S), tolerance = 1e-10))) {
    stop(
      "Random-effects covariance matrix for component '", component,
      "', term ", i, " must be symmetric.",
      call. = FALSE
    )
  }

  tryCatch(
    chol(S),
    error = function(e) {
      stop(
        "Random-effects covariance matrix for component '", component,
        "', term ", i, " is not positive definite.",
        call. = FALSE
      )
    }
  )

  S
}

#' Generate Ordinary Random-Effect Contributions
#'
#' @description
#' Internal helper that simulates Gaussian random coefficients and their
#' contribution to one distributional predictor in \code{sim_unitreg()}.
#'
#' @param formula Component formula.
#' @param data Simulated data frame.
#' @param Sigma Random-effect covariance specification.
#' @param component Distributional component name.
#'
#' @return A list containing the predictor contribution and generated random
#'   effects.
#'
#' @examples
#' d <- data.frame(id = factor(rep(1:3, each = 2)), x = rnorm(6))
#' unitregTMB:::.sim_unitreg_random_contribution(
#'   ~ x + (1 | id), d, matrix(0.2^2), "mu"
#' )
#'
#' @noRd
.sim_unitreg_random_contribution <- function(
    formula,
    data,
    Sigma,
    component) {

  if (is.null(formula)) {
    return(list(
      contribution = rep(0, nrow(data)),
      effects = list()
    ))
  }

  bars <- reformulas::findbars(formula)

  if (!length(bars)) {
    return(list(
      contribution = rep(0, nrow(data)),
      effects = list()
    ))
  }

  ans <- rep(0, nrow(data))
  effects <- vector("list", length(bars))

  for (i in seq_along(bars)) {
    bar <- bars[[i]]
    grp_expr <- bar[[3L]]

    grp <- tryCatch(
      eval(grp_expr, envir = data, enclos = environment(formula)),
      error = function(e) {
        stop(
          "Could not evaluate grouping expression '",
          paste(deparse(grp_expr), collapse = " "),
          "' for component '", component, "'.",
          call. = FALSE
        )
      }
    )

    if (anyNA(grp)) {
      stop(
        "Grouping expression for component '", component,
        "' produced missing values.",
        call. = FALSE
      )
    }

    grp <- factor(grp)
    lev <- levels(grp)
    idx <- as.integer(grp)
    lhs_expr <- bar[[2L]]

    lhs_formula <- stats::as.formula(
      paste("~", paste(deparse(lhs_expr), collapse = " ")),
      env = environment(formula)
    )

    X_re <- stats::model.matrix(lhs_formula, data)
    q <- ncol(X_re)

    S <- .sim_unitreg_sigma_for_term(
      Sigma = Sigma,
      i = i,
      n_terms = length(bars),
      q = q,
      component = component
    )

    E <- matrix(
      stats::rnorm(length(lev) * q),
      nrow = length(lev),
      ncol = q
    )

    B <- E %*% chol(S)
    colnames(B) <- colnames(X_re)
    rownames(B) <- lev

    ans <- ans + rowSums(
      X_re * B[idx, , drop = FALSE]
    )

    effects[[i]] <- B
  }

  names(effects) <- make.unique(vapply(
    bars,
    function(bar) paste(deparse(bar[[3L]]), collapse = ""),
    character(1)
  ))

  list(
    contribution = ans,
    effects = effects
  )
}

#' Prepare an Adjacency Matrix for Spatial Simulation
#'
#' @description
#' Internal helper that validates and aligns an adjacency matrix with the region
#' levels used by a proper CAR field in \code{sim_unitreg()}.
#'
#' @param adjacency Square adjacency matrix.
#' @param region_levels Character vector of region levels.
#'
#' @return A binary symmetric adjacency matrix aligned to \code{region_levels}.
#'
#' @examples
#' W <- matrix(c(0, 1, 1, 0), 2, 2)
#' dimnames(W) <- list(c("A", "B"), c("A", "B"))
#' unitregTMB:::.sim_unitreg_prepare_adjacency(W, c("A", "B"))
#'
#' @noRd
.sim_unitreg_prepare_adjacency <- function(adjacency, region_levels) {
  if (is.null(adjacency)) {
    stop("A spatial simulation requires 'adjacency'.", call. = FALSE)
  }

  W <- as.matrix(adjacency)

  if (nrow(W) != ncol(W)) {
    stop("'adjacency' must be square.", call. = FALSE)
  }

  if (any(!is.finite(W)) || any(W < 0)) {
    stop(
      "'adjacency' must contain finite non-negative values.",
      call. = FALSE
    )
  }

  if (any(abs(diag(W)) > 1e-12)) {
    stop("The diagonal of 'adjacency' must be zero.", call. = FALSE)
  }

  if (!isTRUE(all.equal(W, t(W), tolerance = 1e-12))) {
    stop("'adjacency' must be symmetric.", call. = FALSE)
  }

  W[W > 0] <- 1
  rn <- rownames(W)
  cn <- colnames(W)

  if (!is.null(rn) && !is.null(cn)) {
    if (anyDuplicated(rn) || anyDuplicated(cn)) {
      stop("Row and column names of 'adjacency' must be unique.", call. = FALSE)
    }

    if (!identical(rn, cn)) {
      stop(
        "Row and column names of 'adjacency' must be identical.",
        call. = FALSE
      )
    }

    missing <- setdiff(region_levels, rn)

    if (length(missing)) {
      stop(
        "Spatial region(s) absent from 'adjacency': ",
        paste(missing, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    W <- W[region_levels, region_levels, drop = FALSE]
  } else {
    if (nrow(W) != length(region_levels)) {
      stop(
        "Unnamed 'adjacency' has ", nrow(W),
        " rows, but the spatial variable has ",
        length(region_levels), " levels.",
        call. = FALSE
      )
    }

    dimnames(W) <- list(region_levels, region_levels)
  }

  if (any(rowSums(W) <= 0)) {
    stop("Proper CAR simulation requires every region to have at least one neighbor.", call. = FALSE)
  }

  W
}

#' Generate a Proper CAR Contribution
#'
#' @description
#' Internal helper that simulates a proper CAR Gaussian field using the same
#' precision structure as the \code{unitregTMB} TMB objective.
#'
#' @param formula One-sided spatial formula.
#' @param data Simulated data frame.
#' @param adjacency Adjacency matrix.
#' @param scale Positive CAR scale parameter.
#' @param rho Proper CAR dependence parameter.
#' @param component Distributional component name.
#'
#' @return A list containing the observation-level contribution, region effects,
#'   and precision matrix.
#'
#' @examples
#' d <- data.frame(region = factor(rep(c("A", "B"), each = 2)))
#' W <- matrix(c(0, 1, 1, 0), 2, 2)
#' dimnames(W) <- list(c("A", "B"), c("A", "B"))
#' unitregTMB:::.sim_unitreg_spatial_contribution(
#'   ~ region, d, W, 0.3, 0.5, "mu"
#' )
#'
#' @noRd
.sim_unitreg_spatial_contribution <- function(formula, data, adjacency, scale, rho,component) {

  if (is.null(formula)) {
    return(list(
      contribution = rep(0, nrow(data)),
      effects = NULL
    ))
  }

  vars <- all.vars(formula)

  if (length(vars) != 1L) {
    stop(
      "Spatial formula for component '", component,
      "' must identify exactly one region variable.",
      call. = FALSE
    )
  }

  region_var <- vars[1L]

  if (!region_var %in% names(data)) {
    stop(
      "Spatial variable '", region_var,
      "' is missing from the simulated data.",
      call. = FALSE
    )
  }

  if (anyNA(data[[region_var]])) {
    stop(
      "Spatial variable '", region_var,
      "' contains missing values.",
      call. = FALSE
    )
  }

  region <- factor(data[[region_var]])
  region_levels <- levels(region)
  W <- .sim_unitreg_prepare_adjacency(adjacency, region_levels)

  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop(
      "'spatial_scale_", component,
      "' must be a positive finite number.",
      call. = FALSE
    )
  }

  rho_upper <- 1 - 1e-8

  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) ||  rho < 0 || rho >= rho_upper) {
    stop(
      "'spatial_rho_", component,
      "' must be in [0, ", format(rho_upper, digits = 10), ").",
      call. = FALSE
    )
  }

  D <- diag(rowSums(W))
  Q <- D - rho * W

  R <- tryCatch(
    chol(Q),
    error = function(e) {
      stop(
        "The proper CAR precision matrix is not positive definite for component '",
        component, "'.",
        call. = FALSE
      )
    }
  )

  z <- stats::rnorm(length(region_levels))
  u <- as.numeric(scale * backsolve(R, z))
  names(u) <- region_levels

  list(
    contribution = u[as.integer(region)],
    effects = u,
    precision = Q
  )
}

#' Simulate Bounded Data for unitregTMB
#'
#' @description
#' Simulates independent, clustered, longitudinal, or areal bounded responses
#' from distributions supported by \code{unitregTMB}. Fixed effects, ordinary
#' Gaussian random effects, zero/one augmentation, and proper CAR spatial
#' effects may be specified independently in the \code{mu}, \code{phi},
#' \code{p0}, and \code{p1} predictors.
#'
#' @param n_id Positive integer giving the number of primary clusters or
#'   individuals.
#' @param n_rep Positive integer giving the number of observations per primary
#'   cluster.
#' @param family A \code{unitregTMBFamily} object.
#' @param formula_mu,formula_phi,formula_p0,formula_p1 One-sided formulas for
#'   the four distributional predictors. Ordinary random effects use
#'   lme4-style terms such as \code{(1 + time | id)}. \code{NULL},
#'   \code{~0}, or \code{~-1} disables \code{p0} or \code{p1}.
#' @param beta_mu,beta_phi,beta_p0,beta_p1 Fixed-effect coefficient vectors on
#'   the corresponding link scales.
#' @param Sigma Backward-compatible alias for \code{Sigma_mu}.
#' @param Sigma_mu,Sigma_phi,Sigma_p0,Sigma_p1 Random-effect covariance
#'   matrices. For multiple random-effects terms, supply a list with one
#'   covariance matrix per term.
#' @param spatial,spatial_phi,spatial_p0,spatial_p1 Optional one-sided formulas
#'   identifying the areal region for proper CAR effects.
#' @param adjacency Adjacency matrix shared by requested CAR effects.
#' @param spatial_scale_mu,spatial_scale_phi,spatial_scale_p0,spatial_scale_p1
#'   Positive CAR scale parameters.
#' @param spatial_rho_mu,spatial_rho_phi,spatial_rho_p0,spatial_rho_p1 Proper
#'   CAR dependence parameters. Values must be at least zero and below
#'   \code{1 - 1e-8}, matching the current TMB parameterization.
#' @param cov_config Named list defining generated covariates. Supported
#'   distributions are \code{"normal"}, \code{"binomial"}, \code{"uniform"},
#'   \code{"poisson"}, \code{"factor"}, and \code{"categorical"}.
#' @param seed Optional finite numeric random-number seed.
#'
#' @details
#' Continuous observations are generated only for observations allocated to the
#' continuous mixture component. No arbitrary clipping such as \code{1e-10} is
#' applied. Exact floating-point zero or one returned by a continuous random
#' generator is moved only to machine-representable interior values.
#'
#' The proper CAR field uses precision
#' \deqn{
#' Q = D - \rho W
#' }
#' and covariance \eqn{\mathrm{scale}^2 Q^{-1}}, matching the current
#' \code{unitregTMB} spatial likelihood.
#'
#' @return A data frame containing \code{id}, generated covariates, and response
#'   \code{Y}. Generated ordinary and spatial latent effects, linear predictors,
#'   and distributional parameters are stored as attributes.
#'
#' @examples
#' \donttest{
#' d <- sim_unitreg(
#'   n_id = 20,
#'   n_rep = 4,
#'   family = beta_fam(model_for = "mean"),
#'   formula_mu = ~ x + (1 | id),
#'   formula_phi = ~ 1,
#'   beta_mu = c(-0.5, 0.8),
#'   beta_phi = log(20),
#'   Sigma_mu = matrix(0.30^2),
#'   cov_config = list(x = list(dist = "normal")),
#'   seed = 123
#' )
#'
#' head(d)
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom stats as.formula make.link model.matrix rbinom rnorm rpois runif
#' @export
sim_unitreg <- function(
    n_id = 100,
    n_rep = 1,
    family,
    formula_mu = ~ 1,
    formula_phi = ~ 1,
    formula_p0 = ~ 0,
    formula_p1 = ~ 0,
    beta_mu,
    beta_phi,
    beta_p0 = NULL,
    beta_p1 = NULL,
    Sigma = NULL,
    Sigma_mu = Sigma,
    Sigma_phi = NULL,
    Sigma_p0 = NULL,
    Sigma_p1 = NULL,
    spatial = NULL,
    spatial_phi = NULL,
    spatial_p0 = NULL,
    spatial_p1 = NULL,
    adjacency = NULL,
    spatial_scale_mu = 0.30,
    spatial_scale_phi = 0.30,
    spatial_scale_p0 = 0.30,
    spatial_scale_p1 = 0.30,
    spatial_rho_mu = 0.50,
    spatial_rho_phi = 0.50,
    spatial_rho_p0 = 0.50,
    spatial_rho_p1 = 0.50,
    cov_config = list(),
    seed = NULL) {

  if (!is.numeric(n_id) || length(n_id) != 1L || !is.finite(n_id) ||
      n_id < 1 || n_id != as.integer(n_id) ||
      !is.numeric(n_rep) || length(n_rep) != 1L || !is.finite(n_rep) ||
      n_rep < 1 || n_rep != as.integer(n_rep)) {
    stop("'n_id' and 'n_rep' must be positive integers.", call. = FALSE)
  }

  n_id <- as.integer(n_id)
  n_rep <- as.integer(n_rep)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("'seed' must be NULL or one finite number.", call. = FALSE)
    }
    set.seed(seed)
  }

  if (!inherits(family, "unitregTMBFamily")) {
    stop("'family' must be a unitregTMBFamily object.", call. = FALSE)
  }

  if (!is.list(cov_config)) {
    stop("'cov_config' must be a list.", call. = FALSE)
  }

  formula_args <- list(
    formula_mu = formula_mu,
    formula_phi = formula_phi,
    formula_p0 = formula_p0,
    formula_p1 = formula_p1,
    spatial = spatial,
    spatial_phi = spatial_phi,
    spatial_p0 = spatial_p0,
    spatial_p1 = spatial_p1
  )

  required <- c("formula_mu", "formula_phi")

  for (nm in names(formula_args)) {
    f <- formula_args[[nm]]

    if (is.null(f)) {
      if (nm %in% required) {
        stop("'", nm, "' cannot be NULL.", call. = FALSE)
      }
      next
    }

    if (!inherits(f, "formula") || length(f) != 2L) {
      stop("'", nm, "' must be a one-sided formula or NULL.", call. = FALSE)
    }
  }

  n <- n_id * n_rep
  data_sim <- data.frame(
    id = factor(rep(seq_len(n_id), each = n_rep))
  )

  formulas <- list(
    mu = formula_mu,
    phi = formula_phi,
    p0 = formula_p0,
    p1 = formula_p1
  )

  spatial_formulas <- list(
    mu = spatial,
    phi = spatial_phi,
    p0 = spatial_p0,
    p1 = spatial_p1
  )

  fixed_forms <- lapply(formulas, function(f) {
    if (is.null(f)) return(~ 0)
    reformulas::nobars(f)
  })

  all_vars <- unique(unlist(lapply(
    c(formulas, spatial_formulas),
    function(f) if (is.null(f)) character() else all.vars(f)
  )))

  all_vars <- setdiff(all_vars, "id")

  spatial_vars <- unique(unlist(lapply(
    spatial_formulas,
    function(f) if (is.null(f)) character() else all.vars(f)
  )))

  for (v in all_vars) {
    if (v %in% names(data_sim)) next

    conf <- cov_config[[v]]

    if (v %in% spatial_vars && is.null(conf)) {
      stop(
        "Spatial variable '", v,
        "' must have a 'cov_config' entry. Use a factor/categorical ",
        "distribution with levels matching the adjacency matrix.",
        call. = FALSE
      )
    }

    data_sim[[v]] <- .sim_unitreg_generate_covariate(
      n = n,
      conf = conf,
      name = v
    )
  }

  betas <- list(
    mu = beta_mu,
    phi = beta_phi,
    p0 = beta_p0,
    p1 = beta_p1
  )

  sigmas <- list(
    mu = Sigma_mu,
    phi = Sigma_phi,
    p0 = Sigma_p0,
    p1 = Sigma_p1
  )

  spatial_scales <- list(
    mu = spatial_scale_mu,
    phi = spatial_scale_phi,
    p0 = spatial_scale_p0,
    p1 = spatial_scale_p1
  )

  spatial_rhos <- list(
    mu = spatial_rho_mu,
    phi = spatial_rho_phi,
    p0 = spatial_rho_p0,
    p1 = spatial_rho_p1
  )

  eta <- setNames(vector("list", 4L), .unitreg_components)
  re_store <- setNames(vector("list", 4L), .unitreg_components)
  sp_store <- setNames(vector("list", 4L), .unitreg_components)

  active <- c(
    mu = TRUE,
    phi = TRUE,
    p0 = !.sim_unitreg_is_zero_formula(formula_p0),
    p1 = !.sim_unitreg_is_zero_formula(formula_p1)
  )

  for (comp in .unitreg_components) {
    if (!active[[comp]]) {
      eta[[comp]] <- rep(0, n)
      next
    }

    b <- betas[[comp]]

    if (!is.numeric(b) || !length(b) || any(!is.finite(b))) {
      stop(
        "Finite fixed-effect coefficients are required for active component '",
        comp, "'.",
        call. = FALSE
      )
    }

    X <- stats::model.matrix(
      fixed_forms[[comp]],
      data_sim
    )

    if (ncol(X) != length(b)) {
      stop(
        "Length of beta_", comp, " (", length(b),
        ") does not match the ", comp,
        " model matrix (", ncol(X), " columns).",
        call. = FALSE
      )
    }

    eta_comp <- as.vector(X %*% b)

    re_part <- .sim_unitreg_random_contribution(
      formula = formulas[[comp]],
      data = data_sim,
      Sigma = sigmas[[comp]],
      component = comp
    )

    eta_comp <- eta_comp + re_part$contribution
    re_store[[comp]] <- re_part$effects

    sp_part <- .sim_unitreg_spatial_contribution(
      formula = spatial_formulas[[comp]],
      data = data_sim,
      adjacency = adjacency,
      scale = spatial_scales[[comp]],
      rho = spatial_rhos[[comp]],
      component = comp
    )

    eta_comp <- eta_comp + sp_part$contribution
    sp_store[[comp]] <- sp_part$effects
    eta[[comp]] <- eta_comp
  }

  mu_link <- stats::make.link(family$link_r_name)
  phi_link <- stats::make.link(family$phi_link_r_name)

  mu <- .unitreg_open_unit_guard(
    mu_link$linkinv(eta$mu),
    "simulated mu"
  )

  phi <- as.numeric(phi_link$linkinv(eta$phi))

  if (any(!is.finite(phi)) || any(phi <= 0)) {
    stop("Simulated phi values must be positive and finite.", call. = FALSE)
  }

  probs <- .unitreg_probabilities(
    eta_p0 = eta$p0,
    eta_p1 = eta$p1,
    has_p0 = active[["p0"]],
    has_p1 = active[["p1"]]
  )

  prob_chk <- .unitreg_validate_aug_probs(probs$p0, probs$p1)
  p0 <- prob_chk$p0
  p1 <- prob_chk$p1
  pc <- prob_chk$pc

  u <- stats::runif(n)
  is0 <- u < p0
  is1 <- !is0 & u < (p0 + p1)
  isc <- !(is0 | is1)

  Y <- numeric(n)
  Y[is0] <- 0
  Y[is1] <- 1

  if (any(isc)) {
    y_cont <- get_random_continuous(
      family_code = family$family_code,
      n = sum(isc),
      mu = mu[isc],
      phi = phi[isc],
      tau = family$tau
    )

    if (length(y_cont) != sum(isc)) {
      stop(
        "The continuous random generator returned an incompatible length.",
        call. = FALSE
      )
    }

    Y[isc] <- .unitreg_open_unit_guard(
      y_cont,
      "simulated continuous response"
    )
  }

  data_sim$Y <- Y

  attr(data_sim, "random_effects") <- re_store
  attr(data_sim, "spatial_effects") <- sp_store
  attr(data_sim, "linear_predictors") <- eta
  attr(data_sim, "parameters") <- list(
    mu = mu,
    phi = phi,
    p0 = p0,
    p1 = p1,
    pc = pc
  )

  data_sim
}

