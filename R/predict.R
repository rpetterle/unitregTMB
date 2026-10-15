#' @title Predict from a Fitted unitregTMB Model
#'
#' @description
#' Computes conditional or population-level predictions from a fitted
#' \code{unitregTMB} model. Predictions are available for the \code{mu},
#' \code{phi}, \code{p0}, and \code{p1} distributional predictors and support
#' ordinary random effects and proper CAR spatial effects.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param newdata Optional data frame containing values at which predictions
#'   are required. If \code{NULL}, predictions are evaluated at the observations
#'   used to fit the model.
#' @param type Prediction type. \code{"response"} returns the augmented
#'   conditional mean for mean-parameterized families and the modeled
#'   response-scale location parameter for modal and quantile families.
#'   \code{"parameter"} returns the four distributional parameters.
#'   \code{"mu"}, \code{"phi"}, \code{"p0"}, and \code{"p1"} return a single
#'   parameter. \code{"link"} returns the linear predictor for
#'   \code{component}.
#' @param component Distributional component used when \code{type = "link"}:
#'   \code{"mu"}, \code{"phi"}, \code{"p0"}, or \code{"p1"}.
#' @param re.form \code{NULL} includes fitted ordinary and spatial latent
#'   effects. \code{NA} or \code{~0} excludes them and returns population-level
#'   predictions.
#' @param allow.new.levels Logical. If \code{TRUE}, unseen grouping levels or
#'   spatial regions in \code{newdata} receive zero latent-effect contribution.
#' @param offset Optional numeric offset for the \code{mu} predictor when
#'   \code{newdata} is supplied and the desired offset is not already specified
#'   through an \code{offset()} term in the model formula.
#' @param ... Additional arguments, currently unused.
#'
#' @details
#' For a mean-parameterized continuous distribution with zero/one augmentation,
#' \code{type = "response"} returns
#' \deqn{
#' p_1 + p_c\mu,\qquad p_c = 1-p_0-p_1.
#' }
#' For modal and quantile parameterizations, \code{mu} is the modeled mode or
#' quantile of the continuous component rather than its expectation.
#' Consequently, \code{type = "response"} returns this modeled response-scale
#' functional.
#'
#' With \code{re.form = NULL}, predictions at the original data use the fitted
#' conditional modes of ordinary random effects and proper CAR effects. With
#' \code{re.form = NA} or \code{~0}, all latent-effect contributions are set to
#' zero.
#'
#' @return A numeric vector for \code{type = "response"}, \code{"link"},
#'   \code{"mu"}, \code{"phi"}, \code{"p0"}, or \code{"p1"}. For
#'   \code{type = "parameter"}, a data frame with columns \code{mu},
#'   \code{phi}, \code{p0}, and \code{p1} is returned.
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
#' head(predict(fit, type = "mu"))
#' head(predict(fit, type = "parameter"))
#' head(predict(fit, type = "response", re.form = NA))
#' head(predict(fit, type = "link", component = "phi"))
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom stats as.formula delete.response make.link model.frame
#' @importFrom stats model.matrix model.offset na.pass plogis terms
#' @export
predict.unitregTMB <- function(object, newdata = NULL,  type = c("response", "link", "parameter", "mu", "phi", "p0", "p1"),
                               component = c("mu", "phi", "p0", "p1"), re.form = NULL, allow.new.levels = FALSE, offset = NULL, ...) {

  if (!inherits(object, "unitregTMB")) {
    stop("'object' must inherit from 'unitregTMB'.", call. = FALSE)
  }

  type <- match.arg(type)
  component <- match.arg(component)
  include_latent <- .unitreg_include_latent(re.form)

  if (length(allow.new.levels) != 1L ||
      !is.logical(allow.new.levels) ||
      is.na(allow.new.levels)) {
    stop("'allow.new.levels' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(newdata) && !is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  if (is.null(newdata) && !is.null(offset)) {
    stop("'offset' can only be supplied together with 'newdata'.", call. = FALSE)
  }

  if (type == "link" && !.unitreg_active_component(object, component)) {
    stop(
      "Component '", component, "' is not active in this model.",
      call. = FALSE
    )
  }

  if (is.null(newdata)) {
    eta <- setNames(vector("list", length(.unitreg_components)), .unitreg_components)

    for (comp in .unitreg_components) {
      eta[[comp]] <- .unitreg_fitted_predictor(
        object = object,
        component = comp,
        include_latent = include_latent
      )
    }

    params <- .unitreg_parameters_from_eta(
      object = object,
      eta_mu = eta$mu,
      eta_phi = eta$phi,
      eta_p0 = eta$p0,
      eta_p1 = eta$p1
    )
  } else {
    n_new <- nrow(newdata)
    eta <- setNames(vector("list", length(.unitreg_components)), .unitreg_components)

    for (comp in .unitreg_components) {
      if (!.unitreg_active_component(object, comp)) {
        eta[[comp]] <- rep(0, n_new)
        next
      }

      eta_comp <- .unitreg_fixed_predictor_newdata(
        object = object,
        component = comp,
        newdata = newdata,
        offset = if (comp == "mu") offset else NULL
      )

      if (include_latent) {
        eta_comp <- eta_comp +
          .unitreg_random_contribution_newdata(
            object = object,
            component = comp,
            newdata = newdata,
            allow.new.levels = allow.new.levels
          ) +
          .unitreg_spatial_contribution_newdata(
            object = object,
            component = comp,
            newdata = newdata,
            allow.new.levels = allow.new.levels
          )
      }

      eta[[comp]] <- eta_comp
    }

    params <- .unitreg_parameters_from_eta(
      object = object,
      eta_mu = eta$mu,
      eta_phi = eta$phi,
      eta_p0 = eta$p0,
      eta_p1 = eta$p1
    )
  }

  if (type == "link") {
    return(as.numeric(eta[[component]]))
  }

  if (type %in% .unitreg_components) {
    return(as.numeric(params[[type]]))
  }

  if (type == "parameter") {
    return(data.frame(
      mu = params$mu,
      phi = params$phi,
      p0 = params$p0,
      p1 = params$p1,
      check.names = FALSE
    ))
  }

  as.numeric(
    .unitreg_response_from_parameters(
      object = object,
      mu = params$mu,
      p0 = params$p0,
      p1 = params$p1,
      pc = params$pc
    )
  )
}

#' @keywords internal
#' @noRd
.unitreg_model_role <- function(object) {
  role <- NULL

  if (!is.null(object$family_object$model_for) &&
      length(object$family_object$model_for)) {
    role <- as.character(object$family_object$model_for)[1L]
  }

  if (!is.null(role) && role %in% c("mean", "mode", "quantile")) {
    return(role)
  }

  code <- NULL

  if (!is.null(object$family_object$family_code) &&
      length(object$family_object$family_code)) {
    code <- as.integer(object$family_object$family_code)[1L]
  } else if (!is.null(object$obj$env$data$family) &&
             length(object$obj$env$data$family)) {
    code <- as.integer(object$obj$env$data$family)[1L]
  }

  if (is.null(code) || is.na(code)) return("unknown")
  if (code >= 0L && code <= 4L) return("mean")
  if (code >= 5L && code <= 8L) return("mode")
  if (code >= 9L && code <= 15L) return("quantile")
  "unknown"
}

#' @keywords internal
#' @noRd
.unitreg_response_from_parameters <- function(object, mu, p0, p1, pc = NULL) {
  role <- .unitreg_model_role(object)

  if (identical(role, "mean")) {
    if (is.null(pc)) {
      pc <- 1 - p0 - p1
      tol <- 100 * .Machine$double.eps
      if (any(pc < -tol, na.rm = TRUE)) {
        stop("Invalid zero/one-inflation probabilities.", call. = FALSE)
      }
      pc[pc < 0 & pc >= -tol] <- 0
    }
    return(p1 + pc * mu)
  }

  mu
}

#' @keywords internal
#' @noRd
.unitreg_include_latent <- function(re.form) {
  if (is.null(re.form)) return(TRUE)

  if (length(re.form) == 1L &&
      is.atomic(re.form) &&
      is.na(re.form)) {
    return(FALSE)
  }

  if (inherits(re.form, "formula")) {
    txt <- gsub("[[:space:]]+", "", paste(deparse(re.form), collapse = ""))
    if (txt %in% c("~0", "~1-1", "~0+0")) return(FALSE)
  }

  stop(
    "Only re.form = NULL and re.form = NA or ~0 are currently supported.",
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
.unitreg_prediction_formula <- function(object, component) {
  switch(
    component,
    mu = object$formula,
    phi = object$phi.formula,
    p0 = object$p0.formula,
    p1 = object$p1.formula,
    stop("Unknown component: ", component, call. = FALSE)
  )
}

#' @keywords internal
#' @noRd
.unitreg_xlev_for_formula <- function(object, formula) {
  if (is.null(formula)) return(list())

  vars <- all.vars(formula)
  src <- object$model

  if (is.null(src) || !is.data.frame(src)) {
    src <- object$data
  }

  if (is.null(src) || !is.data.frame(src)) return(list())

  vars <- intersect(vars, names(src))
  out <- lapply(src[vars], function(x) if (is.factor(x)) levels(x) else NULL)
  out[!vapply(out, is.null, logical(1))]
}

#' @keywords internal
#' @noRd
.unitreg_fixed_predictor_newdata <- function(
    object,
    component,
    newdata,
    offset = NULL) {

  form <- .unitreg_prediction_formula(object, component)

  if (is.null(form)) {
    return(rep(0, nrow(newdata)))
  }

  fixed_form <- reformulas::nobars(form)
  tt <- stats::terms(fixed_form)

  if (component == "mu") {
    tt <- stats::delete.response(tt)
  }

  mf <- stats::model.frame(
    tt,
    data = newdata,
    na.action = stats::na.pass,
    xlev = .unitreg_xlev_for_formula(object, fixed_form)
  )

  X <- stats::model.matrix(tt, mf)
  b <- object$model_coef[[component]]

  if (is.null(b)) {
    return(rep(0, nrow(newdata)))
  }

  missing_cols <- setdiff(names(b), colnames(X))

  if (length(missing_cols)) {
    stop(
      "The model matrix for component '", component,
      "' is missing coefficient columns: ",
      paste(missing_cols, collapse = ", "),
      ". Check factor levels and transformations in 'newdata'.",
      call. = FALSE
    )
  }

  X <- X[, names(b), drop = FALSE]
  eta <- as.vector(X %*% b)

  if (component == "mu") {
    off_formula <- stats::model.offset(mf)

    if (!is.null(off_formula)) {
      eta <- eta + as.numeric(off_formula)
    }

    if (!is.null(offset)) {
      if (!is.numeric(offset)) {
        stop("'offset' must be numeric.", call. = FALSE)
      }

      if (length(offset) == 1L) {
        offset <- rep(offset, nrow(newdata))
      }

      if (length(offset) != nrow(newdata)) {
        stop(
          "'offset' must have length 1 or nrow(newdata).",
          call. = FALSE
        )
      }

      eta <- eta + as.numeric(offset)
    }
  }

  eta
}

#' @keywords internal
#' @noRd
.unitreg_fitted_predictor <- function(object, component, include_latent) {
  n <- object$nobs

  if (!.unitreg_active_component(object, component)) {
    return(rep(0, n))
  }

  dat <- object$obj$env$data
  X <- dat[[paste0("X_", component)]]
  b <- object$model_coef[[component]]

  if (is.null(X) || is.null(b) || nrow(X) != n || ncol(X) != length(b)) {
    stop(
      "Fixed-effect information is inconsistent for component '",
      component, "'.",
      call. = FALSE
    )
  }

  eta <- as.vector(X %*% b)

  if (component == "mu") {
    off <- dat$offset_mu
    if (is.null(off)) off <- rep(0, n)
    if (length(off) != n) {
      stop("Stored mu offset has an incompatible length.", call. = FALSE)
    }
    eta <- eta + as.numeric(off)
  }

  if (include_latent && .unitreg_has_re(object, component)) {
    Z <- dat[[paste0("Z_", component)]]
    u <- object$random_effects_modes[[component]]

    if (is.null(u) || !length(u)) {
      par_full <- .unitreg_full_par(object)
      u <- if (is.null(par_full)) NULL else par_full[[paste0("re_coef_", component)]]
    }

    if (is.null(Z) || is.null(u) || ncol(Z) != length(u)) {
      stop(
        "Random-effect information is inconsistent for component '",
        component, "'.",
        call. = FALSE
      )
    }

    eta <- eta + as.vector(Z %*% as.numeric(u))
  }

  if (include_latent && .unitreg_has_spatial(object, component)) {
    Zs <- dat[[paste0("Z_spatial_", component)]]
    eff <- object[[paste0("spatial_effects_", component)]]
    us <- NULL

    if (is.data.frame(eff) && "Estimate" %in% names(eff)) {
      us <- as.numeric(eff$Estimate)
    }

    if (is.null(us) || !length(us)) {
      par_full <- .unitreg_full_par(object)
      us <- if (is.null(par_full)) NULL else par_full[[paste0("u_spatial_", component)]]
    }

    if (is.null(Zs) || is.null(us) || ncol(Zs) != length(us)) {
      stop(
        "Spatial-effect information is inconsistent for component '",
        component, "'.",
        call. = FALSE
      )
    }

    eta <- eta + as.vector(Zs %*% as.numeric(us))
  }

  eta
}

#' @keywords internal
#' @noRd
.unitreg_random_contribution_newdata <- function(
    object,
    component,
    newdata,
    allow.new.levels = FALSE) {

  if (!.unitreg_has_re(object, component)) {
    return(rep(0, nrow(newdata)))
  }

  form <- .unitreg_prediction_formula(object, component)
  bars <- reformulas::findbars(form)
  info <- .unitreg_re_info(object, component)

  if (!length(bars) ||
      is.null(info) ||
      !isTRUE(info$has_random_effects)) {
    stop(
      "Random-effect metadata are inconsistent for component '",
      component, "'.",
      call. = FALSE
    )
  }

  re <- random_effects(object, component = component)

  if (is.null(re) || length(re) != length(bars)) {
    stop(
      "Conditional random effects could not be reconstructed for component '",
      component, "'.",
      call. = FALSE
    )
  }

  ans <- rep(0, nrow(newdata))

  for (i in seq_along(bars)) {
    bar <- bars[[i]]
    grp_expr <- bar[[3L]]

    grp <- tryCatch(
      eval(grp_expr, envir = newdata, enclos = environment(form)),
      error = function(e) {
        stop(
          "Could not evaluate grouping expression '",
          paste(deparse(grp_expr), collapse = " "),
          "' in newdata for component '", component, "'.",
          call. = FALSE
        )
      }
    )

    grp_chr <- as.character(grp)
    known_levels <- as.character(info$group_levels[[i]])
    idx <- match(grp_chr, known_levels)
    unknown <- is.na(idx)

    if (any(unknown) && !isTRUE(allow.new.levels)) {
      bad <- unique(grp_chr[unknown])
      stop(
        "New grouping level(s) in component '", component, "': ",
        paste(bad, collapse = ", "),
        ". Use allow.new.levels = TRUE for zero contribution from unseen levels.",
        call. = FALSE
      )
    }

    lhs_expr <- bar[[2L]]
    lhs_formula <- stats::as.formula(
      paste("~", paste(deparse(lhs_expr), collapse = " ")),
      env = environment(form)
    )

    mf_re <- stats::model.frame(
      lhs_formula,
      data = newdata,
      na.action = stats::na.pass
    )

    X_re <- stats::model.matrix(lhs_formula, mf_re)
    B <- as.matrix(re[[i]])

    if (is.null(colnames(B))) {
      stop(
        "Random-effect coefficient names are unavailable for component '",
        component, "'.",
        call. = FALSE
      )
    }

    missing_re <- setdiff(colnames(B), colnames(X_re))

    if (length(missing_re)) {
      stop(
        "Could not match random-coefficient columns in newdata for component '",
        component, "': ",
        paste(missing_re, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    X_re <- X_re[, colnames(B), drop = FALSE]
    contrib <- rep(0, nrow(newdata))
    known <- !unknown

    if (any(known)) {
      Bi <- B[idx[known], , drop = FALSE]
      contrib[known] <- rowSums(X_re[known, , drop = FALSE] * Bi)
    }

    ans <- ans + contrib
  }

  ans
}

#' @keywords internal
#' @noRd
.unitreg_spatial_contribution_newdata <- function(
    object,
    component,
    newdata,
    allow.new.levels = FALSE) {

  if (!.unitreg_has_spatial(object, component)) {
    return(rep(0, nrow(newdata)))
  }

  info <- object[[paste0("spatial_info_", component)]]

  if (is.null(info) ||
      is.null(info$region_var) ||
      is.null(info$region_levels)) {
    stop(
      "Spatial metadata are incomplete for component '",
      component, "'.",
      call. = FALSE
    )
  }

  region_var <- info$region_var

  if (!region_var %in% names(newdata)) {
    stop(
      "Spatial variable '", region_var,
      "' is missing from newdata for component '", component, "'.",
      call. = FALSE
    )
  }

  region <- as.character(newdata[[region_var]])
  levels_sp <- as.character(info$region_levels)
  idx <- match(region, levels_sp)
  unknown <- is.na(idx)

  if (any(unknown) && !isTRUE(allow.new.levels)) {
    bad <- unique(region[unknown])
    stop(
      "New spatial region(s) in component '", component, "': ",
      paste(bad, collapse = ", "),
      ". Use allow.new.levels = TRUE for zero contribution from unseen regions.",
      call. = FALSE
    )
  }

  eff <- object[[paste0("spatial_effects_", component)]]

  if (!is.null(eff) &&
      is.data.frame(eff) &&
      "Estimate" %in% names(eff)) {
    u <- as.numeric(eff$Estimate)
    names(u) <- rownames(eff)
    u <- u[levels_sp]
  } else {
    par_full <- .unitreg_full_par(object)
    u <- par_full[[paste0("u_spatial_", component)]]
  }

  if (is.null(u) || length(u) != length(levels_sp)) {
    stop(
      "Spatial conditional modes could not be recovered for component '",
      component, "'.",
      call. = FALSE
    )
  }

  ans <- rep(0, nrow(newdata))
  known <- !unknown
  ans[known] <- as.numeric(u[idx[known]])
  ans
}

#' @keywords internal
#' @noRd
.unitreg_probabilities <- function(eta_p0, eta_p1, has_p0, has_p1) {
  n <- max(length(eta_p0), length(eta_p1))

  if (!has_p0 && !has_p1) {
    return(list(
      p0 = rep(0, n),
      p1 = rep(0, n),
      pc = rep(1, n)
    ))
  }

  if (has_p0 && has_p1) {
    mm <- pmax(0, eta_p0, eta_p1)
    ebase <- exp(-mm)
    e0 <- exp(eta_p0 - mm)
    e1 <- exp(eta_p1 - mm)
    den <- ebase + e0 + e1

    return(list(
      p0 = e0 / den,
      p1 = e1 / den,
      pc = ebase / den
    ))
  }

  if (has_p0) {
    return(list(
      p0 = stats::plogis(eta_p0),
      p1 = rep(0, n),
      pc = stats::plogis(-eta_p0)
    ))
  }

  list(
    p0 = rep(0, n),
    p1 = stats::plogis(eta_p1),
    pc = stats::plogis(-eta_p1)
  )
}

#' @keywords internal
#' @noRd
.unitreg_parameters_from_eta <- function(
    object,
    eta_mu,
    eta_phi,
    eta_p0,
    eta_p1) {

  mu <- object$link.mu$linkinv(eta_mu)

  phi_link_name <- object$family_object$phi_link_r_name
  if (is.null(phi_link_name) || !length(phi_link_name)) {
    phi_link_name <- "log"
  }

  phi <- stats::make.link(as.character(phi_link_name)[1L])$linkinv(eta_phi)

  probs <- .unitreg_probabilities(
    eta_p0 = eta_p0,
    eta_p1 = eta_p1,
    has_p0 = isTRUE(object$has_p0),
    has_p1 = isTRUE(object$has_p1)
  )

  list(
    mu = as.numeric(mu),
    phi = as.numeric(phi),
    p0 = as.numeric(probs$p0),
    p1 = as.numeric(probs$p1),
    pc = as.numeric(probs$pc)
  )
}
