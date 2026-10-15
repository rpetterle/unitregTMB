#' Create a unitregTMB Family Object
#'
#' @description
#' Internal constructor used by the exported family functions in
#' \code{unitregTMB}. It validates the location link, the dispersion/shape link,
#' the TMB family code, the modeled location parameter, and the quantile level
#' when applicable.
#'
#' @param name Character string giving the display name of the family.
#' @param tmb_id Integer family code used by the TMB objective function.
#' @param model_for Character string specifying whether the modeled location
#'   parameter is the \code{"mean"}, \code{"mode"}, or \code{"quantile"}.
#' @param link Character string specifying the link for the location parameter.
#' @param tau Optional quantile level in \eqn{(0,1)}.
#' @param phi_link Character string specifying the link for the dispersion or
#'   shape parameter. Supported values are \code{"log"} and \code{"logit"}.
#'
#' @return An object of class \code{unitregTMBFamily}.
#'
#' @examples
#' unitregTMB:::create_family_object(
#'   name = "Beta (mean)",
#'   tmb_id = 0L,
#'   model_for = "mean",
#'   link = "logit"
#' )
#'
#' @noRd
create_family_object <- function(name, tmb_id, model_for, link, tau = NULL, phi_link = "log") {
  valid_links <- c("logit", "probit", "cloglog", "cauchit")
  if (!is.character(link) || length(link) != 1L || is.na(link) || !link %in% valid_links) {
    stop("Unsupported location link. Use: logit, probit, cloglog, or cauchit.", call. = FALSE)
  }

  valid_model_for <- c("mean", "mode", "quantile")
  if (!is.character(model_for) || length(model_for) != 1L || is.na(model_for) || !model_for %in% valid_model_for) {
    stop("'model_for' must be one of: mean, mode, quantile.", call. = FALSE)
  }

  valid_phi_links <- c("log", "logit")
  if (!is.character(phi_link) || length(phi_link) != 1L || is.na(phi_link) || !phi_link %in% valid_phi_links) {
    stop("Unsupported phi link. Use: log or logit.", call. = FALSE)
  }

  if (!is.numeric(tmb_id) || length(tmb_id) != 1L || !is.finite(tmb_id) || tmb_id != as.integer(tmb_id)) {
    stop("'tmb_id' must be a single finite integer.", call. = FALSE)
  }

  if (model_for == "quantile") {
    if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau >= 1) {
      stop("'tau' must be a single finite number in (0, 1) for quantile models.", call. = FALSE)
    }
  } else if (!is.null(tau)) {
    stop("'tau' must be NULL unless model_for = 'quantile'.", call. = FALSE)
  }

  link_id <- switch(link, logit = 0L, probit = 1L, cloglog = 2L, cauchit = 3L)
  name <- sub(" \\(quantile\\)", "", name)

  structure(list(name = name, 
                 family_code = as.integer(tmb_id),
                 model_for = model_for,
                 link_r_name = link,
                 link_code = as.integer(link_id),
                 phi_link_r_name = phi_link,
                 tau = tau),
                 class = "unitregTMBFamily")
}

#' Print a unitregTMB Family Object
#'
#' @description
#' Prints the distribution name, modeled location parameter, location link,
#' dispersion/shape link, and quantile level when applicable.
#'
#' @param x An object of class \code{unitregTMBFamily}.
#' @param ... Additional arguments, currently unused.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' fam <- vasicek(model_for = "quantile", tau = 0.5)
#' print(fam)
#' @noRd
print.unitregTMBFamily <- function(x, ...) {
  cat("\nFamily:", x$name)
  cat("\nParameterization:", x$model_for)
  if (!is.null(x$tau)) cat("\nQuantile level:", x$tau)
  cat("\nLocation link:", x$link_r_name)
  cat("\nPhi link:", x$phi_link_r_name)
  cat("\n")
  invisible(x)
}

#' Family Objects for unitregTMB Models
#'
#' @description
#' Constructors for the continuous distributions implemented in
#' \code{unitregTMB}. Depending on the family, the location parameter can
#' represent the mean, mode, or a selected quantile of the continuous component.
#'
#' @param model_for Character string specifying the modeled location parameter.
#'   Available choices depend on the family and include \code{"mean"},
#'   \code{"mode"}, and \code{"quantile"}.
#' @param link Character string specifying the link for the modeled location
#'   parameter. Supported links are \code{"logit"}, \code{"probit"},
#'   \code{"cloglog"}, and \code{"cauchit"}.
#' @param tau Numeric quantile level in \eqn{(0,1)} used when
#'   \code{model_for = "quantile"}. The default is \code{0.5}.
#'
#' @details
#' The Beta family supports mean and mode parameterizations. Simplex and Bessel
#' support the mean parameterization. Vasicek supports mean and quantile
#' parameterizations and uses a logit link for its dispersion/shape parameter.
#' Unit-Gamma supports mean and mode parameterizations. Kumaraswamy and
#' Unit-Gompertz support mode and quantile parameterizations. Unit-Weibull,
#' Johnson SB, ASHW, and Unit-Birnbaum-Saunders are quantile-parameterized.
#'
#' @return An object of class \code{unitregTMBFamily} containing the family
#'   code, parameterization, location link, dispersion/shape link, and quantile
#'   level when applicable.
#'
#' @examples
#' beta_fam(model_for = "mean")
#' beta_fam(model_for = "mode", link = "probit")
#' simplex()
#' vasicek(model_for = "mean")
#' vasicek(model_for = "quantile", tau = 0.25)
#' unitgamma(model_for = "mean")
#' unitgamma(model_for = "mode")
#' bessel()
#' kumaraswamy(model_for = "mode")
#' kumaraswamy(model_for = "quantile", tau = 0.5)
#' unitgompertz(model_for = "mode")
#' unitgompertz(model_for = "quantile", tau = 0.75)
#' unitweibull(tau = 0.5)
#' johnsonsb(tau = 0.5)
#' ashw(tau = 0.5)
#' ubs(tau = 0.5)
#'
#' @name unitregTMB-families
NULL

#' @rdname unitregTMB-families
#' @export
beta_fam <- function(model_for = c("mean", "mode"), link = "logit") {
  model_target <- match.arg(model_for)
  tmb_id <- switch(model_target, mean = 0L, mode = 5L)
  create_family_object(paste0("Beta (", model_target, ")"), tmb_id, model_target, link)
}

#' @rdname unitregTMB-families
#' @export
simplex <- function(model_for = "mean", link = "logit") {
  model_target <- match.arg(model_for, "mean")
  create_family_object("Simplex (mean)", 1L, model_target, link)
}

#' @rdname unitregTMB-families
#' @export
vasicek <- function(model_for = c("mean", "quantile"), link = "logit", tau = 0.5) {
  model_target <- match.arg(model_for)
  if (model_target == "mean") {
    create_family_object("Vasicek (mean)", 2L, model_target, link, phi_link = "logit")
  } else {
    create_family_object("Vasicek (quantile)", 10L, model_target, link, tau = tau, phi_link = "logit")
  }
}

#' @rdname unitregTMB-families
#' @export
unitgamma <- function(model_for = c("mean", "mode"), link = "logit") {
  model_target <- match.arg(model_for)
  tmb_id <- switch(model_target, mean = 3L, mode = 7L)
  create_family_object(paste0("Unit-gamma (", model_target, ")"), tmb_id, model_target, link)
}

#' @rdname unitregTMB-families
#' @export
bessel <- function(model_for = "mean", link = "logit") {
  model_target <- match.arg(model_for, "mean")
  create_family_object("Bessel (mean)", 4L, model_target, link)
}

#' @rdname unitregTMB-families
#' @export
kumaraswamy <- function(model_for = c("mode", "quantile"), link = "logit", tau = 0.5) {
  model_target <- match.arg(model_for)
  if (model_target == "mode") {
    create_family_object("Kumaraswamy (mode)", 6L, model_target, link)
  } else {
    create_family_object("Kumaraswamy (quantile)", 9L, model_target, link, tau = tau)
  }
}

#' @rdname unitregTMB-families
#' @export
unitgompertz <- function(model_for = c("mode", "quantile"), link = "logit", tau = 0.5) {
  model_target <- match.arg(model_for)
  if (model_target == "mode") {
    create_family_object("Unit-Gompertz (mode)", 8L, model_target, link)
  } else {
    create_family_object("Unit-Gompertz (quantile)", 12L, model_target, link, tau = tau)
  }
}

#' @rdname unitregTMB-families
#' @export
unitweibull <- function(model_for = "quantile", link = "logit", tau = 0.5) {
  model_target <- match.arg(model_for, "quantile")
  create_family_object("Unit-Weibull (quantile)", 11L, model_target, link, tau = tau)
}

#' @rdname unitregTMB-families
#' @export
johnsonsb <- function(model_for = "quantile", link = "logit", tau = 0.5) {
  model_target <- match.arg(model_for, "quantile")
  create_family_object("Johnson SB (quantile)", 13L, model_target, link, tau = tau)
}

#' @rdname unitregTMB-families
#' @export
ashw <- function(model_for = "quantile", link = "logit", tau = 0.5) {
  model_target <- match.arg(model_for, "quantile")
  create_family_object("ASHW (quantile)", 14L, model_target, link, tau = tau)
}

#' @rdname unitregTMB-families
#' @export
ubs <- function(model_for = "quantile", link = "logit", tau = 0.5) {
  model_target <- match.arg(model_for, "quantile")
  create_family_object("Unit-Birnbaum-Saunders (quantile)", 15L, model_target, link, tau = tau)
}
