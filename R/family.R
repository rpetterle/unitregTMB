#' @name unitregTMB-families
#' @title Family Objects for unitregTMB Models
#' @description 
#' Functions to specify the distribution family and link function for \code{unitregTMB} models.
#' 
#' @param model_for Character string indicating the parameterization to be used 
#'   (e.g., \code{"mean"}, \code{"mode"}, or \code{"quantile"}).
#' @param link Character string specifying the link function for the location parameter. 
#'   Supported links are \code{"logit"}, \code{"probit"}, \code{"cloglog"}, and \code{"cauchit"}.
#' @param tau Numeric value between 0 and 1 indicating the quantile to be modeled 
#'   (only applicable when \code{model_for = "quantile"}). Default is 0.5.
#' 
#' @return An object of class \code{unitregTMBFamily} containing the family details and link functions.
NULL

create_family_object <- function(name, tmb_id, link, tau = NULL, phi_link = "log") {
  link_id <- switch(link,
                    "logit"   = 0L,
                    "probit"  = 1L,
                    "cloglog" = 2L,
                    "cauchit" = 3L,
    stop("Link function '", link, "' is not supported. Use: logit, probit, cloglog, or cauchit.", call. = FALSE))

  valid_phi_links <- c("log", "logit")

  if (length(phi_link) != 1L || !phi_link %in% valid_phi_links) {
    stop("Unsupported phi link: ", phi_link, ".",call. = FALSE)
  }

  name <- sub(" \\(quantile\\)", "", name)

  structure(
    list(
      name = name,
      family_code = as.integer(tmb_id),
      link_r_name = link,
      link_code = as.integer(link_id),
      phi_link_r_name = phi_link,
      tau = tau
    ),
    class = "unitregTMBFamily"
  )
}


#' @export
print.unitregTMBFamily <- function(x, ...) {
  if (!is.null(x$tau)) {
    cat("\nFamily:", x$name, paste0("(tau = ", x$tau, ")"))
  } else {
    cat("\nFamily:", x$name)
  }
  cat("\nLink function:", x$link_r_name)
  cat("\n")
  invisible(x)
}

#' @rdname unitregTMB-families
#' @export
beta_fam <- function(model_for = c("mean", "mode"), link = "logit") {
  model_target <- match.arg(model_for)
  tmb_id <- switch(model_target, "mean" = 0, "mode" = 5)
  create_family_object(paste0("Beta (", model_target, ")"), tmb_id, link)
}

#' @rdname unitregTMB-families
#' @export
simplex <- function(model_for = c("mean"), link = "logit") {
  model_target <- match.arg(model_for)
  create_family_object("Simplex (mean)", 1, link)
}

#' @rdname unitregTMB-families
#' @export
vasicek <- function(model_for = c("mean", "quantile"), link = "logit", tau = 0.5) {
  model_target <- match.arg(model_for)
  current_tau <- NULL
  if (model_target == "mean") {
    tmb_id <- 2
  } else { 
    if (!is.numeric(tau) || tau <= 0 || tau >= 1) stop("tau must be in (0,1).", call. = FALSE)
    tmb_id <- 10
    current_tau <- tau
  }
  create_family_object(paste0("Vasicek (", model_target, ")"), tmb_id, link, current_tau, phi_link = "logit")
}

#' @rdname unitregTMB-families
#' @export
unitgamma <- function(model_for = c("mean", "mode"), link = "logit") {
  model_target <- match.arg(model_for)
  tmb_id <- switch(model_target, "mean" = 3, "mode" = 7)
  create_family_object(paste0("Unit-Gamma (", model_target, ")"), tmb_id, link)
}

#' @rdname unitregTMB-families
#' @export
bessel <- function(model_for = c("mean"), link = "logit") {
  model_target <- match.arg(model_for)
  create_family_object("Bessel (mean)", 4, link)
}

#' @rdname unitregTMB-families
#' @export
kumaraswamy <- function(model_for = c("mode", "quantile"), link = "logit", tau = 0.5) {
  model_target <- match.arg(model_for)
  current_tau <- NULL
  if (model_target == "mode") {
    tmb_id <- 6
  } else { 
    if (!is.numeric(tau) || tau <= 0 || tau >= 1) stop("tau must be in (0,1).", call. = FALSE)
    tmb_id <- 9
    current_tau <- tau
  }
  create_family_object(paste0("Kumaraswamy (", model_target, ")"), tmb_id, link, current_tau)
}

#' @rdname unitregTMB-families
#' @export
unitgompertz <- function(model_for = c("mode", "quantile"), link = "logit", tau = 0.5) {
  model_target <- match.arg(model_for)
  current_tau <- NULL
  if (model_target == "mode") {
    tmb_id <- 8
  } else { 
    if (!is.numeric(tau) || tau <= 0 || tau >= 1) stop("tau must be in (0,1).", call. = FALSE)
    tmb_id <- 12
    current_tau <- tau
  }
  create_family_object(paste0("Unit-Gompertz (", model_target, ")"), tmb_id, link, current_tau)
}

#' @rdname unitregTMB-families
#' @export
unitweibull <- function(model_for = c("quantile"), link = "logit", tau = 0.5) {
  model_target <- match.arg(model_for)
  if (!is.numeric(tau) || tau <= 0 || tau >= 1) stop("tau must be in (0,1).", call. = FALSE)
  create_family_object("Unit-Weibull (quantile)", 11, link, tau)
}

#' @rdname unitregTMB-families
#' @export
johnsonsb <- function(model_for = c("quantile"), link = "logit", tau = 0.5) {
  model_target <- match.arg(model_for)
  if (!is.numeric(tau) || tau <= 0 || tau >= 1) stop("tau must be in (0,1).", call. = FALSE)
  create_family_object("Johnson SB (quantile)", 13, link, tau)
}

#' @rdname unitregTMB-families
#' @export
ashw <- function(model_for = c("quantile"), link = "logit", tau = 0.5) {
  model_target <- match.arg(model_for)
  if (!is.numeric(tau) || tau <= 0 || tau >= 1) stop("tau must be in (0,1).", call. = FALSE)
  create_family_object("ASHW (quantile)", 14, link, tau)
}

#' @rdname unitregTMB-families
#' @export
ubs <- function(model_for = c("quantile"), link = "logit", tau = 0.5) {
  model_target <- match.arg(model_for)
  if (!is.numeric(tau) || tau <= 0 || tau >= 1) stop("tau must be in (0,1).", call. = FALSE)
  create_family_object("Unit-Birnbaum-Saunders (quantile)", 15, link, tau)
}
