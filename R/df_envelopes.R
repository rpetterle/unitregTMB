#' Extract Simulated-Envelope Data from unitregTMB Models
#'
#' @description
#' Computes simulated residual envelopes for one or more fitted
#' \code{unitregTMB} models through the \code{regcore} envelope machinery and
#' combines the resulting coordinates in long format for subsequent plotting or
#' comparison.
#'
#' @param object A fitted object of class \code{unitregTMB}.
#' @param ... Additional fitted \code{unitregTMB} objects.
#' @param nsim Positive integer giving the number of bootstrap simulations used
#'   to construct each envelope.
#' @param halfnormal Logical. If \code{TRUE}, a half-normal envelope is
#'   constructed.
#' @param resid.type Character string specifying the residual type. Available
#'   choices are \code{"quantile"} and \code{"cox-snell"}.
#' @param ncpus Positive integer giving the number of outer bootstrap workers.
#' @param show.progress Logical. If \code{TRUE}, progress information is
#'   displayed.
#'
#' @details
#' Model labels include the distribution name and parameterization, when
#' available, without repeating the parameterization. They can also include
#' the quantile level, zero/one augmentation, ordinary random-effects
#' structure, and spatial structure. The actual simulation and envelope
#' construction are delegated to \code{regcore::plot_envelope}.
#'
#' @return A \code{data.frame} containing the envelope coordinates returned by
#'   \code{regcore}, including \code{Theoretical}, \code{Observed},
#'   \code{Lower}, \code{Median}, \code{Upper}, and a factor \code{Model}
#'   identifying the fitted model.
#'
#' @examples
#' \donttest{
#' data("wqi", package = "unitregTMB")
#'
#' fit <- unitregTMB(y ~ location,
#'                   family = unitgamma(model_for = "mean"),
#'                   data = wqi)
#'
#' env <- regcore::df_envelopes(fit, nsim = 5,
#'                              ncpus = 1,
#'                              show.progress = FALSE)
#'
#' head(env)
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom regcore df_envelopes
#' @method df_envelopes unitregTMB
#' @export
df_envelopes.unitregTMB <- function(object, ..., nsim = 50, halfnormal = TRUE,
                                    resid.type = c("quantile", "cox-snell"),
                                    ncpus = 1, show.progress = TRUE) {
  
  resid.type <- match.arg(resid.type)
  
  if (length(nsim) != 1L || !is.finite(nsim) || nsim < 1 ||
      nsim != as.integer(nsim)) {
    stop("'nsim' must be a positive integer.", call. = FALSE)
  }
  
  if (length(ncpus) != 1L || !is.finite(ncpus) || ncpus < 1 ||
      ncpus != as.integer(ncpus)) {
    stop("'ncpus' must be a positive integer.", call. = FALSE)
  }
  
  if (length(halfnormal) != 1L || !is.logical(halfnormal) ||
      is.na(halfnormal)) {
    stop("'halfnormal' must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (length(show.progress) != 1L || !is.logical(show.progress) ||
      is.na(show.progress)) {
    stop("'show.progress' must be TRUE or FALSE.", call. = FALSE)
  }
  
  nsim <- as.integer(nsim)
  ncpus <- as.integer(ncpus)
  
  models <- list(object, ...)
  ok <- vapply(models, inherits, logical(1), what = "unitregTMB")
  
  if (!all(ok)) {
    stop("All supplied objects must inherit from 'unitregTMB'.", call. = FALSE)
  }
  
  model_label <- function(model) {
    fam <- if (!is.null(model$family) && length(model$family)) {
      as.character(model$family)[1L]
    } else {
      "unitregTMB"
    }
    
    fam_obj <- model$family_object
    model_for <- NULL
    
    if (!is.null(fam_obj$model_for) && length(fam_obj$model_for)) {
      model_for <- as.character(fam_obj$model_for)[1L]
    } else if (!is.null(fam_obj$family_code) && length(fam_obj$family_code)) {
      code <- as.integer(fam_obj$family_code)[1L]
      if (is.finite(code)) {
        if (code %in% 0:4) model_for <- "mean"
        if (code %in% 5:8) model_for <- "mode"
        if (code %in% 9:15) model_for <- "quantile"
      }
    }
    
    details <- character()
    
    role_in_family <- FALSE
    if (!is.null(model_for) && length(model_for) == 1L &&
        !is.na(model_for) && nzchar(model_for)) {
      role_pattern <- paste0("\\(\\s*", model_for, "\\s*\\)")
      role_in_family <- grepl(role_pattern, fam, ignore.case = TRUE, perl = TRUE)
      
      if (!role_in_family) {
        details <- c(details, model_for)
      }
    }
    
    if (!is.null(model$tau) && length(model$tau) == 1L &&
        is.finite(model$tau)) {
      details <- c(details, paste0("tau=", format(model$tau, trim = TRUE)))
    }
    
    aug <- c(if (isTRUE(model$has_p0)) "p0" else NULL,
             if (isTRUE(model$has_p1)) "p1" else NULL)
    if (length(aug)) {
      details <- c(details, paste0("inflation=", paste(aug, collapse = "+")))
    }
    
    has_re <- any(
      isTRUE(model$has_random_effects_mu),
      isTRUE(model$has_random_effects_phi),
      isTRUE(model$has_random_effects_p0),
      isTRUE(model$has_random_effects_p1)
    )
    if (has_re) details <- c(details, "mixed")
    
    has_sp <- any(
      isTRUE(model$has_spatial_effects_mu),
      isTRUE(model$has_spatial_effects_phi),
      isTRUE(model$has_spatial_effects_p0),
      isTRUE(model$has_spatial_effects_p1)
    )
    if (has_sp) details <- c(details, "spatial")
    
    if (length(details)) {
      paste0(fam, " [", paste(details, collapse = ", "), "]")
    } else {
      fam
    }
  }
  
  labels <- make.unique(vapply(models, model_label, character(1)))
  
  out <- lapply(seq_along(models), function(i) {
    if (show.progress) {
      cat("\nEnvelope for model: ", labels[i], "\n", sep = "")
    }
    
    env <- regcore::plot_envelope(
      models[[i]],
      plot = FALSE,
      nsim = nsim,
      halfnormal = halfnormal,
      resid.type = resid.type,
      ncpus = ncpus,
      show.progress = show.progress
    )
    
    if (is.null(env) || is.null(env$data)) {
      stop(
        "Envelope data could not be extracted for model '",
        labels[i],
        "'.",
        call. = FALSE
      )
    }
    
    required <- c("Theoretical", "Observed", "Lower", "Median", "Upper")
    missing <- setdiff(required, names(env$data))
    
    if (length(missing)) {
      stop(
        "Envelope data for model '",
        labels[i],
        "' are missing required columns: ",
        paste(missing, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    
    z <- env$data
    z$Model <- factor(labels[i], levels = labels)
    z
  })
  
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}
