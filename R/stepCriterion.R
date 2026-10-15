#' Backward Fixed-Effect Selection for unitregTMB Models
#'
#' @description
#' Performs backward selection of fixed-effect terms in one distributional
#' predictor of a fitted \code{unitregTMB} model. Ordinary random-effect terms,
#' spatial effects, and all non-selected distributional predictors are
#' preserved.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param criterion Selection criterion: \code{"p-value"}, \code{"AIC"}, or
#'   \code{"BIC"}.
#' @param p.value.threshold Significance threshold used when
#'   \code{criterion = "p-value"}.
#' @param component Distributional predictor in which fixed effects are
#'   selected: \code{"mu"}, \code{"phi"}, \code{"p0"}, or \code{"p1"}.
#' @param trace Logical; print the backward-selection path.
#' @param ... Additional arguments passed to \code{update.unitregTMB} when
#'   candidate models are refitted.
#'
#' @details
#' Candidate fixed-effect terms are restricted using
#' \code{\link[stats]{drop.scope}}, so marginality is respected when
#' interactions are present. Random-effect and spatial terms are not candidates
#' for removal.
#'
#' With \code{criterion = "p-value"}, one-coefficient terms are evaluated by
#' Wald tests and multi-coefficient terms by joint Wald chi-squared tests using
#' the fixed-effect covariance matrix. A multi-coefficient term is not assigned
#' an approximate joint p-value when its covariance block is unavailable or
#' singular.
#'
#' With AIC or BIC, a candidate is eligible only when the refit converges and
#' uses the same number of observations as the current model. Comparisons across
#' different analysis samples are not performed.
#'
#' This function implements automated exploratory model selection. The reported
#' post-selection coefficient tests and standard errors do not account for the
#' selection process itself.
#'
#' @return The final fitted \code{unitregTMB} model.
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
#' fit_aic <- stepCriterion(
#'   fit,
#'   component = "mu",
#'   criterion = "AIC",
#'   trace = FALSE
#' )
#'
#' formula(fit_aic)
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom regcore stepCriterion
#' @importFrom stats AIC BIC drop.scope model.matrix pchisq terms update
#' @importFrom stats update.formula vcov
#' @importFrom reformulas nobars
#' @method stepCriterion unitregTMB
#' @rdname stepCriterion_unitregTMB
#' @export
stepCriterion.unitregTMB <- function(object, criterion = c("p-value", "AIC", "BIC"), p.value.threshold = 0.05,
                                     component = c("mu", "phi", "p0", "p1"), trace = TRUE, ...) {

  if (!inherits(object, "unitregTMB")) {
    stop("'object' must inherit from 'unitregTMB'.", call. = FALSE)
  }

  criterion <- match.arg(criterion)
  component <- match.arg(component)

  if (length(trace) != 1L || !is.logical(trace) || is.na(trace)) {
    stop("'trace' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!.unitreg_active_component(object, component)) {
    stop(
      "Component '", component, "' is not active in this model.",
      call. = FALSE
    )
  }

  if (!is.numeric(p.value.threshold) ||
      length(p.value.threshold) != 1L ||
      !is.finite(p.value.threshold) ||
      p.value.threshold <= 0 ||
      p.value.threshold >= 1) {
    stop(
      "'p.value.threshold' must be a single number in (0, 1).",
      call. = FALSE
    )
  }

  dots <- list(...)

  get_formula <- function(model) {
    formula(model, component = component)
  }

  refit_with_formula <- function(model, new_formula) {
    args <- c(list(object = model), dots)

    if (component == "mu") {
      args$formula. <- new_formula
    } else {
      args[[paste0(component, ".formula")]] <- new_formula
    }

    do.call(stats::update, args)
  }

  drop_term_formula <- function(formula, term) {
    rhs <- if (component == "mu") {
      paste(". ~ . -", term)
    } else {
      paste("~ . -", term)
    }

    stats::update.formula(formula, rhs)
  }

  removable_terms <- function(full_formula) {
    fixed_formula <- reformulas::nobars(full_formula)
    sc <- stats::drop.scope(stats::terms(fixed_formula))
    as.character(sc)
  }

  valid_refit <- function(model, reference) {
    if (!inherits(model, "unitregTMB")) return(FALSE)
    if (!isTRUE(as.integer(model$convergence)[1L] == 0L)) return(FALSE)
    if (!identical(as.integer(model$nobs), as.integer(reference$nobs))) {
      return(FALSE)
    }

    ll <- tryCatch(as.numeric(stats::logLik(model)), error = function(e) NA_real_)
    is.finite(ll)
  }

  joint_term_pvalues <- function(model, full_formula, candidates) {
    if (!length(candidates)) {
      return(setNames(numeric(), character()))
    }

    fixed_formula <- reformulas::nobars(full_formula)
    trm <- stats::terms(fixed_formula)
    term_labels <- attr(trm, "term.labels")

    X <- stats::model.matrix(
      fixed_formula,
      data = model$model
    )

    assign <- attr(X, "assign")
    b <- model$model_coef[[component]]
    se <- model$model_std[[paste0("std.", component)]]

    if (is.null(se)) {
      se <- rep(NA_real_, length(b))
    }

    V_all <- tryCatch(
      stats::vcov(model),
      error = function(e) NULL
    )

    V_comp <- NULL

    if (!is.null(V_all)) {
      nm <- paste0(component, ".", names(b))

      if (all(nm %in% rownames(V_all))) {
        V_comp <- V_all[nm, nm, drop = FALSE]
        dimnames(V_comp) <- list(names(b), names(b))
      }
    }

    out <- setNames(rep(NA_real_, length(candidates)), candidates)

    for (term in candidates) {
      j <- match(term, term_labels)

      if (is.na(j)) next

      cols <- colnames(X)[assign == j]
      idx <- match(cols, names(b))
      idx <- idx[!is.na(idx)]

      if (!length(idx)) next

      bb <- as.numeric(b[idx])
      df <- length(bb)
      stat <- NA_real_

      if (df == 1L) {
        ss <- as.numeric(se[idx])

        if (length(ss) == 1L && is.finite(ss) && ss > 0) {
          stat <- (bb / ss)^2
        }
      } else if (!is.null(V_comp)) {
        Vb <- V_comp[names(b)[idx], names(b)[idx], drop = FALSE]

        stat <- tryCatch(
          {
            R <- chol(Vb)
            z <- backsolve(R, bb, transpose = TRUE)
            sum(z^2)
          },
          error = function(e) NA_real_
        )
      }

      if (is.finite(stat) && df > 0L) {
        out[term] <- stats::pchisq(
          stat,
          df = df,
          lower.tail = FALSE
        )
      }
    }

    out
  }

  criterion_value <- function(model) {
    value <- if (criterion == "AIC") {
      as.numeric(stats::AIC(model))
    } else {
      as.numeric(stats::BIC(model))
    }

    if (!length(value) || !is.finite(value)) Inf else value
  }

  current_model <- object
  current_formula <- get_formula(current_model)
  step_no <- 1L

  if (trace) {
    cat(
      "\nBackward selection in component '",
      component,
      "' using ",
      criterion,
      "\n",
      sep = ""
    )
    cat(strrep("-", 60L), "\n", sep = "")
  }

  repeat {
    term_labels <- removable_terms(current_formula)

    if (!length(term_labels)) {
      if (trace) cat("No removable fixed-effect terms remain.\n")
      break
    }

    if (criterion == "p-value") {
      pvals <- joint_term_pvalues(
        current_model,
        current_formula,
        term_labels
      )

      finite <- is.finite(pvals)

      if (!all(finite)) {
        bad <- names(pvals)[!finite]
        warning(
          "Joint Wald p-values could not be computed for: ",
          paste(bad, collapse = ", "),
          ".",
          call. = FALSE
        )
      }

      if (!any(finite)) {
        warning(
          "No finite Wald p-values could be computed; selection stopped.",
          call. = FALSE
        )
        break
      }

      worst <- names(which.max(ifelse(finite, pvals, -Inf)))
      max_p <- pvals[worst]

      if (max_p <= p.value.threshold) {
        if (trace) {
          cat(
            "All testable remaining terms satisfy p <= ",
            format(p.value.threshold),
            ".\n",
            sep = ""
          )
        }
        break
      }

      new_formula <- drop_term_formula(
        current_formula,
        worst
      )

      if (trace) {
        cat(
          "Step ", step_no, ": removing '", worst,
          "' (joint Wald p = ", format(max_p, digits = 4), ")\n",
          sep = ""
        )
      }

      candidate <- tryCatch(
        refit_with_formula(current_model, new_formula),
        error = function(e) e
      )

      if (inherits(candidate, "error")) {
        warning(
          "Refitting failed after removing '", worst, "': ",
          conditionMessage(candidate),
          call. = FALSE
        )
        break
      }

      if (!identical(as.integer(candidate$nobs), as.integer(current_model$nobs))) {
        warning(
          "Removing '", worst,
          "' changed the analysis sample; selection stopped.",
          call. = FALSE
        )
        break
      }

      if (!valid_refit(candidate, current_model)) {
        warning(
          "The model obtained after removing '", worst,
          "' did not converge successfully; selection stopped.",
          call. = FALSE
        )
        break
      }

      current_model <- candidate
      current_formula <- get_formula(current_model)
      step_no <- step_no + 1L
    } else {
      current_value <- criterion_value(current_model)

      if (!is.finite(current_value)) {
        stop(
          "The current model has a non-finite ", criterion, ".",
          call. = FALSE
        )
      }

      candidates <- lapply(term_labels, function(term) {
        new_formula <- drop_term_formula(
          current_formula,
          term
        )

        fit <- tryCatch(
          refit_with_formula(current_model, new_formula),
          error = function(e) NULL
        )

        if (is.null(fit) || !valid_refit(fit, current_model)) {
          return(list(
            term = term,
            value = Inf,
            fit = NULL,
            formula = new_formula
          ))
        }

        list(
          term = term,
          value = criterion_value(fit),
          fit = fit,
          formula = new_formula
        )
      })

      vals <- vapply(candidates, `[[`, numeric(1), "value")
      best_idx <- which.min(vals)
      best <- candidates[[best_idx]]

      if (!is.finite(best$value) || best$value >= current_value) {
        if (trace) {
          cat("No removal improves ", criterion, ".\n", sep = "")
        }
        break
      }

      if (trace) {
        cat(
          "Step ", step_no, ": removing '", best$term,
          "' (", criterion, ": ",
          format(current_value, digits = 6), " -> ",
          format(best$value, digits = 6), ")\n",
          sep = ""
        )
      }

      current_model <- best$fit
      current_formula <- get_formula(current_model)
      step_no <- step_no + 1L
    }
  }

  if (trace) {
    cat(strrep("-", 60L), "\n", sep = "")
    cat("Final formula for component '", component, "':\n", sep = "")
    print(current_formula)
  }

  current_model
}
