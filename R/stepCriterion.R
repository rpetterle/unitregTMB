#' @name stepCriterion.unitregTMB
#'
#' @title Stepwise Variable Selection for unitregTMB Models
#' 
#' @description Performs backward variable selection for the location (mu) submodel 
#' based on p-values from the Wald test or Information Criteria (AIC, BIC). Random 
#' effects and parameters for other submodels (phi, p0, p1) are kept intact during 
#' the selection process.
#' 
#' @param object A fitted object of class \code{unitregTMB}.
#' @param criterion Character string specifying the selection criterion. 
#'        Options are \code{"p-value"} (default), \code{"AIC"}, or \code{"BIC"}.
#' @param p.value.threshold Numeric. The significance level threshold to keep a 
#'        variable in the model when \code{criterion = "p-value"}. Default is 0.05.
#' @param ... Additional arguments passed to methods.
#' 
#' @return A fitted \code{unitregTMB} model object representing the final selected model.
#' 
#' @examples
#' \donttest{
#' # Assuming 'fit' is a full unitregTMB model with multiple covariates:
#' # final_model <- stepCriterion(fit, criterion = "p-value", p.value.threshold = 0.05)
#' #
#' # # Or using the Akaike Information Criterion (AIC):
#' # final_model_aic <- stepCriterion(fit, criterion = "AIC")
#' }
#' 
#' @importFrom stats terms update.formula pchisq update
#' @importFrom reformulas nobars
#' @export
stepCriterion.unitregTMB <- function(object, criterion = c("p-value", "AIC", "BIC"), 
                                     p.value.threshold = 0.05, ...) {
  
  criterion <- match.arg(criterion)
  
  cat(sprintf("\nStarting Backward Elimination based on: %s\n", criterion))
  cat(rep("-", 50), "\n", sep = "")
  
  current_model <- object
  formula_mu <- current_model$formula
  
  keep_stepping <- TRUE
  step_count <- 1
  
  while(keep_stepping) {
    fixed_form <- reformulas::nobars(formula_mu)
    if (is.null(fixed_form)) fixed_form <- ~ 1
    
    terms_obj <- stats::terms(fixed_form)
    term_labels <- attr(terms_obj, "term.labels")
    
    if (length(term_labels) == 0) {
      cat("Only the intercept (and/or random effects) remains. Stopping.\n")
      break
    }
    
    if (criterion == "p-value") {
      est <- current_model$model_coef$mu
      se <- current_model$model_std$std.mu
      
      p_vals <- stats::pchisq((est / se)^2, df = 1, lower.tail = FALSE)
      p_vals <- p_vals[names(p_vals) != "(Intercept)"]
      
      if (length(p_vals) == 0) break
      
      max_p <- max(p_vals, na.rm = TRUE)
      worst_var <- names(p_vals)[which.max(p_vals)]
      
      if (max_p > p.value.threshold) {
        
        sorted_terms <- term_labels[order(nchar(term_labels), decreasing = TRUE)]
        worst_var_clean <- worst_var
        for (tl in sorted_terms) {
          if (startsWith(worst_var, tl) || grepl(tl, worst_var, fixed = TRUE)) {
            worst_var_clean <- tl
            break
          }
        }
        
        cat(sprintf("Step %d: Removing '%s' (p-value = %.4f)\n", step_count, worst_var_clean, max_p))
        
        new_formula <- stats::update.formula(formula_mu, paste(". ~ . -", worst_var_clean))
        
        if (length(attr(stats::terms(reformulas::nobars(new_formula)), "term.labels")) == length(term_labels)) {
          cat("Warning: Could not parse variable removal. Stopping to prevent infinite loop.\n")
          break
        }
        
        formula_mu <- new_formula
        current_model <- stats::update(current_model, formula. = formula_mu)
        step_count <- step_count + 1
        
      } else {
        cat(sprintf("All remaining variables are significant (p <= %.2f). Stopping.\n", p.value.threshold))
        keep_stepping <- FALSE
      }
      
    } else {
      current_crit <- if(criterion == "AIC") current_model$AIC else current_model$BIC
      
      best_new_crit <- Inf
      var_to_remove <- NULL
      best_new_model <- NULL
      
      for (var in term_labels) {
        test_formula <- stats::update.formula(formula_mu, paste(". ~ . -", var))
        
        test_model <- tryCatch(
          stats::update(current_model, formula. = test_formula, evaluate = TRUE), 
          error = function(e) NULL
        )
        
        if (!is.null(test_model)) {
          if (test_model$nobs != current_model$nobs) {
             warning(sprintf("Dropping '%s' changed the number of observations (due to NAs). Criteria comparison may be invalid.", var), call. = FALSE)
          }
          
          test_crit <- if(criterion == "AIC") test_model$AIC else test_model$BIC
          
          if (test_crit < best_new_crit) {
            best_new_crit <- test_crit
            var_to_remove <- var
            best_new_model <- test_model
          }
        }
      }
      
      if (best_new_crit < current_crit) {
        cat(sprintf("Step %d: Removing '%s' (Improved %s from %.2f to %.2f)\n", 
                    step_count, var_to_remove, criterion, current_crit, best_new_crit))
        
        new_formula <- stats::update.formula(formula_mu, paste(". ~ . -", var_to_remove))
        
        if (length(attr(stats::terms(reformulas::nobars(new_formula)), "term.labels")) == length(term_labels)) {
           cat("Warning: Could not parse variable removal. Stopping.\n")
           break
        }
        
        formula_mu <- new_formula
        current_model <- best_new_model
        step_count <- step_count + 1
        
      } else {
        cat(sprintf("No removal improves the %s. Stopping.\n", criterion))
        keep_stepping <- FALSE
      }
    }
  }
  
  cat(rep("-", 50), "\n", sep = "")
  cat("Final Model Formula:\n")
  print(current_model$formula)
  
  return(current_model)
}
