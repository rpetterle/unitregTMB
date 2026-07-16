#' @name summary_coef.unitregTMB
#' 
#' @title Summarize Coefficients from unitregTMB Models
#' 
#' @description 
#' Extracts, formats, and compiles fixed-effects coefficients and their corresponding 
#' standard errors from one or multiple fitted \code{unitregTMB} models. This function is 
#' particularly useful for generating side-by-side model comparisons.
#' 
#' @param object An object of class \code{unitregTMB}.
#' @param ... Additional \code{unitregTMB} model objects for comparison.
#' @param component A character string or vector specifying which model components to extract. 
#'   Options include \code{"all"}, \code{"mu"} (location), \code{"phi"} (precision/dispersion), 
#'   \code{"p0"} (zero-inflation), and \code{"p1"} (one-inflation). Default is \code{"mu"}.
#' @param digits An integer specifying the number of decimal places to be used in the output. Default is 3.
#' @param print_tex A logical value. If \code{TRUE}, the summary prepares attributes to be 
#'   printed as a LaTeX table. Default is \code{FALSE}.
#' 
#' @return An object of class \code{summary_coef.unitregTMB} containing a list of compiled 
#'   coefficients, standard errors, and metadata for each provided model.
#' 
#' @examples
#' \donttest{
#' # Assuming 'da' is your dataset:
#' # fit1 <- unitregTMB(Y ~ educ, data = da, family = vasicek())
#' # fit2 <- unitregTMB(Y ~ educ, data = da, family = kumaraswamy())
#' # 
#' # # Console Output Comparison
#' # summary_coef(fit1, fit2, component = "all", digits = 2)
#' # 
#' # # LaTeX Table Output Comparison
#' # summary_coef(fit1, fit2, component = "mu", print_tex = TRUE)
#' }
#' 
#' @rdname summary_coef_unitregTMB
#' @export
summary_coef.unitregTMB <- function(object, ..., component = "mu", digits = 3, print_tex = FALSE) {
  
  models <- list(object, ...)
  
  is_valid_model <- sapply(models, inherits, "unitregTMB")
  if (!all(is_valid_model)) {
    stop("All unnamed arguments must be 'unitregTMB' objects. Check if you misspelled an argument name.", call. = FALSE)
  }
  
  if (length(models) > 1) {
    num_coefs <- sapply(models, function(m) sum(sapply(m$model_coef, length)))
    if (length(unique(num_coefs)) > 1) {
      warning("Models have a different number of coefficients. Comparison may not be straightforward.", call. = FALSE)
    }
  }
  
  results <- lapply(models, function(model) {
    
    coefs_mu <- model$model_coef$mu
    coefs_phi <- model$model_coef$phi
    coefs_p0 <- model$model_coef$p0
    coefs_p1 <- model$model_coef$p1
    
    ses_mu <- model$model_std$std.mu
    ses_phi <- model$model_std$std.phi
    ses_p0 <- model$model_std$std.p0
    ses_p1 <- model$model_std$std.p1
    
    clean_and_prefix <- function(coef_vec, prefix) {
      if (is.null(coef_vec)) return(NULL)
      names(coef_vec) <- paste0(prefix, trimws(names(coef_vec)))
      return(coef_vec)
    }
    
    if (!is.null(coefs_mu)) names(coefs_mu) <- trimws(names(coefs_mu))
    coefs_phi <- clean_and_prefix(coefs_phi, "phi_")
    coefs_p0 <- clean_and_prefix(coefs_p0, "p0_")
    coefs_p1 <- clean_and_prefix(coefs_p1, "p1_")
    
    if (!is.null(ses_mu)) names(ses_mu) <- trimws(names(ses_mu))
    ses_phi <- clean_and_prefix(ses_phi, "phi_")
    ses_p0 <- clean_and_prefix(ses_p0, "p0_")
    ses_p1 <- clean_and_prefix(ses_p1, "p1_")
    
    all_coefs <- c(coefs_mu, coefs_phi, coefs_p0, coefs_p1)
    all_std_errors <- c(ses_mu, ses_phi, ses_p0, ses_p1)
    
    if (length(all_coefs) == 0) {
      return(list(model_name = model$family, coefficients = numeric(0), std_errors = numeric(0), coef_names = character(0)))
    }
    
    all_names <- names(all_coefs)
    all_std_errors <- all_std_errors[all_names]
    
    is_random_effect <- grepl("log_chol_re", all_names)
    fixed_indices <- !is_random_effect
    indices_to_keep <- logical(length(all_names))
    
    if ("all" %in% component) {
      indices_to_keep <- fixed_indices
    } else {
      if ("mu" %in% component) {
        indices_to_keep <- indices_to_keep | (!grepl("^(phi_|p0_|p1_)", all_names) & fixed_indices)
      }
      if ("phi" %in% component) {
        indices_to_keep <- indices_to_keep | (grepl("^phi_", all_names) & fixed_indices)
      }
      if ("p0" %in% component) {
        indices_to_keep <- indices_to_keep | (grepl("^p0_", all_names) & fixed_indices)
      }
      if ("p1" %in% component) {
        indices_to_keep <- indices_to_keep | (grepl("^p1_", all_names) & fixed_indices)
      }
    }
    
    list(
      model_name   = model$family,
      tau          = model$tau,
      coefficients = all_coefs[indices_to_keep],
      std_errors   = all_std_errors[indices_to_keep],
      coef_names   = all_names[indices_to_keep]
    )
  })
  
  class(results) <- "summary_coef.unitregTMB"
  
  attr(results, "digits") <- digits
  attr(results, "print_tex") <- print_tex 
  
  return(results)
}

#' @rdname summary_coef_unitregTMB
#' @param x An object of class \code{summary_coef.unitregTMB}.
#' @export
print.summary_coef.unitregTMB <- function(x, digits = NULL, print_tex = NULL, ...) {
  
  if (is.null(digits)) {
    digits <- attr(x, "digits")
    if (is.null(digits)) digits <- 3
  }
  
  if (is.null(print_tex)) {
    print_tex <- attr(x, "print_tex")
    if (is.null(print_tex)) print_tex <- FALSE
  }
  
  capitalize <- function(s) paste0(toupper(substring(s, 1, 1)), substring(s, 2))
  
  model_col_names <- sapply(x, function(m) {
    fam_name <- if (is.list(m$model_name)) m$model_name$name else m$model_name
    
    if (is.null(fam_name) || length(fam_name) == 0) return("Unknown")
    if (grepl("\\(", fam_name)) return(fam_name)
    
    family_parts <- strsplit(as.character(fam_name), "_")[[1]]
    dist_name <- capitalize(family_parts[1])
    parameter_str <- NULL
    
    if (!is.null(m$tau) && is.numeric(m$tau)) {
      parameter_str <- paste0("tau = ", m$tau)
    } else if (length(family_parts) > 1) {
      parameter_str <- family_parts[2]
    }
    
    if (!is.null(parameter_str)) {
      paste0(dist_name, " (", parameter_str, ")")
    } else {
      dist_name
    }
  })
  
  list_of_name_vectors <- lapply(x, `[[`, "coef_names")
  cleaned_list_of_vectors <- lapply(list_of_name_vectors, function(name_vec) {
    if (is.null(name_vec)) return(character(0))
    trimws(name_vec)
  })
  all_coef_names <- Reduce(union, cleaned_list_of_vectors)
  
  if (length(all_coef_names) == 0) {
    cat("No coefficients to display.\n")
    return(invisible(x))
  }
  
  results_matrix <- matrix("---",
                           nrow = length(all_coef_names),
                           ncol = length(x),
                           dimnames = list(all_coef_names, model_col_names))
  
  for (j in 1:length(x)) {
    model_data <- x[[j]]
    if (length(model_data$coefficients) == 0) next
    
    for (i in 1:length(model_data$coefficients)) {
      coef_name <- trimws(model_data$coef_names[i])
      coef_val <- model_data$coefficients[i]
      se_val <- model_data$std_errors[i]
      
      if (is.na(se_val)) {
        formatted_str <- sprintf(paste0("%.", digits, "f (---)"), coef_val)
      } else {
        formatted_str <- sprintf(paste0("%.", digits, "f (%.", digits, "f)"), coef_val, se_val)
      }
      results_matrix[coef_name, j] <- formatted_str
    }
  }
  
  if (print_tex) {
    cat("\n% --- LaTeX Table generated by unitregTMB ---\n")
    cat("\\begin{table}[ht]\n\\centering\n")
    
    align_str <- paste0("l", paste(rep("c", ncol(results_matrix)), collapse = ""))
    cat(sprintf("\\begin{tabular}{%s}\n\\hline\n", align_str))
    
    col_headers <- colnames(results_matrix)
    cat("\\textbf{Predictor} & ", paste(sprintf("\\textbf{%s}", col_headers), collapse = " & "), " \\\\ \\hline\n")
    
    row_names <- rownames(results_matrix)
    row_names <- gsub("_", "\\\\_", row_names) 
    
    for (i in 1:nrow(results_matrix)) {
      row_vals <- results_matrix[i, ]
      cat(sprintf("%s & %s \\\\ \n", row_names[i], paste(row_vals, collapse = " & ")))
    }
    
    cat("\\hline\n\\end{tabular}\n\\end{table}\n")
    
  } else {
    cat("Regression Coefficients (Standard Errors):\n")
    col_width_first <- max(nchar(rownames(results_matrix))) + 2
    col_width_other <- max(sapply(c(results_matrix, colnames(results_matrix)), nchar)) + 2
    
    cat(format("Predictor", width = col_width_first, justify = "left"))
    for (name in colnames(results_matrix)) {
      cat(format(name, width = col_width_other, justify = "right"))
    }
    cat("\n")
    
    for (i in 1:nrow(results_matrix)) {
      cat(format(rownames(results_matrix)[i], width = col_width_first, justify = "left"))
      for (j in 1:ncol(results_matrix)) {
        cat(format(results_matrix[i, j], width = col_width_other, justify = "right"))
      }
      cat("\n")
    }
  }
  
  invisible(x)
}
