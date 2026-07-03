#' Extract Model Equations in LaTeX Format for unitregTMB
#' 
#' @description Generates formatted LaTeX equations for fitted unitregTMB models, 
#' supporting both symbolic and numeric representation. It automatically identifies 
#' random effects and assigns proper indices to longitudinal/repeated measures.
#' 
#' @param model A fitted unitregTMB object.
#' @param mode Character; "symbolic" (default) for theoretical equations with greek letters, 
#'        or "numeric" for fitted equations with estimated coefficients.
#' @param component Character; specifies which submodel to extract: "mu", "phi", "p0", "p1", or "all".
#' @param link_style Character; "name" (default) to print the link function as text (e.g., probit, cauchit), 
#'        or "inverse" to print the mathematical inverse function (e.g., \\Phi^{-1}).
#' @param grouping_var Character; a substring matching the name of the repeated measures variable 
#'        (e.g., "tooth"). This variable will receive the compound index {ij}. Defaults to NULL.
#' @param digits Integer; number of decimal places for numeric mode.
#' @param custom_labels Named character vector; pairs of old variable names and desired new labels.
#' @param index_labels Logical; if TRUE (default), appends the domain bounds (e.g., for i = 1,...,N).
#' 
#' @return A character string (invisibly) containing the LaTeX code for the equations. 
#'   The code is also printed to the console.
#' 
#' @examples
#' \donttest{
#' # Assuming 'fit' is a fitted unitregTMB model:
#' # Extract the theoretical symbolic equation for the mean (mu)
#' # extract_equations(fit, mode = "symbolic", component = "mu")
#' 
#' # Extract the numeric fitted equation for all components with custom labels
#' # eq <- extract_equations(fit, mode = "numeric", component = "all", 
#' #                         custom_labels = c("educHS+" = "HighSchool"))
#' }
#' @export
extract_equations <- function(model, mode = c("symbolic", "numeric"), 
                              component = c("mu", "phi", "p0", "p1", "all"), 
                              link_style = c("name", "inverse"),
                              grouping_var = NULL, digits = 2, 
                              custom_labels = NULL, index_labels = TRUE) {
  
  mode <- match.arg(mode)
  component <- match.arg(component)
  link_style <- match.arg(link_style)
  
  sd_rep <- if (!is.null(model$sd_report)) model$sd_report else model$sdreport
  if (is.null(sd_rep) && mode == "numeric") {
    stop("Model must be fitted with standard errors (sdreport) to extract numeric equations.")
  }
  
  link_name <- "logit"
  if (!is.null(model$link.mu) && !is.null(model$link.mu$name)) {
    link_name <- model$link.mu$name[1]
  } else if (is.list(model$family)) {
    if (!is.null(model$family$link_r_name)) link_name <- model$family$link_r_name[1]
    else if (!is.null(model$family$link)) link_name <- model$family$link[1]
  }
  
  has_re <- isTRUE(model$has_random_effects_mu) || isTRUE(model$obj$env$data$has_random_effects_mu == 1)
  
  lhs_idx <- if (has_re && !is.null(grouping_var)) "{ij}" else "i"
  mu_symbol <- sprintf("\\mu_{%s}", lhs_idx)
  
  get_link_latex <- function(lnk, param_name) {
    if (is.null(lnk) || length(lnk) != 1) lnk <- "logit"
    
    if (link_style == "name") {
      if (lnk == "log") return(sprintf("\\log(%s)", param_name))
      return(sprintf("\\mathrm{%s}(%s)", lnk, param_name))
    } else {
      switch(lnk,
             "logit"   = sprintf("\\log\\left(\\frac{%s}{1 - %s}\\right)", param_name, param_name),
             "probit"  = sprintf("\\Phi^{-1}(%s)", param_name),
             "cloglog" = sprintf("\\log(-\\log(1-%s))", param_name),
             "cauchit" = sprintf("\\tan(\\pi(%s - 0.5))", param_name),
             "log"     = sprintf("\\log(%s)", param_name),
             sprintf("\\mathrm{%s}(%s)", lnk, param_name) # Fallback
      )
    }
  }
  
  get_idx <- function(v) {
    if (has_re && !is.null(grouping_var) && grepl(grouping_var, v)) return("{ij}")
    return("i")
  }
  
  chunk_terms <- function(terms_vec, chunk_size = 3, is_numeric = FALSE) {
    if (length(terms_vec) == 0) return("")
    chunks <- split(terms_vec, ceiling(seq_along(terms_vec) / chunk_size))
    res <- ""
    for (i in seq_along(chunks)) {
      if (is_numeric) {
        line <- paste(chunks[[i]], collapse = " ")
        if (i == 1) res <- line else res <- paste0(res, " \\nonumber \\\\ \n  & \\quad ", line)
      } else {
        line <- paste(chunks[[i]], collapse = " + ")
        if (i == 1) res <- line else res <- paste0(res, " + \\nonumber \\\\ \n  & \\quad + ", line)
      }
    }
    return(res)
  }

  comps_to_print <- if (component == "all") c("mu", "phi", "p0", "p1") else component

  lines_out <- c()
  
  if (mode == "symbolic") {
    
    if ("mu" %in% comps_to_print) {
      fix_form <- reformulas::nobars(model$formula)
      mt <- stats::terms(fix_form)
      term_labels <- attr(mt, "term.labels")
      
      re_names <- character(0)
      if (has_re && !is.null(model$re_info_mu$cnms)) {
        re_names <- model$re_info_mu$cnms[[1]]
      }
      
      rhs_mu <- c()
      b_idx <- 0
      
      if (attr(mt, "intercept") == 1) {
        if ("(Intercept)" %in% re_names) {
          rhs_mu <- c(rhs_mu, sprintf("(\\beta_{%d} + u_{%di})", b_idx, b_idx))
        } else {
          rhs_mu <- c(rhs_mu, sprintf("\\beta_{%d}", b_idx))
        }
        b_idx <- b_idx + 1
      }
      
      for (var in term_labels) {
        print_var <- var
        if (!is.null(custom_labels) && var %in% names(custom_labels)) print_var <- custom_labels[[var]]
        safe_var <- gsub("_", "\\\\_", print_var)
        
        v_idx <- get_idx(var)
        
        if (var %in% re_names) {
          rhs_mu <- c(rhs_mu, sprintf("(\\beta_{%d} + u_{%di})~\\texttt{%s}_{%s}", b_idx, b_idx, safe_var, v_idx))
        } else {
          rhs_mu <- c(rhs_mu, sprintf("\\beta_{%d}~\\texttt{%s}_{%s}", b_idx, safe_var, v_idx))
        }
        b_idx <- b_idx + 1
      }
      
      lhs_mu <- get_link_latex(link_name, mu_symbol)
      lines_out <- c(lines_out, sprintf("  %s &= %s, \\nonumber", lhs_mu, chunk_terms(rhs_mu, chunk_size = 3, is_numeric = FALSE)))
    }
    
    components_sym <- list(
      phi = list(name = "phi", lhs = sprintf("\\log(\\phi_{%s})", lhs_idx), form = model$phi.formula, sym = "\\psi"),
      p0  = list(name = "p0",  lhs = sprintf("\\mathrm{logit}(p_{0%s})", lhs_idx), form = model$p0.formula, sym = "\\gamma"),
      p1  = list(name = "p1",  lhs = sprintf("\\mathrm{logit}(p_{1%s})", lhs_idx), form = model$p1.formula, sym = "\\delta")
    )
    
    for (comp in intersect(names(components_sym), comps_to_print)) {
      form <- components_sym[[comp]]$form
      if (!is.null(form) && !identical(form, ~0)) {
        mt_c <- stats::terms(form)
        tls <- attr(mt_c, "term.labels")
        sym <- components_sym[[comp]]$sym
        
        rhs_c <- c()
        b_c <- 0
        if (attr(mt_c, "intercept") == 1) {
          rhs_c <- c(rhs_c, sprintf("%s_{%d}", sym, b_c))
          b_c <- b_c + 1
        }
        
        for (var in tls) {
          print_var <- var
          if (!is.null(custom_labels) && var %in% names(custom_labels)) print_var <- custom_labels[[var]]
          safe_var <- gsub("_", "\\\\_", print_var)
          v_idx <- get_idx(var)
          
          rhs_c <- c(rhs_c, sprintf("%s_{%d}~\\texttt{%s}_{%s}", sym, b_c, safe_var, v_idx))
          b_c <- b_c + 1
        }
        lines_out <- c(lines_out, sprintf("  %s &= %s \\nonumber", components_sym[[comp]]$lhs, chunk_terms(rhs_c, chunk_size = 4, is_numeric = FALSE)))
      }
    }
    
  } else {
    summ <- summary(sd_rep, "fixed")
    components_num <- list(
      mu  = list(prefix = "beta_mu", lhs = get_link_latex(link_name, mu_symbol), mat = "X_mu"),
      phi = list(prefix = "beta_phi", lhs = sprintf("\\log(\\phi_{%s})", lhs_idx), mat = "X_phi"),
      p0  = list(prefix = "beta_p0", lhs = sprintf("\\mathrm{logit}(p_{0%s})", lhs_idx), mat = "X_p0"),
      p1  = list(prefix = "beta_p1", lhs = sprintf("\\mathrm{logit}(p_{1%s})", lhs_idx), mat = "X_p1")
    )
    
    for (comp in intersect(names(components_num), comps_to_print)) {
      prefix <- components_num[[comp]]$prefix
      idx <- grep(paste0("^", prefix), rownames(summ))
      
      if (length(idx) > 0) {
        lhs <- components_num[[comp]]$lhs
        beta_vals <- summ[idx, 1]
        mat_name <- components_num[[comp]]$mat
        
        term_names <- rep(prefix, length(beta_vals))
        if (!is.null(model$obj$env$data[[mat_name]])) {
          real_names <- colnames(model$obj$env$data[[mat_name]])
          if (!is.null(real_names) && length(real_names) == length(beta_vals)) {
            term_names <- real_names
          }
        }
        
        if (!is.null(custom_labels)) {
          matches <- match(term_names, names(custom_labels))
          term_names[!is.na(matches)] <- custom_labels[matches[!is.na(matches)]]
        }
        term_names <- gsub("_", "\\\\_", term_names) 
        
        rhs_parts <- character(length(beta_vals))
        for (k in 1:length(beta_vals)) {
          val <- beta_vals[k]
          term <- term_names[k]
          
          val_str <- sprintf(paste0("%.", digits, "f"), abs(val))
          sign_str <- if (val >= 0) "+" else "-"
          
          if (term == "(Intercept)") {
            rhs_parts[k] <- sprintf(paste0("%.", digits, "f"), val)
          } else {
            v_idx <- get_idx(term)
            rhs_parts[k] <- sprintf("%s %s~\\texttt{%s}_{%s}", sign_str, val_str, term, v_idx)
          }
        }
        
        if (comp == "mu" && has_re) {
          rhs_parts <- c(rhs_parts, "+ \\mathbf{z}_i^\\top \\mathbf{u}")
        }
        
        rhs_parts[1] <- sub("^\\+ ", "", rhs_parts[1])
        lines_out <- c(lines_out, sprintf("  %s &= %s \\nonumber", lhs, chunk_terms(rhs_parts, chunk_size = 3, is_numeric = TRUE)))
      }
    }
  }
  
  if (has_re && ("mu" %in% comps_to_print || component == "all")) {
    lines_out <- c(lines_out, "  & \\quad \\text{where } \\mathbf{u} \\sim N(\\mathbf{0}, \\mathbf{\\Sigma}) \\nonumber")
  }
  
  if (index_labels) {
    n_i <- if (has_re && !is.null(model$re_info_mu)) model$re_info_mu$n_re_levels_list[1] else model$nobs
    
    if (has_re && !is.null(grouping_var)) {
      n_j <- round(model$nobs / n_i)
      lbl <- sprintf("  & \\quad \\text{for } i = 1, \\dots, %s \\quad \\text{and } j = 1, \\dots, %s \\nonumber", n_i, n_j)
    } else {
      lbl <- sprintf("  & \\quad \\text{for } i = 1, \\dots, %s \\nonumber", n_i)
    }
    lines_out <- c(lines_out, lbl)
  }
  
  final_tex <- paste0(
    "\\begin{align}\n",
    paste(lines_out, collapse = " \\\\\n"), "\n",
    "\\end{align}\n"
  )
  
  cat(final_tex)
  invisible(final_tex)
}
