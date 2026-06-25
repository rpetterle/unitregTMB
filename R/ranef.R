#' Extract Random Effects (BLUPs) from a unitregTMB Object
#'
#' @param object A fitted object of class "unitregTMB".
#' @param ... Additional arguments (not used).
#' @return A list (one element per random term/group) 
#'         containing data.frames of the random effects.
#' @export
#' 
ranef.unitregTMB <- function(object, ...) {
  
  # 1. Check if the sd_report exists and if random effects are present
  if (!isTRUE(object$has_random_effects_mu)) {
    warning("The model was fitted without random effects for mu.")
    return(NULL)
  }
  
  if (is.null(object$sd_report)) {
    stop("sd_report not found in the object. ",
         "Please re-fit the model with 'control = unitregTMB.control(sd.report = TRUE)'.")
  }
  
  # 2. Extract the summary of "random" parameters
  blups_matrix <- tryCatch(
    summary(object$sd_report, "random"),
    error = function(e) NULL
  )
  
  if (is.null(blups_matrix) || !is.matrix(blups_matrix)) {
    stop("Failed to extract the 'random' matrix from sd_report. ",
         "The output was not the expected matrix.")
  }
  
  # 3. Get the vector of estimates (BLUPs)
  blup_estimates_vector <- blups_matrix[, "Estimate"]
  
  # 4. Get structural information from the 'unitregTMB' object
  re_info <- object$re_info_mu
  
  if (is.null(re_info) || !isTRUE(re_info$has_random_effects)) {
    stop("The 're_info_mu' object is missing or malformed.")
  }
  
  # 5. Reshape the BLUPs into the standard ranef list format
  output_list <- list()
  current_u_start_idx <- 1 
  
  for (t_idx in seq_len(re_info$n_re_terms)) {
    
    group_var_name <- re_info$group_var_names[t_idx]  
    group_levels <- re_info$group_levels[[t_idx]]     
    n_levels <- re_info$n_re_levels_list[t_idx]       
    n_comp <- re_info$re_term_n_components[t_idx]     
    comp_names <- re_info$cnms[[t_idx]]               
    
    expected_length <- n_levels * n_comp 
    
    if (current_u_start_idx + expected_length - 1 > length(blup_estimates_vector)) {
      stop(paste("Dimension mismatch while processing term:", group_var_name))
    }
    
    # -------------------------------------------------------------------------
    # CRITICAL FIX: Matrix Reshaping
    # KhatriRao and TMB C++ loop group the vector by Level, then Component.
    # e.g., [Level1_Int, Level1_Slope, Level2_Int, Level2_Slope, ...]
    # We must fill the matrix row by row (byrow = TRUE) to align components to columns.
    # -------------------------------------------------------------------------
    term_blups <- blup_estimates_vector[current_u_start_idx:(current_u_start_idx + expected_length - 1)]
    
    blups_mat <- matrix(term_blups, nrow = n_levels, ncol = n_comp, byrow = TRUE)
    
    blups_df_this_group <- as.data.frame(blups_mat)
    colnames(blups_df_this_group) <- comp_names
    rownames(blups_df_this_group) <- group_levels
    
    output_list[[group_var_name]] <- blups_df_this_group
    
    current_u_start_idx <- current_u_start_idx + expected_length
  }
  
  return(output_list)
}
