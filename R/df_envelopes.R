#' Extract Envelope Data into a Data Frame for unitregTMB Models
#'
#' @description 
#' Method for the \code{unitcore::df_envelopes} generic. Extracts the coordinates of 
#' simulated envelopes from one or more fitted \code{unitregTMB} models. The result 
#' is a tidy data frame in long format, highly optimized for creating faceted 
#' graphics using \code{ggplot2}.
#'
#' @param object An object of class \code{unitregTMB}.
#' @param ... Additional objects of class \code{unitregTMB} to be compared.
#' @param nsim Number of simulations used to compute the envelope. Default is 50.
#' @param halfnormal logical. If \code{TRUE} (default), a half-normal plot is produced. 
#'        If \code{FALSE}, a normal plot is produced. Ignored if \code{resid.type = "cox-snell"}.
#' @param resid.type Type of residuals to be used: \code{"quantile"} (default) or \code{"cox-snell"}.
#' @param ncpus Number of cores to use for parallel computing. Default is 1.
#' @param show.progress logical. Should the progress bar and messages be displayed? Default is \code{TRUE}.
#'
#' @return A \code{data.frame} containing the simulated envelope coordinates 
#'         (\code{Theoretical}, \code{Observed}, \code{Lower}, \code{Median}, \code{Upper}) 
#'         for all provided models. A factor column \code{Model} is included to distinguish 
#'         the source of each set of residuals, properly ordered.
#'         
#' @importFrom unitcore df_envelopes plot_envelope
#' @export
#'
#' @examples
#' \donttest{
#' # Assuming 'da' is your dataset
#' # fit_vas <- unitregTMB(Y ~ x, family = vasicek(), data = da)
#' # fit_bet <- unitregTMB(Y ~ x, family = beta_fam(), data = da)
#'
#' # Extract tidy data frame for ggplot2
#' # df_plot <- df_envelopes(fit_vas, fit_bet, nsim = 19)
#' 
#' # Example ggplot2 code:
#' # library(ggplot2)
#' # ggplot(df_plot, aes(x = Theoretical)) +
#' #   geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey85", alpha = 0.5) +
#' #   geom_point(aes(y = Observed)) +
#' #   facet_wrap(~ Model)
#' }
df_envelopes.unitregTMB <- function(object, ..., nsim = 50, halfnormal = TRUE, 
                                    resid.type = c("quantile", "cox-snell"), 
                                    ncpus = 1, show.progress = TRUE) {
  
  resid.type <- match.arg(resid.type)
  models <- list(object, ...)
  
  is_valid_model <- sapply(models, inherits, "unitregTMB")
  if (!all(is_valid_model)) {
    stop("All unnamed arguments must be 'unitregTMB' objects.", call. = FALSE)
  }
  
  model_names <- sapply(models, function(m) {
    if (is.list(m$family)) {
      as.character(m$family$name)
    } else if (!is.null(m$family)) {
      as.character(m$family)
    } else {
      "unitregTMB model"
    }
  })
  
  if (any(duplicated(model_names))) {
    model_names <- make.unique(model_names)
  }
  
  df_list <- lapply(seq_along(models), function(i) {
    mod <- models[[i]]
    m_name <- model_names[i]
    
    if (show.progress) {
      cat(sprintf("\nExtracting envelope data for: %s\n", m_name))
    }
    
    env_res <- plot_envelope(mod, plot = FALSE, nsim = nsim, halfnormal = halfnormal, 
                             resid.type = resid.type, ncpus = ncpus, show.progress = show.progress)
    
    if (is.null(env_res) || is.null(env_res$data)) {
      stop(sprintf("Failed to extract envelope data for model: %s", m_name), call. = FALSE)
    }
    
    df <- env_res$data
    df$Model <- factor(m_name, levels = model_names)
    
    return(df)
  })
  
  df_final <- do.call(rbind, df_list)
  rownames(df_final) <- NULL
  
  return(df_final)
}