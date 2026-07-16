#' @name plot_ranef.unitregTMB
#' 
#' @title Plot Conditional Modes (Caterpillar Plot) for unitregTMB
#' 
#' @description 
#' Produces a caterpillar plot of the conditional modes (BLUPs) of the random effects 
#' from a fitted \code{unitregTMB} model, along with their confidence intervals.
#' 
#' @param object A fitted \code{unitregTMB} model.
#' @param level Confidence level for the error bars. Default is 0.95.
#' @param ... Additional graphical parameters passed to \code{plot}.
#' 
#' @return Invisibly returns a \code{data.frame} containing the sorted estimated random effects 
#'   (BLUPs), their standard errors, and the lower and upper bounds of the confidence intervals.
#' 
#' @examples
#' \donttest{
#' # Assuming 'da' is your dataset and 'id' is a grouping variable:
#' # fit <- unitregTMB(Y ~ educ + (1 | id), data = da, family = vasicek())
#' # 
#' # # Generate the caterpillar plot
#' # plot_ranef(fit, level = 0.95)
#' # 
#' # # You can also save the data used for the plot
#' # re_data <- plot_ranef(fit)
#' # head(re_data)
#' }
#' 
#' @importFrom stats qnorm
#' @importFrom graphics plot abline segments points grid
#' @rdname plot_ranef_unitregTMB
#' @export
plot_ranef.unitregTMB <- function(object, level = 0.95, ...) {
  
  if (!object$has_random_effects_mu) {
    stop("This model does not contain random effects.")
  }
  
  if (is.null(object$sd_report)) {
    stop("Standard errors (sd_report) are required but not available in this object.")
  }
  
  sdr <- summary(object$sd_report, "random")
  u_idx <- grepl("^u_mu", rownames(sdr))
  
  blups <- sdr[u_idx, "Estimate"]
  ses   <- sdr[u_idx, "Std. Error"]
  
  ord <- order(blups)
  blups_ord <- blups[ord]
  ses_ord   <- ses[ord]
  
  alpha <- (1 - level) / 2
  z_val <- stats::qnorm(1 - alpha)
  lwr <- blups_ord - z_val * ses_ord
  upr <- blups_ord + z_val * ses_ord
  
  n_re <- length(blups_ord)
  
  graphics::plot(blups_ord, 1:n_re, 
                 xlim = range(c(lwr, upr), na.rm = TRUE), 
                 type = "n",
                 xlab = "Conditional Modes (BLUPs)", 
                 ylab = "Group Levels (Sorted)",
                 main = paste("Random Effects -", object$family),
                 yaxt = "n", bty = "n", ...)
  
  graphics::grid(nx = NULL, ny = NA, col = "gray92", lty = 1)
  
  graphics::abline(v = 0, lty = 2, col = "#D55E00", lwd = 1.5)
  
  graphics::segments(x0 = lwr, y0 = 1:n_re, x1 = upr, y1 = 1:n_re, 
                     col = "gray50", lwd = 1.2)
  
  graphics::points(blups_ord, 1:n_re, pch = 16, col = "#0072B2", cex = 0.9)
  
  invisible(data.frame(Estimate = blups_ord, Std.Error = ses_ord, Lower = lwr, Upper = upr))
}
