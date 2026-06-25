#' Plot Conditional Modes (Caterpillar Plot) for unitregTMB
#' 
#' @param object A fitted \code{unitregTMB} model.
#' @param level Confidence level for the error bars. Default is 0.95.
#' @param ... Additional graphical parameters passed to \code{plot}.
#' 
#' @importFrom stats qnorm
#' @importFrom graphics plot abline segments points
#' @export
plot_ranef <- function(object, level = 0.95, ...) {
  
  if (!object$has_random_effects_mu) {
    stop("This model does not contain random effects.")
  }
  
  if (is.null(object$sd_report)) {
    stop("Standard errors (sd_report) are required but not available in this object.")
  }
  
  # Extração dos Modos Condicionais e Erros Padrão
  sdr <- summary(object$sd_report, "random")
  u_idx <- grepl("^u_mu", rownames(sdr))
  
  blups <- sdr[u_idx, "Estimate"]
  ses   <- sdr[u_idx, "Std. Error"]
  
  # Ordenação para o efeito "Caterpillar"
  ord <- order(blups)
  blups_ord <- blups[ord]
  ses_ord   <- ses[ord]
  
  # Intervalos de Confiança
  alpha <- (1 - level) / 2
  z_val <- stats::qnorm(1 - alpha)
  lwr <- blups_ord - z_val * ses_ord
  upr <- blups_ord + z_val * ses_ord
  
  n_re <- length(blups_ord)
  
  # Configuração do gráfico Base R (estilo limpo)
  graphics::plot(blups_ord, 1:n_re, 
                 xlim = range(c(lwr, upr), na.rm = TRUE), 
                 type = "n",
                 xlab = "Conditional Modes (BLUPs)", 
                 ylab = "Group Levels (Sorted)",
                 main = paste("Random Effects -", object$family),
                 yaxt = "n", bty = "n", ...)
  
  # Grelha de fundo opcional para facilitar a leitura
  graphics::grid(nx = NULL, ny = NA, col = "gray92", lty = 1)
  
  # Linha de referência no zero
  graphics::abline(v = 0, lty = 2, col = "#D55E00", lwd = 1.5)
  
  # Segmentos de incerteza (barras de erro horizontais)
  graphics::segments(x0 = lwr, y0 = 1:n_re, x1 = upr, y1 = 1:n_re, 
                     col = "gray50", lwd = 1.2)
  
  # Pontos dos BLUPs
  graphics::points(blups_ord, 1:n_re, pch = 16, col = "#0072B2", cex = 0.9)
  
  invisible(data.frame(Estimate = blups_ord, Std.Error = ses_ord, Lower = lwr, Upper = upr))
}