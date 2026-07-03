#' @name BetaMean
#' @aliases BetaMean dbeta_mean pbeta_mean qbeta_mean rbeta_mean
#' 
#' @title Beta Distribution (Mean Parameterization)
#' 
#' @description 
#' Density, distribution function, quantile function, and random number generation for the 
#' Beta distribution parameterized by the mean \code{mu} and precision \code{phi}.
#' 
#' @author Ricardo Rasmussen Petterle
#' 
#' @references 
#' Ferrari, S. L. P., & Cribari-Neto, F. (2004). Beta regression for modeling rates and proportions. 
#' \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu vector of mean parameters (0 < mu < 1).
#' @param phi vector of precision parameters (phi > 0).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#' 
#' @return 
#' \code{dbeta_mean} gives the density, \code{pbeta_mean} gives the distribution function,
#' \code{qbeta_mean} gives the quantile function, and \code{rbeta_mean} generates random deviates.
#' 
#' @details 
#' The Beta distribution with mean \eqn{\mu} and precision \eqn{\phi} has shape parameters:
#' \deqn{\alpha = \mu \phi}
#' \deqn{\beta = (1 - \mu) \phi}
#' where \eqn{0 < \mu < 1} and \eqn{\phi > 0}.
#' 
#' @examples
#' # 1. Set parameters
#' mu_val <- 0.6
#' phi_val <- 10
#' n_sim <- 1000
#' 
#' # 2. Generate random values
#' set.seed(123)
#' x <- rbeta_mean(n = n_sim, mu = mu_val, phi = phi_val)
#' 
#' # 3. Visual Diagnostics
#' oldpar <- par(mfrow = c(1, 3))
#' 
#' # Panel A: Histogram + Density
#' hist(x, prob = TRUE, main = "Density", col = "lightgray", border = "white")
#' curve(dbeta_mean(x, mu = mu_val, phi = phi_val), add = TRUE, col = "blue", lwd = 2)
#' legend("topleft", legend = "Theoretical", col = "blue", lwd = 2, bty = "n")
#' 
#' # Panel B: ECDF + CDF
#' plot(ecdf(x), main = "CDF", col = "black", lwd = 1)
#' curve(pbeta_mean(x, mu = mu_val, phi = phi_val), add = TRUE, col = "red", lwd = 2)
#' legend("topleft", legend = "Theoretical", col = "red", lwd = 2, bty = "n")
#' 
#' # Panel C: Quantile vs Quantile (P-P plot equivalent logic)
#' probs <- seq(0, 1, length.out = 100)
#' q_theo <- qbeta_mean(probs, mu = mu_val, phi = phi_val)
#' plot(probs, q_theo, type = "l", col = "red", lwd = 2, 
#'      main = "Quantile Function", xlab = "Probability", ylab = "Quantile")
#' points(probs, quantile(x, probs), col = "black", pch = 20, cex = 0.5)
#' legend("topleft", legend = c("Theoretical", "Empirical"), 
#'        col = c("red", "black"), lty = c(1, NA), pch = c(NA, 20), bty = "n")
#' 
#' par(oldpar)
#' 
NULL

#' @rdname BetaMean
#' @export
dbeta_mean <- function(x, mu, phi, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_dbeta_mean(x, mu, phi, log[1L])
}

#' @rdname BetaMean
#' @export
pbeta_mean <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_pbeta_mean(q, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname BetaMean
#' @export
qbeta_mean <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_qbeta_mean(p, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname BetaMean
#' @export
rbeta_mean <- function(n, mu, phi) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0)
  cpp_rbeta_mean(n, mu, phi)
}
