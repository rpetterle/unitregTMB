#' @name UnitWeibull
#' @aliases UnitWeibull duweibull puweibull quweibull ruweibull
#' 
#' @title The Unit-Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, and random number generation for the 
#' Unit-Weibull distribution parameterized by a quantile \code{mu} at level \code{tau} and shape \code{phi}.
#' 
#' @author Ricardo Rasmussen Petterle
#' 
#' @references 
#' Mazucheli, J., Menezes, A. F. B., & Ghitany, M. E. (2018). The unit-Weibull distribution and associated inference. 
#' \emph{Journal of Applied Statistics}, \bold{45}(8), 1351--1368.
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu quantile parameter (0 < mu < 1).
#' @param phi shape parameter (phi > 0).
#' @param tau the quantile level (0 < tau < 1). Defaults to 0.5 (median).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}.
#' 
#' @examples
#' # 1. Set parameters (Median parameterization)
#' mu_med <- 0.6
#' phi_val <- 2.5
#' 
#' # 2. Random generation
#' set.seed(123)
#' x <- ruweibull(1000, mu_med, phi_val, tau = 0.5)
#' 
#' # 3. Visualization
#' hist(x, prob = TRUE, main = "Unit-Weibull (Median)", col = "lightgray")
#' curve(duweibull(x, mu_med, phi_val, tau = 0.5), add = TRUE, col = "red", lwd = 2)
#' 
#' # 4. Check P(X <= mu) = tau
#' p_mu <- puweibull(mu_med, mu_med, phi_val, tau = 0.5)
#' cat("P(X <= mu) should be 0.5:", p_mu, "\n")
#' 
NULL

#' @rdname UnitWeibull
#' @export
duweibull <- function(x, mu, phi, tau = 0.5, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_duweibull(x, mu, phi, tau, log[1L])
}

#' @rdname UnitWeibull
#' @export
puweibull <- function(q, mu, phi, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_puweibull(q, mu, phi, tau, lower.tail[1L], log.p[1L])
}

#' @rdname UnitWeibull
#' @export
quweibull <- function(p, mu, phi, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_quweibull(p, mu, phi, tau, lower.tail[1L], log.p[1L])
}

#' @rdname UnitWeibull
#' @export
ruweibull <- function(n, mu, phi, tau = 0.5) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_ruweibull(n, mu, phi, tau)
}
