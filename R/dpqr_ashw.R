#' @name ASHW
#' @aliases ASHW dashw pashw qashw rashw
#' 
#' @title Arcsecant Hyperbolic Weibull (ASHW) Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, and random number generation for the 
#' Arcsecant Hyperbolic Weibull (ASHW) distribution parameterized by a quantile \code{mu} at level \code{tau}.
#' 
#' @author Ricardo Rasmussen Petterle
#' 
#' @references 
#' Kieschnick, R., & McCullough, B. D. (2003). Regression analysis of variates observed on (0, 1): percentages, proportions and fractions. 
#' \emph{Statistical Modelling}, \bold{3}(3), 193--213.
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
#' mu_med <- 0.5
#' phi_val <- 1.5
#' 
#' # 2. Random generation
#' set.seed(123)
#' x <- rashw(1000, mu_med, phi_val, tau = 0.5)
#' 
#' # 3. Visualization
#' hist(x, prob = TRUE, main = "ASHW (Median 0.5)", col = "lightgray")
#' curve(dashw(x, mu_med, phi_val, tau = 0.5), add = TRUE, col = "red", lwd = 2)
#' 
#' # 4. Check P(X <= mu) = tau
#' p_mu <- pashw(mu_med, mu_med, phi_val, tau = 0.5)
#' cat("P(X <= mu) should be 0.5:", p_mu, "\n")
#' 
NULL

#' @rdname ASHW
#' @export
dashw <- function(x, mu, phi, tau = 0.5, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_dashw(x, mu, phi, tau, log[1L])
}

#' @rdname ASHW
#' @export
pashw <- function(q, mu, phi, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_pashw(q, mu, phi, tau, lower.tail[1L], log.p[1L])
}

#' @rdname ASHW
#' @export
qashw <- function(p, mu, phi, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_qashw(p, mu, phi, tau, lower.tail[1L], log.p[1L])
}

#' @rdname ASHW
#' @export
rashw <- function(n, mu, phi, tau = 0.5) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_rashw(n, mu, phi, tau)
}
