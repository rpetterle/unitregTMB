#' @name Vasicek
#' @aliases Vasicek dvasicek_mean pvasicek_mean qvasicek_mean rvasicek_mean dvasicek_quantile pvasicek_quantile qvasicek_quantile rvasicek_quantile
#' 
#' @title Vasicek Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, and random number generation for the 
#' Vasicek distribution parameterized either by the Mean or a given Quantile.
#' 
#' @author Ricardo Rasmussen Petterle
#' 
#' @references 
#' Mazucheli, J., Alves, B., Korkmaz, M. C., & Leiva, V. (2022). Vasicek Quantile and Mean 
#' Regression Models for Bounded Data: New Formulation, Mathematical Derivations, and 
#' Numerical Applications. \emph{Mathematics}, \bold{10}(9), 1389.
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mean or quantile parameter, depending on the function used (0 < mu < 1).
#' @param phi shape parameter (0 < phi < 1). See Details.
#' @param tau the quantile level (0 < tau < 1). Required only for quantile parameterizations.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}.
#' 
#' @details 
#' The Vasicek distribution is bounded on the unit interval. In this implementation, the 
#' shape parameter \code{phi} corresponds directly to the standard \eqn{\theta \in (0,1)} parameter 
#' from the literature. The mapping from the unconstrained real line is handled upstream 
#' by the \code{unitregTMB} model via a logit link.
#' 
#' @examples
#' # 1. Mean Parameterization
#' mu_m <- 0.5
#' phi_m <- 0.2 # Corresponds directly to theta = 0.2
#' x1 <- rvasicek_mean(1000, mu_m, phi_m)
#' hist(x1, prob = TRUE, main = "Vasicek (Mean)", col = "lightgray", border = "white")
#' curve(dvasicek_mean(x, mu_m, phi_m), add = TRUE, col = "red", lwd = 2)
#' 
#' # 2. Quantile Parameterization (Median)
#' mu_q <- 0.7
#' phi_q <- 0.5 # Corresponds directly to theta = 0.5
#' x2 <- rvasicek_quantile(1000, mu_q, phi_q, tau = 0.5)
#' hist(x2, prob = TRUE, main = "Vasicek (Median)", col = "lightblue", border = "white")
#' curve(dvasicek_quantile(x, mu_q, phi_q, tau = 0.5), add = TRUE, col = "blue", lwd = 2)
#' 
NULL

#' @rdname Vasicek
#' @export
dvasicek_mean <- function(x, mu, phi, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, phi < 1)
  cpp_dvasicek_mean(x, mu, phi, log[1L])
}

#' @rdname Vasicek
#' @export
pvasicek_mean <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, phi < 1)
  cpp_pvasicek_mean(q, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname Vasicek
#' @export
qvasicek_mean <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, phi < 1)
  cpp_qvasicek_mean(p, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname Vasicek
#' @export
rvasicek_mean <- function(n, mu, phi) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0, phi < 1)
  cpp_rvasicek_mean(n, mu, phi)
}

#' @rdname Vasicek
#' @export
dvasicek_quantile <- function(x, mu, phi, tau = 0.5, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, phi < 1, tau > 0, tau < 1)
  cpp_dvasicek_quantile(x, mu, phi, tau, log[1L])
}

#' @rdname Vasicek
#' @export
pvasicek_quantile <- function(q, mu, phi, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, phi < 1, tau > 0, tau < 1)
  cpp_pvasicek_quantile(q, mu, phi, tau, lower.tail[1L], log.p[1L])
}

#' @rdname Vasicek
#' @export
qvasicek_quantile <- function(p, mu, phi, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, phi < 1, tau > 0, tau < 1)
  cpp_qvasicek_quantile(p, mu, phi, tau, lower.tail[1L], log.p[1L])
}

#' @rdname Vasicek
#' @export
rvasicek_quantile <- function(n, mu, phi, tau = 0.5) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0, phi < 1, tau > 0, tau < 1)
  cpp_rvasicek_quantile(n, mu, phi, tau)
}
