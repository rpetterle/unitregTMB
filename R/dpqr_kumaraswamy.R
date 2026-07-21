#' @name Kumaraswamy
#' @aliases Kumaraswamy dkum_mode pkum_mode qkum_mode rkum_mode dkum_quantile pkum_quantile qkum_quantile rkum_quantile
#' 
#' @title Kumaraswamy Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, and random number generation for the 
#' Kumaraswamy distribution parameterized by the Mode or a Quantile.
#' 
#' @author Ricardo Rasmussen Petterle
#' 
#' @references 
#' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes. 
#' \emph{Journal of Hydrology}, \bold{46}(1-2), 79--88.
#' 
#' Jones, M. C., (2009). Kumaraswamy's distribution: A beta-type distribution with some tractability advantages. Statistical Methodology, 6(1), 70-81.
#'
#' Menezes, A. F., Mazucheli, J. and Chakraborty, S. (2021). A collection of parametric modal regression models 
#' for bounded data, Journal of Biopharmaceutical Statistics 31(4): 490–506.
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu location parameter (Mode or Quantile), depending on the function (0 < mu < 1).
#' @param phi shape parameter 'a' (phi > 0, note that phi > 1 for Mode parametrization).
#' @param tau the quantile level to be modeled (only for quantile parametrization).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}.
#' 
#' @details 
#' The Kumaraswamy distribution has PDF:
#' \deqn{f(x) = a b x^{a-1} (1 - x^a)^{b-1}}
#' 
NULL

#' @rdname Kumaraswamy
#' @export
dkum_mode <- function(x, mu, phi, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 1)
  cpp_dkum_mode(x, mu, phi, log[1L])
}

#' @rdname Kumaraswamy
#' @export
pkum_mode <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 1)
  cpp_pkum_mode(q, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname Kumaraswamy
#' @export
qkum_mode <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 1)
  cpp_qkum_mode(p, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname Kumaraswamy
#' @export
rkum_mode <- function(n, mu, phi) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 1)
  cpp_rkum_mode(n, mu, phi)
}

#' @rdname Kumaraswamy
#' @export
dkum_quantile <- function(x, mu, phi, tau = 0.5, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_dkum_quantile(x, mu, phi, tau, log[1L])
}

#' @rdname Kumaraswamy
#' @export
pkum_quantile <- function(q, mu, phi, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_pkum_quantile(q, mu, phi, tau, lower.tail[1L], log.p[1L])
}

#' @rdname Kumaraswamy
#' @export
qkum_quantile <- function(p, mu, phi, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_qkum_quantile(p, mu, phi, tau, lower.tail[1L], log.p[1L])
}

#' @rdname Kumaraswamy
#' @export
rkum_quantile <- function(n, mu, phi, tau = 0.5) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_rkum_quantile(n, mu, phi, tau)
}
