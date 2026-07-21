#' @name BetaMode
#' @aliases BetaMode dbeta_mode pbeta_mode qbeta_mode rbeta_mode
#' 
#' @title Beta Distribution (Mode Parameterization) 
#' 
#' @description 
#' Density, distribution function, quantile function, and random number generation for the 
#' Beta distribution parameterized by the mode \code{mu} and precision \code{phi}.
#' 
#' @author Ricardo Rasmussen Petterle
#' 
#' @details 
#' The Beta distribution with mode \eqn{\mu} and precision \eqn{\phi} has shape parameters:
#' \deqn{\alpha = \mu\phi + 1}
#' \deqn{\beta = (1 - \mu)\phi + 1}
#' where \eqn{0 < \mu < 1} and \eqn{\phi > 0}.
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu vector of mode parameters (0 < mu < 1).
#' @param phi vector of precision parameters (phi > 0).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}.
#' 
#' @examples
#' # 1. Set parameters
#' mu_val <- 0.7
#' phi_val <- 10
#' 
#' # 2. Generate random values
#' set.seed(123)
#' x <- rbeta_mode(1000, mu_val, phi_val)
#' 
#' # 3. Visualization
#' hist(x, prob = TRUE, main = "Beta (Mode 0.7)", col = "lightgray")
#' curve(dbeta_mode(x, mu_val, phi_val), add = TRUE, col = "red", lwd = 2)
#' abline(v = mu_val, col = "blue", lty = 2) # Mode line
#' 
NULL

#' @rdname BetaMode
#' @export
dbeta_mode <- function(x, mu, phi, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_dbeta_mode(x, mu, phi, log[1L])
}

#' @rdname BetaMode
#' @export
pbeta_mode <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_pbeta_mode(q, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname BetaMode
#' @export
qbeta_mode <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_qbeta_mode(p, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname BetaMode
#' @export
rbeta_mode <- function(n, mu, phi) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0)
  cpp_rbeta_mode(n, mu, phi)
}
