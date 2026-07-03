#' @name Simplex
#' @aliases Simplex dsimplex psimplex qsimplex rsimplex
#' 
#' @title Simplex Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, and random number generation for the 
#' Simplex distribution parameterized by the mean \code{mu} and dispersion \code{phi}.
#' 
#' @author Ricardo Rasmussen Petterle
#' 
#' @references 
#' Barndorff-Nielsen, O. E., & Jørgensen, B. (1991). Some parametric models on the simplex. 
#' \emph{Journal of Multivariate Analysis}, \bold{39}(1), 106--116.
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mean parameter (0 < mu < 1).
#' @param phi dispersion parameter (phi > 0).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}.
#' 
#' @details 
#' The Simplex distribution is parameterized such that \eqn{\mu} is the mean 
#' and \eqn{\phi = \sigma^2} is the dispersion parameter. The variance of the 
#' distribution is approximately \eqn{Var(Y) \approx \phi \mu^3 (1-\mu)^3}.
#' 
#' @examples
#' # 1. Set parameters
#' mu_val <- 0.5
#' phi_val <- 1.5
#' 
#' # 2. Generate random values
#' set.seed(123)
#' x <- rsimplex(1000, mu_val, phi_val)
#' 
#' # 3. Visualization
#' hist(x, prob = TRUE, main = "Simplex (Mean 0.5)", col = "lightgreen")
#' curve(dsimplex(x, mu_val, phi_val), add = TRUE, col = "darkgreen", lwd = 2)
#' 
NULL

#' @rdname Simplex
#' @export
dsimplex <- function(x, mu, phi, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_dsimplex(x, mu, phi, log[1L])
}

#' @rdname Simplex
#' @export
psimplex <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_psimplex(q, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname Simplex
#' @export
qsimplex <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_qsimplex(p, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname Simplex
#' @export
rsimplex <- function(n, mu, phi) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0)
  cpp_rsimplex(n, mu, phi)
}
