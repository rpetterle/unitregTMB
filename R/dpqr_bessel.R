#' @name Bessel
#' @aliases Bessel dbessel pbessel qbessel rbessel
#' 
#' @title Bessel Distribution
#'
#' @description 
#' Density, distribution function, quantile function and random generation for 
#' the Bessel regression distribution with mean \code{mu} and precision \code{phi}.
#'
#' @author Ricardo Rasmussen Petterle
#'
#' @references 
#' Barreto-Souza, W., Mayrink, V. D., & Simas, A. B. (2022). 
#' Bessel regression model: a robust alternative to beta regression. 
#' \emph{Journal of Applied Statistics}, \bold{49}(1), 1--20.
#'
#' @param x,q vector of quantiles (0 < x, q < 1).
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu mean parameter (0 < mu < 1).
#' @param phi precision parameter (phi > 0).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @details
#' The probability density function for the Bessel distribution is given by:
#' \deqn{f(y \mid \mu, \phi) = \frac{\mu(1-\mu) \phi e^{\phi}}{\pi [y(1-y)]^{3/2} \zeta_{\mu}(y)} K_1(\phi \zeta_{\mu}(y))}
#' where \eqn{\zeta_{\mu}(y) = \sqrt{1 + \frac{(y-\mu)^2}{y(1-y)}}} and \eqn{K_1} is the modified Bessel function of the second kind.
#'
#' The random generation utilizes a highly efficient C++ implementation of the Inverse Gaussian mixture property proposed by Giner and Smyth (2016).
#'
#' @examples
#' # 1. Density Evaluation
#' dbessel(0.5, mu = 0.5, phi = 10)
#' 
#' # 2. Cumulative Probability
#' pbessel(0.5, mu = 0.5, phi = 10)
#' 
#' # 3. Random Generation & Histogram
#' set.seed(123)
#' r <- rbessel(1000, mu = 0.6, phi = 5)
#' hist(r, prob = TRUE, breaks = 30, col = "grey", main = "Bessel(0.6, 5)")
#' curve(dbessel(x, mu = 0.6, phi = 5), add = TRUE, col = "red", lwd = 2)
#' 
NULL

#' @rdname Bessel
#' @export
dbessel <- function(x, mu, phi, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_dbessel(x, mu, phi, log[1L])
}

#' @rdname Bessel
#' @export
pbessel <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_pbessel(q, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname Bessel
#' @export
qbessel <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_qbessel(p, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname Bessel
#' @export
rbessel <- function(n, mu, phi) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0)
  cpp_rbessel(n, mu, phi)
}
