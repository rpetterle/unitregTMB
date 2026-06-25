#' @name UnitGamma
#' @aliases UnitGamma dugamma_mean pugamma_mean qugamma_mean rugamma_mean dugamma_mode pugamma_mode qugamma_mode rugamma_mode
#' 
#' @title The Unit-Gamma Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, and random number generation for the 
#' Unit-Gamma distribution parameterized by the Mean or Mode.
#' 
#' @author Ricardo Rasmussen Petterle
#' 
#' @references 
#' Grassia, A. (1977). On a family of distributions with argument between 0 and 1 obtained by transformation of the gamma and beta families. 
#' \emph{Australian Journal of Statistics}, \bold{19}(2), 108--114.
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mean or mode parameter, depending on the function used (0 < mu < 1).
#' @param phi shape parameter (phi > 0).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}.
#' 
#' @details 
#' The Unit-Gamma distribution is obtained by the transformation \eqn{Y = \exp(-X)}, where \eqn{X \sim \Gamma(\text{shape}=\phi, \text{rate}=\beta)}.
#' 
#' @examples
#' # 1. Mean Parameterization
#' mu_m <- 0.5
#' phi_m <- 2.0
#' x1 <- rugamma_mean(1000, mu_m, phi_m)
#' hist(x1, prob = TRUE, main = "Unit-Gamma (Mean)", col = "lightgray", border = "white")
#' curve(dugamma_mean(x, mu_m, phi_m), add = TRUE, col = "red", lwd = 2)
#' 
#' # 2. Mode Parameterization
#' mu_d <- 0.7
#' phi_d <- 5.0
#' x2 <- rugamma_mode(1000, mu_d, phi_d)
#' hist(x2, prob = TRUE, main = "Unit-Gamma (Mode)", col = "lightblue", border = "white")
#' curve(dugamma_mode(x, mu_d, phi_d), add = TRUE, col = "blue", lwd = 2)
#' 
NULL

# --- Mean Parameterization ---

#' @rdname UnitGamma
#' @export
dugamma_mean <- function(x, mu, phi, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_dugamma_mean(x, mu, phi, log[1L])
}

#' @rdname UnitGamma
#' @export
pugamma_mean <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_pugamma_mean(q, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname UnitGamma
#' @export
qugamma_mean <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_qugamma_mean(p, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname UnitGamma
#' @export
rugamma_mean <- function(n, mu, phi) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0)
  cpp_rugamma_mean(n, mu, phi)
}

# --- Mode Parameterization ---

#' @rdname UnitGamma
#' @export
dugamma_mode <- function(x, mu, phi, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_dugamma_mode(x, mu, phi, log[1L])
}

#' @rdname UnitGamma
#' @export
pugamma_mode <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_pugamma_mode(q, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname UnitGamma
#' @export
qugamma_mode <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_qugamma_mode(p, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname UnitGamma
#' @export
rugamma_mode <- function(n, mu, phi) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0)
  cpp_rugamma_mode(n, mu, phi)
}
