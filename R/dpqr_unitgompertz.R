#' @name UnitGompertz
#' @aliases UnitGompertz dugompertz_mode pugompertz_mode qugompertz_mode rugompertz_mode dugompertz_quantile pugompertz_quantile qugompertz_quantile rugompertz_quantile
#' 
#' @title Unit-Gompertz Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, and random number generation for the 
#' Unit-Gompertz distribution parameterized either by the Mode or a given Quantile.
#' 
#' @author Ricardo Rasmussen Petterle
#' 
#' @references 
#' Mazucheli, J., Menezes, A. F. B., & Dey, S. (2019). The unit-Gompertz distribution with applications. 
#' \emph{Statistica}, \bold{79}(1), 25--43.
#' 
#' Menezes, A. F., Mazucheli, J. and Chakraborty, S. (2021). A collection of parametric modal regression models for bounded data.
#' \emph{Journal of Biopharmaceutical Statistics}, \bold{31}(4), 490--506. 
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu mode or quantile parameter, depending on the function used (0 < mu < 1).
#' @param phi shape parameter (phi > 0).
#' @param tau the quantile level (0 < tau < 1). Required only for quantile parameterizations.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}.
#' 
#' @examples
#' # 1. Mode Parameterization
#' mu_mode <- 0.5
#' phi_val <- 2.0
#' x1 <- rugompertz_mode(1000, mu_mode, phi_val)
#' hist(x1, prob = TRUE, main = "Unit-Gompertz (Mode)", col = "lightgray")
#' curve(dugompertz_mode(x, mu_mode, phi_val), add = TRUE, col = "red", lwd = 2)
#' 
#' # 2. Quantile Parameterization (Median)
#' mu_quant <- 0.7
#' phi_val <- 1.5
#' x2 <- rugompertz_quantile(1000, mu_quant, phi_val, tau = 0.5)
#' hist(x2, prob = TRUE, main = "Unit-Gompertz (Median)", col = "lightblue")
#' curve(dugompertz_quantile(x, mu_quant, phi_val, tau = 0.5), add = TRUE, col = "blue", lwd = 2)
#' 
NULL

#' @rdname UnitGompertz
#' @export
dugompertz_mode <- function(x, mu, phi, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_dugompertz_mode(x, mu, phi, log[1L])
}

#' @rdname UnitGompertz
#' @export
pugompertz_mode <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_pugompertz_mode(q, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname UnitGompertz
#' @export
qugompertz_mode <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0)
  cpp_qugompertz_mode(p, mu, phi, lower.tail[1L], log.p[1L])
}

#' @rdname UnitGompertz
#' @export
rugompertz_mode <- function(n, mu, phi) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0)
  cpp_rugompertz_mode(n, mu, phi)
}

#' @rdname UnitGompertz
#' @export
dugompertz_quantile <- function(x, mu, phi, tau = 0.5, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_dugompertz_quantile(x, mu, phi, tau, log[1L])
}

#' @rdname UnitGompertz
#' @export
pugompertz_quantile <- function(q, mu, phi, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_pugompertz_quantile(q, mu, phi, tau, lower.tail[1L], log.p[1L])
}

#' @rdname UnitGompertz
#' @export
qugompertz_quantile <- function(p, mu, phi, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_qugompertz_quantile(p, mu, phi, tau, lower.tail[1L], log.p[1L])
}

#' @rdname UnitGompertz
#' @export
rugompertz_quantile <- function(n, mu, phi, tau = 0.5) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_rugompertz_quantile(n, mu, phi, tau)
}
