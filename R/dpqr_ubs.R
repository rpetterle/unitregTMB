#' @name UBS
#' @aliases UBS dubs pubs qubs rubs
#' 
#' @title Unit-Birnbaum-Saunders (UBS) Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, and random number generation for the 
#' Unit-Birnbaum-Saunders (UBS) distribution parameterized by a quantile \code{mu} at level \code{tau}.
#' 
#' @author Ricardo Rasmussen Petterle
#' 
#' @references 
#' Mazucheli, J., Menezes, A. F. B., & Dey, S. (2018). Unit-Birnbaum-Saunders distribution with applications. 
#' \emph{Chilean Journal of Statistics}, \bold{9}(1), 47--57.
#'
#' Mazucheli, J., Alves, B. and Menezes, A. F. B., (2021). A new quantile regression for modeling bounded data under 
#' a unit Birnbaum-Saunders distribution with applications. Simmetry, (), 1--28.
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
#' # 1. Set parameters (Median)
#' mu_med <- 0.5
#' phi_val <- 0.5
#' 
#' # 2. Generate random values
#' set.seed(123)
#' x <- rubs(1000, mu_med, phi_val, tau = 0.5)
#' 
#' # 3. Visualization
#' hist(x, prob = TRUE, main = "UBS (Median 0.5)", col = "lavender")
#' curve(dubs(x, mu_med, phi_val, tau = 0.5), add = TRUE, col = "purple", lwd = 2)
#' 
NULL

#' @rdname UBS
#' @export
dubs <- function(x, mu, phi, tau = 0.5, log = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_dubs(x, mu, phi, tau, log[1L])
}

#' @rdname UBS
#' @export
pubs <- function(q, mu, phi, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_pubs(q, mu, phi, tau, lower.tail[1L], log.p[1L])
}

#' @rdname UBS
#' @export
qubs <- function(p, mu, phi, tau = 0.5, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_qubs(p, mu, phi, tau, lower.tail[1L], log.p[1L])
}

#' @rdname UBS
#' @export
rubs <- function(n, mu, phi, tau = 0.5) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0, tau > 0, tau < 1)
  cpp_rubs(n, mu, phi, tau)
}
