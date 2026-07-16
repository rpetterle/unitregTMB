#' @title Vuong Test for Non-Nested unitregTMB Models
#'
#' @description Performs the Vuong likelihood ratio test for model selection 
#' between two fitted non-nested models of class \code{unitregTMB}.
#'
#' @param object1,object2 Objects of class \code{unitregTMB} containing the fitted models.
#' @param alternative Indicates the alternative hypothesis and must be one
#' of \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}. You can
#' specify just the initial letter of the value.
#' @param ... Additional arguments passed to methods.
#'
#' @details The statistic of Vuong likelihood ratio test for compare two
#' non-nested regression models is defined by
#' \deqn{T = \frac{1}{\widehat{\omega}^2\,\sqrt{n}}\,\sum_{i=1}^{n}\,
#' \log\frac{f(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\theta}})}{
#' g(y_i \mid \boldsymbol{x}_i,\widehat{\boldsymbol{\gamma}})}}
#' where
#' \deqn{\widehat{\omega}^2 = \frac{1}{n}\,\sum_{i=1}^{n}\,\left(\log \frac{f(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\theta}})}{g(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\gamma}})}\right)^2 - \left[\frac{1}{n}\,\sum_{i=1}^{n}\,\left(\log \frac{f(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\theta}})}{ g(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\gamma}})}\right)\right]^2}
#' is an estimator for the variance of the log-likelihood ratio.
#'
#' When \eqn{n \rightarrow \infty} we have that \eqn{T \rightarrow N(0, 1)} in distribution.
#' Therefore, at \eqn{\alpha\%} level of significance, the null hypothesis of
#' the equivalence of the competing models is rejected if \eqn{|T| > z_{\alpha/2}},
#' where \eqn{z_{\alpha/2}} is the \eqn{\alpha/2} quantile of standard normal distribution.
#'
#' In practical terms, \eqn{f(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\theta}})}
#' is better (worse) than \eqn{g(y_i \mid \boldsymbol{x}_i, \widehat{\boldsymbol{\gamma}})}
#' if \eqn{T > z_{\alpha/2}} (or \eqn{T < -z_{\alpha/2}}).
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string with the method used.}
#' \item{data.name}{a character string given the name of families models under comparison.}
#'
#' @author Ricardo Rasmussen Petterle
#'
#' @references
#' Vuong, Q. (1989). Likelihood ratio tests for model selection and
#' non-nested hypotheses. \emph{Econometrica}, \bold{57}(2), 307--333.
#'
#' @examples
#' \donttest{
#' # Simulate independent data
#' sim_data <- sim_unitreg(n_id = 500, family = vasicek(), formula_mu = ~ x, 
#'                         beta_mu = c(0.5, -0.5), beta_phi = 1.5,
#'                         cov_config = list(x = list(dist = "normal")))
#'
#' # Fit two competing models
#' fit_vasicek <- unitregTMB(Y ~ x, data = sim_data, family = vasicek())
#' fit_beta <- unitregTMB(Y ~ x, data = sim_data, family = beta_fam())
#'
#' # Perform Vuong test
#' ans <- vuong_test(fit_vasicek, fit_beta)
#' print(ans)
#' }
#'
#' @importFrom stats var pnorm
#' @rdname vuong_test_unitregTMB
#' @export
vuong_test.unitregTMB <- function(object1, object2, alternative = c("two.sided", "less", "greater"), ...) {
  
  alternative <- match.arg(alternative)
  
  if (!inherits(object2, "unitregTMB")) {
    stop("Both objects must be of class 'unitregTMB'.")
  }
  
  n1 <- object1$nobs
  n2 <- object2$nobs
  if (n1 != n2) {
    stop("Models must be fitted to the same number of observations.")
  }
  
  tau_1 <- object1$tau
  tau_2 <- object2$tau
  if (!is.null(tau_1) && !is.null(tau_2) && tau_1 != tau_2) {
    warning("Comparison is being performed between models fitted for different quantiles!", call. = FALSE)
  }
  
  if (object1$has_random_effects_mu || object2$has_random_effects_mu) {
    warning("One or both models contain random effects. The Vuong test utilizes conditional pointwise log-likelihoods, which may not properly reflect the marginal likelihood equivalence.", call. = FALSE)
  }
  
  get_pointwise_ll <- function(model) {
    rep <- model$obj$report()
    
    if (!is.null(rep$conditional_log_lik_unweighted)) {
      return(rep$conditional_log_lik_unweighted)
    }
    
    if (!is.null(rep$ll_i)) return(rep$ll_i)
    if (!is.null(rep$nll_i)) return(-rep$nll_i)
    
    stop(paste("Pointwise log-likelihoods not found in the TMB report for model", 
               if(is.list(model$family)) model$family$name else model$family, 
               ". Please ensure 'conditional_log_lik_unweighted' is exported via REPORT() in the C++ template."), 
         call. = FALSE)
  }
  
  ll_1 <- get_pointwise_ll(object1)
  ll_2 <- get_pointwise_ll(object2)
  
  om2 <- (n1 - 1) / n1 * stats::var(ll_1 - ll_2)
  lr <- sum(ll_1 - ll_2)
  
  if (om2 <= 1e-12) {
    warning("Models appear to be virtually identical; variance of log-likelihood ratio is near zero.", call. = FALSE)
    Tstat <- 0
  } else {
    Tstat <- (1 / sqrt(n1)) * lr / sqrt(om2)
  }
  names(Tstat) <- "T_LR"
  
  pvalue <- switch(alternative,
                   "two.sided" = 2 * stats::pnorm(abs(Tstat), lower.tail = FALSE),
                   "greater"   = stats::pnorm(Tstat, lower.tail = FALSE),
                   "less"      = stats::pnorm(Tstat, lower.tail = TRUE)
  )
  
  fam1 <- if(is.list(object1$family)) object1$family$name else object1$family
  fam2 <- if(is.list(object2$family)) object2$family$name else object2$family
  
  method <- paste0("Vuong likelihood ratio test for non-nested models (",
                   fam1, " versus ", fam2, ")")
  
  out <- list(statistic = Tstat, p.value = pvalue,
              alternative = alternative,
              method = method,
              data.name = paste0(fam1, " versus ", fam2))
  class(out) <- "htest"
  
  return(out)
}

#' @title Pairwise Vuong Test for Multiple unitregTMB Models
#'
#' @description Calculates pairwise comparisons between multiple fitted models, 
#' performing the Vuong test for objects of class \code{unitregTMB}.
#'
#' @param ... \code{unitregTMB} objects separated by commas.
#' @param lt A list containing one or more \code{unitregTMB} objects.
#' @param p.adjust.method A character string specifying the method for multiple
#' testing adjustment; almost always one of \code{p.adjust.methods}. Can be abbreviated.
#' @param alternative Indicates the alternative hypothesis and must be one
#' of \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}. Can be abbreviated.
#'
#' @return An object of class \code{"pairwise_htest"}
#'
#' @seealso \code{\link{vuong_test}}, \code{\link[stats]{p.adjust}}
#'
#' @examples
#' \donttest{
#' sim_data <- sim_unitreg(n_id = 500, family = vasicek(), formula_mu = ~ x, 
#'                         beta_mu = c(0.5, -0.5), beta_phi = 1.5,
#'                         cov_config = list(x = list(dist = "normal")))
#'
#' fit_vas <- unitregTMB(Y ~ x, data = sim_data, family = vasicek())
#' fit_bet <- unitregTMB(Y ~ x, data = sim_data, family = beta_fam())
#' fit_kum <- unitregTMB(Y ~ x, data = sim_data, family = kumaraswamy())
#'
#' ans <- pairwise_vuong_test(fit_vas, fit_bet, fit_kum)
#' print(ans)
#' }
#'
#' @importFrom stats pairwise.table p.adjust.methods
#' @rdname pairwise_vuong_test_unitregTMB
#' @export
pairwise_vuong_test.unitregTMB <- function(..., lt = NULL, p.adjust.method = stats::p.adjust.methods,
                                alternative = c("two.sided", "less", "greater")) {
  
  if (is.null(lt)) {
    lt <- list(...)
  }
  
  if (length(lt) < 2) {
    stop("At least two 'unitregTMB' models are required for pairwise comparison.", call. = FALSE)
  }
  
  if (!all(sapply(lt, inherits, "unitregTMB"))) {
    stop("All elements must be fitted objects of class 'unitregTMB'.", call. = FALSE)
  }
  
  p.adjust.method <- match.arg(p.adjust.method)
  alternative <- match.arg(alternative)
  
  families <- sapply(lt, function(x) {
    if(is.list(x$family)) x$family$name else as.character(x$family)
  })
  
  nfam <- length(families)
  
  compare.levels <- function(i, j) {
    xi <- lt[[as.integer(i)]]
    xj <- lt[[as.integer(j)]]
    vuong_test(object1 = xi, object2 = xj, alternative = alternative)$p.value
  }
  
  pval <- stats::pairwise.table(compare.levels = compare.levels,
                                level.names = seq_along(families),
                                p.adjust.method = p.adjust.method)
  
  colnames(pval) <- families[1L:(nfam - 1)]
  rownames(pval) <- families[2L:nfam]
  
  method <- "Vuong likelihood ratio test for non-nested models"
  dname <- paste(deparse(lt[[1L]]$call, width.cutoff = 60L), collapse = "\n")
  
  out <- list(method = method, data.name = dname, p.value = pval,
              p.adjust.method = p.adjust.method)
  class(out) <- "pairwise.htest"
  
  return(out)
}
