#' @title Predict Method for unitregTMB Fits
#' @description Obtains predictions from a fitted \code{unitregTMB} object.
#' 
#' @param object A fitted \code{unitregTMB} model.
#' @param newdata Optional data frame for which to compute predictions. If omitted, the fitted values are used.
#' @param type Character indicating the type of prediction. 
#'   \code{"response"} (default) returns predictions on the scale of the response variable.
#'   \code{"link"} returns predictions on the scale of the linear predictor for the location parameter (mu).
#'   \code{"parameter"} returns a data.frame with all predicted distribution parameters (\code{mu}, \code{phi}, \code{p0}, \code{p1}).
#' @param re.form (Currently ignored; included for compatibility).
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return A vector or data.frame of predictions.
#' 
#' @examples
#' \donttest{
#' # Assuming 'da' is your dataset:
#' # fit <- unitregTMB(Y ~ educ + refill, data = da, family = vasicek())
#' #
#' # # Predictions on the original dataset
#' # head(predict(fit, type = "response"))
#' # head(predict(fit, type = "parameter"))
#' #
#' # # Predictions on new data
#' # new_data <- data.frame(educ = factor("HS+"), refill = factor("2"))
#' # predict(fit, newdata = new_data, type = "response")
#' }
#' 
#' @importFrom stats model.frame model.matrix model.offset na.pass terms delete.response plogis
#' @export
predict.unitregTMB <- function(object, newdata = NULL, 
                               type = c("response", "link", "parameter"), 
                               re.form = NULL, ...) {
  
  type <- match.arg(type)
  
  if (is.null(newdata)) {
    mu_val  <- fitted(object, type = "mu")
    phi_val <- fitted(object, type = "phi")
    p0_val  <- fitted(object, type = "p0")
    p1_val  <- fitted(object, type = "p1")
    
    if (type == "link") {
      link_obj <- object$link.mu
      return(link_obj$linkfun(mu_val))
    } 
    
    if (type == "parameter") {
      return(data.frame(mu = mu_val, phi = phi_val, p0 = p0_val, p1 = p1_val))
    }
    
    expected_response <- p1_val + (1 - p0_val - p1_val) * mu_val
    return(expected_response)
  }
  
  fixed_form <- reformulas::nobars(object$formula)
  
  tt_mu <- delete.response(terms(fixed_form))
  mf_mu <- model.frame(tt_mu, newdata, na.action = na.pass, xlev = object$terms)
  X_mu  <- model.matrix(tt_mu, mf_mu)
  
  offset_new <- stats::model.offset(mf_mu)
  if (is.null(offset_new)) {
    offset_new <- rep(0, nrow(X_mu))
  }
  
  tt_phi <- terms(object$phi.formula)
  mf_phi <- model.frame(tt_phi, newdata, na.action = na.pass)
  X_phi  <- model.matrix(tt_phi, mf_phi)
  
  eta_mu  <- as.vector(X_mu %*% object$model_coef$mu) + offset_new
  eta_phi <- as.vector(X_phi %*% object$model_coef$phi)
  
  if (object$has_random_effects_mu) {
    tryCatch({
      re_list <- reformulas::findbars(object$formula)
      Z_new_list <- list()
      
      for (i in seq_along(re_list)) {
        bar <- re_list[[i]]
        
        rhs_expr <- bar[[3]]
        grp_vars <- all.vars(rhs_expr)
        orig_levels <- object$re_info_mu$group_levels[[i]]
        
        if (length(grp_vars) > 1) {
          grp_factor <- interaction(newdata[grp_vars], drop = FALSE)
        } else {
          grp_factor <- newdata[[grp_vars]]
        }
        
        grp_factor <- factor(grp_factor, levels = orig_levels)
        
        lhs_expr <- bar[[2]]
        lhs_formula <- stats::as.formula(paste("~", paste(deparse(lhs_expr), collapse = " ")))
        
        mf_re <- stats::model.frame(lhs_formula, newdata, na.action = stats::na.pass)
        X_re <- stats::model.matrix(lhs_formula, mf_re)
        
        A <- Matrix::fac2sparse(grp_factor) 
        B <- Matrix::t(X_re)                
        
        Zt_i <- Matrix::KhatriRao(A, B, make.dimnames = FALSE)
        Z_new_list[[i]] <- Matrix::t(Zt_i)
      }
      
      Z_new <- do.call(cbind, Z_new_list)
      
      par_full <- object$obj$env$parList(object$opt$par)
      u_vec <- par_full$u_mu
      
      if (ncol(Z_new) == length(u_vec)) {
        eta_mu <- eta_mu + as.vector(Z_new %*% u_vec)
      } else {
        warning("Dimension mismatch in random effects for 'newdata'. Predictions will be at the population level (fixed effects only).")
      }
    }, error = function(e) {
      warning("Could not construct random effects for 'newdata'. Predictions will be at the population level.")
    })
  }
  
  n_new <- nrow(newdata)
  eta_p0 <- rep(0, n_new)
  eta_p1 <- rep(0, n_new)
  
  if (object$has_p0) {
    tt_p0 <- terms(object$p0.formula)
    mf_p0 <- model.frame(tt_p0, newdata, na.action = na.pass)
    X_p0  <- model.matrix(tt_p0, mf_p0)
    eta_p0 <- as.vector(X_p0 %*% object$model_coef$p0)
  }
  
  if (object$has_p1) {
    tt_p1 <- terms(object$p1.formula)
    mf_p1 <- model.frame(tt_p1, newdata, na.action = na.pass)
    X_p1  <- model.matrix(tt_p1, mf_p1)
    eta_p1 <- as.vector(X_p1 %*% object$model_coef$p1)
  }
  
  if (object$has_p0 && object$has_p1) {
    denom <- 1 + exp(eta_p0) + exp(eta_p1)
    p0_val <- exp(eta_p0) / denom
    p1_val <- exp(eta_p1) / denom
  } else if (object$has_p0) {
    p0_val <- stats::plogis(eta_p0)
    p1_val <- rep(0, n_new)
  } else if (object$has_p1) {
    p0_val <- rep(0, n_new)
    p1_val <- stats::plogis(eta_p1)
  } else {
    p0_val <- rep(0, n_new)
    p1_val <- rep(0, n_new)
  }
  
  link_obj <- object$link.mu
  mu_val   <- link_obj$linkinv(eta_mu)

  linkobj_phi <- stats::make.link(object$family_object$phi_link_r_name)
  phi_val <- linkobj_phi$linkinv(eta_phi)
  
  if (type == "link") return(eta_mu)
  
  if (type == "parameter") {
    return(data.frame(mu = mu_val, phi = phi_val, p0 = p0_val, p1 = p1_val))
  }
  
  expected_response <- p1_val + (1 - p0_val - p1_val) * mu_val
  return(expected_response)
}
