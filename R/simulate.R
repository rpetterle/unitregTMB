#' @title Simulate Responses from a Fitted unitregTMB Model
#' 
#' @description 
#' Generates simulated response data from a fitted \code{unitregTMB} model. 
#' This function perfectly mimics the data-generating process (DGP) estimated by the model, 
#' including the continuous distribution and any zero/one inflation parameters.
#' 
#' @param object A fitted \code{unitregTMB} object.
#' @param nsim Number of response vectors to simulate. Defaults to 1.
#' @param seed An object specifying if and how the random number generator should be initialized.
#' @param ... Additional optional arguments.
#' 
#' @return A data.frame with \code{nsim} columns containing the simulated responses.
#' 
#' @examples
#' \donttest{
#' # Assuming 'da' is your dataset:
#' # fit <- unitregTMB(Y ~ educ, p0.formula = ~ 1, data = da, family = vasicek())
#' # 
#' # # Simulate 5 new response vectors
#' # sim_data <- simulate(fit, nsim = 5, seed = 123)
#' # head(sim_data)
#' }
#' 
#' @importFrom stats predict runif
#' @export
simulate.unitregTMB <- function(object, nsim = 1, seed = NULL, ...) {
  
  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  
  n <- object$nobs
  tau <- object$tau 
  family_code <- object$obj$env$data$family
  
  preds <- predict(object, type = "parameter")
  mu_hat  <- pmax(pmin(preds$mu, 1 - 1e-6), 1e-6)
  phi_hat <- preds$phi
  p0_hat  <- preds$p0
  p1_hat  <- preds$p1
  
  # --- Vectorized probability normalization (Safety Check) ---
  prob_mat <- cbind(p0_hat, p1_hat, 1 - p0_hat - p1_hat)
  prob_mat[prob_mat < 0] <- 0
  prob_mat <- prob_mat / rowSums(prob_mat)
  
  cum_p0 <- prob_mat[, 1]
  cum_p1 <- cum_p0 + prob_mat[, 2]
  
  sim_list <- vector("list", nsim)
  
  for (i in seq_len(nsim)) {
    
    y_cont_sim <- get_random_continuous(family_code, n, mu_hat, phi_hat, tau)
    y_cont_sim <- pmax(pmin(y_cont_sim, 1 - 1e-6), 1e-6)
    
    # --- Vectorized Sampling for Zero-One Inflation ---
    rand_draws <- stats::runif(n)
    y_sim <- y_cont_sim
    
    y_sim[rand_draws < cum_p0] <- 0
    y_sim[rand_draws >= cum_p0 & rand_draws < cum_p1] <- 1
    
    sim_list[[i]] <- y_sim
  }
  
  out <- as.data.frame(sim_list)
  colnames(out) <- paste0("sim_", seq_len(nsim))
  
  attr(out, "seed") <- seed
  class(out) <- "data.frame"
  
  return(out)
}
