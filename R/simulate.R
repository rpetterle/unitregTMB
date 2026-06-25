#' Simulate Responses from a Fitted unitregTMB Model
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
#' @importFrom stats predict runif
#' @export
simulate.unitregTMB <- function(object, nsim = 1, seed = NULL, ...) {
  
  # Controlo da semente aleatória (Padrão ouro do R)
  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  
  # 1. Extrair os parâmetros preditos pelo modelo ajustado
  n <- object$nobs
  tau <- object$tau 
  family_code <- object$obj$env$data$family
  
  preds <- predict(object, type = "parameter")
  mu_hat  <- pmax(pmin(preds$mu, 1 - 1e-6), 1e-6)
  phi_hat <- preds$phi
  p0_hat  <- preds$p0
  p1_hat  <- preds$p1
  
  # 2. Inicializar a lista que guardará as simulações
  sim_list <- vector("list", nsim)
  
  # 3. Loop de Simulação (Matemática validada no HNP)
  for (i in seq_len(nsim)) {
    
    # A. Gera a parte contínua via Rcpp
    y_cont_sim <- get_random_continuous(family_code, n, mu_hat, phi_hat, tau)
    y_cont_sim <- pmax(pmin(y_cont_sim, 1 - 1e-6), 1e-6)
    
    # B. Aplica a Inflação (ZOI)
    y_sim <- numeric(n)
    for (k in 1:n) {
      cat_prob_vec <- c(p0_hat[k], p1_hat[k], 1 - p0_hat[k] - p1_hat[k])
      
      # Trava numérica de segurança para probabilidades
      if (any(is.na(cat_prob_vec)) || any(cat_prob_vec < 0) || sum(cat_prob_vec) < 0.99 || sum(cat_prob_vec) > 1.01) {
        cat_prob_vec[cat_prob_vec < 0] <- 0
        cat_prob_vec <- cat_prob_vec / sum(cat_prob_vec)
      }
      
      category <- sample(0:2, size = 1, prob = cat_prob_vec)
      y_sim[k] <- if (category == 0) 0 else if (category == 1) 1 else y_cont_sim[k]
    }
    
    sim_list[[i]] <- y_sim
  }
  
  # 4. Retorna no formato clássico do simulate() no R
  out <- as.data.frame(sim_list)
  colnames(out) <- paste0("sim_", seq_len(nsim))
  
  attr(out, "seed") <- seed
  class(out) <- c("data.frame")
  
  return(out)
}