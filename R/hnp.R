#' @export
hnp <- function(object, ...) {
  UseMethod("hnp")
}

#' Half-Normal Plots with Simulated Envelopes for unitregTMB
#' 
#' @description 
#' Produces a (half-)normal probability plot with a simulated envelope for a fitted
#' \code{unitregTMB} model. This diagnostic tool assesses the goodness-of-fit of the 
#' chosen continuous distribution as well as the zero-one inflation (ZOI) components.
#' 
#' @param object A fitted \code{unitregTMB} model.
#' @param nsim Number of simulations used to compute the envelope. Default is 99.
#' @param halfnormal logical. If \code{TRUE} (default), a half-normal plot is produced. 
#'        If \code{FALSE}, a normal plot is produced.
#' @param plot logical. Should the plot be printed? Default is \code{TRUE}.
#' @param level Confidence level of the simulated envelope. Default is 0.95.
#' @param resid.type Type of residuals to be used: \code{"quantile"} (default) or \code{"cox-snell"}.
#' @param ncpus Number of cores to use for parallel computing. Default is 1.
#' @param ... Additional arguments passed to the \code{plot} function.
#' 
#' @return A list (invisibly) containing the observed residuals (\code{res_obs}) and a matrix 
#'         with the lower and upper bounds of the simulated envelope (\code{envelope}).
#' 
#' @importFrom stats qnorm qexp quantile median plogis
#' @importFrom graphics plot lines points grid
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom parallel makeCluster stopCluster parLapply clusterExport clusterEvalQ
#' @export
hnp.unitregTMB <- function(object, nsim = 99, halfnormal = TRUE, plot = TRUE, 
                           level = 0.95, resid.type = c("quantile", "cox-snell"), 
                           ncpus = 1, ...) {
  
  resid.type <- match.arg(resid.type)
  
  # 1. Preparation
  n <- object$nobs
  tau <- object$tau 
  family_code <- object$obj$env$data$family
  
  preds <- predict(object, type = "parameter")
  mu_hat  <- pmax(pmin(preds$mu, 1 - 1e-6), 1e-6)
  phi_hat <- preds$phi
  p0_hat  <- preds$p0
  p1_hat  <- preds$p1
  
  res_sim <- matrix(NA, nrow = n, ncol = nsim)
  
  # 2. Ultra-Fast Refit Settings
  ctrl_fast <- object$control
  ctrl_fast$sd.report <- FALSE 
  ctrl_fast$get.joint.precision <- FALSE
  ctrl_fast$verbose <- FALSE
  
  dt_tmb_base <- list(
    Y = as.numeric(object$obj$env$data$Y), 
    y_class = as.integer(object$obj$env$data$y_class), 
    X_mu = object$obj$env$data$X_mu, 
    X_phi = object$obj$env$data$X_phi, 
    X_p0 = object$obj$env$data$X_p0, 
    X_p1 = object$obj$env$data$X_p1, 
    Z_mu = object$obj$env$data$Z_mu,
    weights = object$obj$env$data$weights,      
    offset_mu = object$obj$env$data$offset_mu,
    family = object$obj$env$data$family, 
    link = object$obj$env$data$link,
    tau = object$obj$env$data$tau,
    has_p0_inflation = object$obj$env$data$has_p0_inflation, 
    has_p1_inflation = object$obj$env$data$has_p1_inflation,
    has_random_effects_mu = object$obj$env$data$has_random_effects_mu,
    re_term_n_components = object$obj$env$data$re_term_n_components, 
    re_term_chol_starts = object$obj$env$data$re_term_chol_starts,
    u_mu_term_starts = object$obj$env$data$u_mu_term_starts, 
    n_re_levels_list_data = object$obj$env$data$n_re_levels_list_data
  )
  
  raw_pars <- object$opt$par
  pars_init_base <- list(
    beta_mu = as.numeric(raw_pars[names(raw_pars) == "beta_mu"]),
    beta_phi = as.numeric(raw_pars[names(raw_pars) == "beta_phi"]),
    beta_p0 = as.numeric(raw_pars[names(raw_pars) == "beta_p0"]),
    beta_p1 = as.numeric(raw_pars[names(raw_pars) == "beta_p1"]),
    u_mu = if(object$has_random_effects_mu) rep(0, ncol(dt_tmb_base$Z_mu)) else numeric(0),
    log_chol_re_mu_combined = as.numeric(raw_pars[names(raw_pars) == "log_chol_re_mu_combined"])
  )
  
  has_re <- object$has_random_effects_mu
  random_args <- if(has_re) "u_mu" else NULL
  
  # ---------------------------------------------------------------------------
  # Função Trabalhadora 
  # ---------------------------------------------------------------------------
  simulate_and_refit <- function(i) {
    dt_tmb <- dt_tmb_base
    pars_init <- pars_init_base
    
    y_cont_sim <- get_random_continuous(family_code, n, mu_hat, phi_hat, tau)
    y_cont_sim <- pmax(pmin(y_cont_sim, 1 - 1e-6), 1e-6)
    
    y_sim <- rep(NA, n)
    for(k in 1:n) {
      cat_prob_vec <- c(p0_hat[k], p1_hat[k], 1 - p0_hat[k] - p1_hat[k])
      if (any(is.na(cat_prob_vec)) || any(cat_prob_vec < 0) || sum(cat_prob_vec) < 0.99 || sum(cat_prob_vec) > 1.01) {
        cat_prob_vec[cat_prob_vec < 0] <- 0
        cat_prob_vec <- cat_prob_vec / sum(cat_prob_vec)
      }
      category <- sample(0:2, size = 1, prob = cat_prob_vec)
      y_sim[k] <- if(category == 0) 0 else if (category == 1) 1 else y_cont_sim[k]
    }
    
    dt_tmb$Y <- y_sim
    y_class_sim <- rep(2L, n)
    if (object$has_p0) y_class_sim[y_sim < 1e-12] <- 0L
    if (object$has_p1) y_class_sim[y_sim > 1 - 1e-12] <- 1L
    dt_tmb$y_class <- as.integer(y_class_sim)
    
    fit_sim <- try({
      unitregTMB.fit(data_tmb = dt_tmb, parameters = pars_init, 
                     random_args = random_args, control = ctrl_fast, silent = TRUE)
    }, silent = TRUE)
    
    if(!inherits(fit_sim, "try-error") && !is.null(fit_sim) && fit_sim$convergence == 0) {
      sim_par <- fit_sim$obj$env$parList(fit_sim$opt$par)
      
      eta_mu <- as.vector(dt_tmb$X_mu %*% sim_par$beta_mu)
      if (!is.null(dt_tmb$offset_mu)) eta_mu <- eta_mu + dt_tmb$offset_mu
      if(has_re) eta_mu <- eta_mu + as.vector(dt_tmb$Z_mu %*% sim_par$u_mu)
      
      sim_mu  <- object$link.mu$linkinv(eta_mu)
      sim_phi <- exp(as.vector(dt_tmb$X_phi %*% sim_par$beta_phi))
      
      eta_p0_sim <- if(object$has_p0) as.vector(dt_tmb$X_p0 %*% sim_par$beta_p0) else rep(0, n)
      eta_p1_sim <- if(object$has_p1) as.vector(dt_tmb$X_p1 %*% sim_par$beta_p1) else rep(0, n)
      
      if (object$has_p0 && object$has_p1) {
        denom_sim <- 1 + exp(eta_p0_sim) + exp(eta_p1_sim)
        sim_p0 <- exp(eta_p0_sim) / denom_sim
        sim_p1 <- exp(eta_p1_sim) / denom_sim
      } else if (object$has_p0) {
        sim_p0 <- stats::plogis(eta_p0_sim); sim_p1 <- rep(0, n)
      } else if (object$has_p1) {
        sim_p0 <- rep(0, n); sim_p1 <- stats::plogis(eta_p1_sim)
      } else {
        sim_p0 <- rep(0, n); sim_p1 <- rep(0, n)
      }
      
      cdf_vals <- get_cdf_continuous(family_code, y_sim, sim_mu, sim_phi, tau)
      if (object$has_p0 || object$has_p1) cdf_vals <- correct_for_inflation(y_sim, cdf_vals, sim_p0, sim_p1)
      
      return(compute_residuals(cdf_vals, resid.type))
    }
    return(rep(NA, n))
  }
  
  # ---------------------------------------------------------------------------
  # 3. Execução
  # ---------------------------------------------------------------------------
  if (ncpus > 1) {
    cat(sprintf("Simulating envelope with %d core(s)...\n", ncpus))
    cl <- parallel::makeCluster(ncpus)
    parallel::clusterEvalQ(cl, { try(library(unitregTMB), silent = TRUE) })
    parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv))
    parallel::clusterExport(cl, varlist = c("dt_tmb_base", "pars_init_base", "simulate_and_refit",
                                            "get_random_continuous", "correct_for_inflation", 
                                            "compute_residuals", "unitregTMB.fit", "object", 
                                            "family_code", "n", "mu_hat", "phi_hat", "tau", 
                                            "p0_hat", "p1_hat", "random_args", "ctrl_fast", "resid.type"), 
                            envir = environment())
    res_list <- parallel::parLapply(cl, 1:nsim, simulate_and_refit)
    parallel::stopCluster(cl)
    
    for (i in 1:nsim) res_sim[, i] <- res_list[[i]]
    
  } else {
    cat("Simulating envelope...\n")
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
    for(i in 1:nsim) {
      res_sim[, i] <- simulate_and_refit(i)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  
  # 4. Envelope Processing
  failed <- apply(res_sim, 2, function(x) any(is.na(x)))
  if(any(failed)) {
    warning(paste(sum(failed), "simulations failed or diverged. Envelope calculated with remaining valid fits."))
    res_sim <- res_sim[, !failed, drop = FALSE]
  }
  if(ncol(res_sim) < 2) stop("Error: Almost all simulations failed. Check your model and data.")
  
  res_obs_raw <- residuals(object, type = resid.type)
  
  # Sorting Logic
  if (resid.type == "cox-snell") {
    res_obs <- sort(res_obs_raw)
    res_sim <- apply(res_sim, 2, sort, na.last = NA)
    n_plot <- length(res_obs)
    res_teo <- stats::qexp((1:n_plot - 0.375) / (n_plot + 0.25))
    dist_lab <- "Exponential"; ylab_txt <- "Ordered absolute Cox-Snell residuals"
    
  } else if (halfnormal) {
    res_obs <- sort(abs(res_obs_raw))
    res_sim <- apply(res_sim, 2, function(x) sort(abs(x), na.last = NA))
    n_plot <- length(res_obs)
    res_teo <- stats::qnorm((1:n_plot + n_plot - 0.125) / (2 * n_plot + 0.5))
    dist_lab <- "Half-Normal"; ylab_txt <- "Randomized quantile residuals"
    
  } else {
    res_obs <- sort(res_obs_raw)
    res_sim <- apply(res_sim, 2, sort, na.last = NA)
    n_plot <- length(res_obs)
    res_teo <- stats::qnorm((1:n_plot - 0.375) / (n_plot + 0.25))
    dist_lab <- "Normal"; ylab_txt <- "Randomized quantile residuals"
  }
  
  alpha <- (1 - level) / 2
  res_lwr <- apply(res_sim, 1, stats::quantile, probs = alpha, na.rm = TRUE)
  res_upr <- apply(res_sim, 1, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
  res_mid <- apply(res_sim, 1, stats::median, na.rm = TRUE)
  
  # =========================================================================
  # 5. Graphical Plotting (Ggplot2 "theme_bw" style using Base R)
  # =========================================================================
  if (plot) {
    Ry <- range(c(res_lwr, res_upr, res_obs), na.rm = TRUE)
    Rx <- range(res_teo, na.rm = TRUE)
    
    # Base vazia com títulos e painel limpo
    graphics::plot(res_teo, res_obs, ylim = Ry, xlim = Rx, type = "n",
                   xlab = "Theoretical quantiles", 
                   ylab = ylab_txt, 
                   main = object$family, 
                   bty = "o", 
                   panel.first = {
                     # Adiciona a grelha cinza claro antes das linhas
                     graphics::grid(col = "gray92", lty = 1, lwd = 1.2)
                   }, ...)
    
    # Limites do Envelope (Linhas Azuis Sólidas)
    graphics::lines(res_teo, res_lwr, lty = 1, col = "#0066CC", lwd = 1.8)
    graphics::lines(res_teo, res_upr, lty = 1, col = "#0066CC", lwd = 1.8)
    
    # Linha Mediana (Preta Tracejada)
    graphics::lines(res_teo, res_mid, lty = 2, col = "black", lwd = 1.5)
    
    # Pontos Observados (Cruzes Pretas)
    graphics::points(res_teo, res_obs, pch = 3, col = "black", cex = 0.8)
  } 
  
  invisible(list(res_obs = res_obs, envelope = cbind(Lower = res_lwr, Median = res_mid, Upper = res_upr)))
}
