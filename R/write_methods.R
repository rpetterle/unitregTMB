#' Generate Methods Section Text for unitregTMB
#' @export
write_methods <- function(model, language = c("en", "pt")) {
  lang <- match.arg(language)
  
  # --- 1. Lógica da Família (Respeitando family$name) ---
  if (inherits(model$family, "unitregTMBFamily") || is.list(model$family)) {
    f_name <- model$family$name
    l_name <- model$family$link_r_name
    link_text_en <- sprintf("a %s link function", l_name)
    link_text_pt <- sprintf("função de ligação %s", l_name)
    if(!is.null(model$family$tau)) f_name <- paste0(f_name, " ($\\tau = ", model$family$tau, "$)")
  } else {
    # Fallback seguro para quando family é apenas um texto (ex: "Beta (mean)")
    f_name <- if (!is.null(model$family)) as.character(model$family) else "Unknown"
    link_text_en <- "its specified link function"
    link_text_pt <- "sua função de ligação especificada"
  }
  
  # --- 2. Checagem de Efeitos Aleatórios ---
  # Garante que lê tanto se estiver na raiz do modelo quanto dentro dos dados do TMB
  has_re <- isTRUE(model$has_random_effects_mu) || isTRUE(model$obj$env$data$has_random_effects_mu == 1)
  
  # --- 3. Checagem de Inflações ---
  has_p0 <- isTRUE(model$has_p0_inflation) || isTRUE(model$obj$env$data$has_p0_inflation == 1)
  has_p1 <- isTRUE(model$has_p1_inflation) || isTRUE(model$obj$env$data$has_p1_inflation == 1)
  
  zoi_text_en <- ""
  zoi_text_pt <- ""
  if (has_p0 && has_p1) {
    zoi_text_en <- " To account for boundary observations, a zero-and-one-inflated mixture structure was incorporated."
    zoi_text_pt <- " Para considerar as observações nos limites, uma estrutura de mistura inflacionada de zeros e uns foi incorporada."
  } else if (has_p0) {
    zoi_text_en <- " To account for exact zeros, a zero-inflated mixture structure was incorporated."
    zoi_text_pt <- " Para considerar zeros exatos, uma estrutura de mistura inflacionada de zeros foi incorporada."
  } else if (has_p1) {
    zoi_text_en <- " To account for exact ones, a one-inflated mixture structure was incorporated."
    zoi_text_pt <- " Para considerar uns exatos, uma estrutura de mistura inflacionada de uns foi incorporada."
  }
  
  # --- 4. Construção do Texto Baseado no Tipo de Modelo ---
  if (lang == "en") {
    
    if (has_re) {
      framework <- "mixed-effects regression framework"
      re_text <- " Unobserved heterogeneity was captured via latent random effects."
      est_method <- "maximum marginal likelihood method via the Laplace approximation"
    } else {
      framework <- "generalized regression framework"
      re_text <- ""
      est_method <- "maximum likelihood method"
    }
    
    text <- sprintf(
      "Data were analyzed using a %s specifically designed for continuous bounded data. The response variable was modeled assuming a %s distribution with %s for the location parameter.%s%s Parameter estimation was performed using the %s, implemented in the Template Model Builder (TMB) automatic differentiation engine. All analyses were conducted in the R statistical environment using the 'unitregTMB' package (Petterle, 2026).",
      framework, f_name, link_text_en, zoi_text_en, re_text, est_method
    )
    
  } else {
    
    if (has_re) {
      framework <- "modelo de regressão de efeitos mistos"
      re_text <- " A heterogeneidade não observada foi capturada via efeitos aleatórios latentes."
      est_method <- "máxima verossimilhança marginal via aproximação de Laplace"
    } else {
      framework <- "modelo de regressão generalizado"
      re_text <- ""
      est_method <- "máxima verossimilhança"
    }
    
    text <- sprintf(
      "Os dados foram analisados utilizando um %s especificamente desenvolvido para dados contínuos limitados. A variável resposta foi modelada assumindo uma distribuição %s com %s para o parâmetro de locação.%s%s A estimação dos parâmetros foi realizada pelo método de %s, implementado através do motor de diferenciação automática Template Model Builder (TMB). Todas as análises foram conduzidas no ambiente R utilizando o pacote 'unitregTMB' (Petterle, 2026).",
      framework, f_name, link_text_pt, zoi_text_pt, re_text, est_method
    )
  }
  
  # --- 5. Impressão ---
  cat("\n------------------------------------------------------------\nSUGGESTED METHODS TEXT:\n------------------------------------------------------------\n")
  cat(strwrap(text, width = 80), sep = "\n")
  cat("\n\n------------------------------------------------------------\nCITATION:\n------------------------------------------------------------\n")
  cat("@Manual{unitregTMB2026,\n  title = {unitregTMB: A Unified Framework for Bounded Data Regression},\n  author = {Ricardo R. Petterle},\n  year = {2026},\n  note = {R package version 0.1.0}\n}\n")
}