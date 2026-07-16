#' @name write_methods.unitregTMB
#' 
#' @title Generate Methods Section Text for unitregTMB
#' 
#' @description 
#' Generates a boilerplate "Methods" section text for academic papers, describing 
#' the statistical model fitted with \code{unitregTMB}. It dynamically detects 
#' the distribution family, link function, presence of zero/one inflation, 
#' and random effects, writing the text in either English or Portuguese.
#' 
#' @param object A fitted \code{unitregTMB} model.
#' @param language Character string indicating the language for the generated text. 
#'        Options are \code{"en"} for English (default) or \code{"pt"} for Portuguese.
#' @param ... Additional arguments passed to methods.
#' 
#' @return A list (invisibly) containing two elements: \code{text} (the generated 
#'    methods paragraph) and \code{bibtex} (the citation string). The text is also 
#'    printed to the console.
#' 
#' @examples
#' \donttest{
#' # Assuming 'fit' is a fitted unitregTMB model:
#' # Generate text in English
#' # write_methods(fit, language = "en")
#' 
#' # Generate text in Portuguese and save to an object
#' # methods_output <- write_methods(fit, language = "pt")
#' # methods_output$text
#' }
#' 
#' @rdname write_methods_unitregTMB
#' @export
write_methods.unitregTMB <- function(object, language = c("en", "pt"), ...) {
  lang <- match.arg(language)
  
  if (inherits(object$family, "unitregTMBFamily") || is.list(object$family)) {
    f_name <- object$family$name
    l_name <- object$family$link_r_name
    link_text_en <- sprintf("a %s link function", l_name)
    link_text_pt <- sprintf("função de ligação %s", l_name)
    if(!is.null(object$family$tau)) f_name <- paste0(f_name, " ($\\tau = ", object$family$tau, "$)")
  } else {
    f_name <- if (!is.null(object$family)) as.character(object$family) else "Unknown"
    link_text_en <- "its specified link function"
    link_text_pt <- "sua função de ligação especificada"
  }
  
  has_re <- isTRUE(object$has_random_effects_mu) || isTRUE(object$obj$env$data$has_random_effects_mu == 1)
  
  has_p0 <- isTRUE(object$has_p0_inflation) || isTRUE(object$obj$env$data$has_p0_inflation == 1)
  has_p1 <- isTRUE(object$has_p1_inflation) || isTRUE(object$obj$env$data$has_p1_inflation == 1)
  
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
  
  wrapped_text <- paste(strwrap(text, width = 80), collapse = "\n")
  
  citation_bib <- "@Manual{unitregTMB2026,
  title = {unitregTMB: Regression Models for Correlated Continuous Bounded Data},
  author = {Ricardo Rasmussen Petterle},
  year = {2026},
  note = {R package version 0.1.0},
  url = {https://github.com/rpetterle/unitregTMB}
}"
  
  cat("\n------------------------------------------------------------\nSUGGESTED METHODS TEXT:\n------------------------------------------------------------\n")
  cat(wrapped_text, "\n")
  cat("\n------------------------------------------------------------\nCITATION:\n------------------------------------------------------------\n")
  cat(citation_bib, "\n")
  
  invisible(list(text = text, bibtex = citation_bib))
}
