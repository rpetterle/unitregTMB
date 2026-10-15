#' Generate Methods-Section Text for a unitregTMB Model
#'
#' @description
#' Generates a concise Methods-section paragraph for an academic manuscript from
#' a fitted \code{unitregTMB} model. The generated text reflects the fitted
#' distribution, mean/mode/quantile parameterization, links, zero/one
#' augmentation, ordinary random effects, proper CAR spatial effects, and
#' estimation method.
#'
#' @param object A fitted \code{unitregTMB} object.
#' @param language Output language: \code{"en"} for English or \code{"pt"} for
#'   Portuguese.
#' @param ... Additional arguments, currently unused.
#'
#' @details
#' When both exact-zero and exact-one components are active, the generated text
#' describes the joint three-category softmax parameterization used by
#' \code{unitregTMB}. When only one boundary component is active, a binary
#' logistic submodel is described.
#'
#' Models containing ordinary or spatial latent effects are described as being
#' estimated by maximum marginal likelihood using the Laplace approximation.
#' Models without latent effects are described as maximum-likelihood fits.
#'
#' The returned text is intended as an editable starting point for a manuscript,
#' not as a substitute for reporting the complete model specification.
#'
#' @return Invisibly returns a list containing \code{text}, \code{bibtex},
#'   \code{parameterization}, \code{random_effect_components}, and
#'   \code{spatial_components}. The text and BibTeX entry are also printed.
#'
#' @examples
#' \donttest{
#' data("wqi", package = "unitregTMB")
#'
#' fit <- unitregTMB(y ~ location + quarter + (1 | id),
#'                   family = unitgamma(model_for = "mean"),
#'                   data = wqi)
#'
#' out <- write_methods(fit, language = "en")
#' out$text
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom regcore write_methods
#' @importFrom utils packageVersion
#' @method write_methods unitregTMB
#' @rdname write_methods_unitregTMB
#' @export
write_methods.unitregTMB <- function(object, language = c("en", "pt"), ...) {

  if (!inherits(object, "unitregTMB")) {
    stop("'object' must inherit from 'unitregTMB'.", call. = FALSE)
  }

  lang <- match.arg(language)
  components <- c("mu", "phi", "p0", "p1")

  value_or <- function(x, y) {
    if (!is.null(x) && length(x) > 0L && !all(is.na(x))) x else y
  }

  object_flag <- function(name) {
    x <- object[[name]]
    !is.null(x) && length(x) > 0L && isTRUE(as.logical(x[1L]))
  }

  has_p0 <- object_flag("has_p0") || object_flag("has_p0_inflation")
  has_p1 <- object_flag("has_p1") || object_flag("has_p1_inflation")

  re_flags <- setNames(
    vapply(
      components,
      function(comp) object_flag(paste0("has_random_effects_", comp)),
      logical(1)
    ),
    components
  )

  spatial_flags <- setNames(
    vapply(
      components,
      function(comp) object_flag(paste0("has_spatial_effects_", comp)),
      logical(1)
    ),
    components
  )

  has_re <- any(re_flags)
  has_spatial <- any(spatial_flags)
  has_latent <- has_re || has_spatial

  family_object <- if (is.list(object$family_object)) {
    object$family_object
  } else {
    list()
  }

  family_name <- value_or(
    family_object$name,
    if (!is.null(object$family) && length(object$family)) {
      as.character(object$family)[1L]
    } else {
      "specified"
    }
  )

  model_for <- value_or(family_object$model_for, NULL)

  if (!is.null(model_for)) {
    model_for <- as.character(model_for)[1L]
  }

  if (is.null(model_for) ||
      !model_for %in% c("mean", "mode", "quantile")) {
    code <- family_object$family_code

    if (!is.null(code) && length(code)) {
      code <- as.integer(code)[1L]
      if (!is.na(code) && code %in% 0:4) model_for <- "mean"
      if (!is.na(code) && code %in% 5:8) model_for <- "mode"
      if (!is.na(code) && code %in% 9:15) model_for <- "quantile"
    }
  }

  if (is.null(model_for) ||
      !model_for %in% c("mean", "mode", "quantile")) {
    model_for <- "location"
  }

  tau <- value_or(object$tau, family_object$tau)
  tau_num <- if (!is.null(tau) && length(tau)) {
    suppressWarnings(as.numeric(tau)[1L])
  } else {
    NA_real_
  }

  mu_link <- if (!is.null(object$link.mu) &&
                 is.list(object$link.mu) &&
                 !is.null(object$link.mu$name)) {
    as.character(object$link.mu$name)[1L]
  } else {
    as.character(value_or(family_object$link_r_name, "specified"))[1L]
  }

  phi_link <- as.character(
    value_or(family_object$phi_link_r_name, "log")
  )[1L]

  parameter_text_en <- switch(
    model_for,
    mean = "conditional mean of the continuous component",
    mode = "conditional mode of the continuous component",
    quantile = if (is.finite(tau_num)) {
      paste0(
        "conditional ", format(tau_num, trim = TRUE),
        "-quantile of the continuous component"
      )
    } else {
      "conditional quantile of the continuous component"
    },
    "location parameter of the continuous component"
  )

  parameter_text_pt <- switch(
    model_for,
    mean = "m\u00e9dia condicional do componente cont\u00ednuo",
    mode = "moda condicional do componente cont\u00ednuo",
    quantile = if (is.finite(tau_num)) {
      paste0(
        "quantil condicional de ordem ",
        format(tau_num, trim = TRUE),
        " do componente cont\u00ednuo"
      )
    } else {
      "quantil condicional do componente cont\u00ednuo"
    },
    "par\u00e2metro de localiza\u00e7\u00e3o do componente cont\u00ednuo"
  )

  labels_en <- c(
    mu = paste0(parameter_text_en, " (mu)"),
    phi = "second distributional parameter (phi)",
    p0 = "zero-mass probability (p0)",
    p1 = "one-mass probability (p1)"
  )

  labels_pt <- c(
    mu = paste0(parameter_text_pt, " (mu)"),
    phi = "segundo par\u00e2metro distribucional (phi)",
    p0 = "probabilidade de massa em zero (p0)",
    p1 = "probabilidade de massa em um (p1)"
  )

  join_labels <- function(x, language) {
    x <- unname(x)

    if (!length(x)) return("")
    if (length(x) == 1L) return(x)

    connector <- if (language == "en") " and " else " e "

    if (length(x) == 2L) {
      return(paste0(x[1L], connector, x[2L]))
    }

    paste0(
      paste(x[-length(x)], collapse = ", "),
      ",",
      connector,
      x[length(x)]
    )
  }

  if (has_p0 && has_p1) {
    augmentation_en <- paste0(
      " Exact-zero and exact-one probabilities were modeled jointly through ",
      "a three-category softmax structure with the continuous component as ",
      "the reference category."
    )
    augmentation_pt <- paste0(
      " As probabilidades de zeros e uns exatos foram modeladas conjuntamente ",
      "por uma estrutura softmax de tr\u00eas categorias, tendo o componente ",
      "cont\u00ednuo como categoria de refer\u00eancia."
    )
  } else if (has_p0) {
    augmentation_en <- paste0(
      " The probability of an exact zero was modeled through a logistic ",
      "boundary-mass submodel."
    )
    augmentation_pt <- paste0(
      " A probabilidade de zero exato foi modelada por um submodelo log\u00edstico ",
      "para a massa no limite."
    )
  } else if (has_p1) {
    augmentation_en <- paste0(
      " The probability of an exact one was modeled through a logistic ",
      "boundary-mass submodel."
    )
    augmentation_pt <- paste0(
      " A probabilidade de um exato foi modelada por um submodelo log\u00edstico ",
      "para a massa no limite."
    )
  } else {
    augmentation_en <- ""
    augmentation_pt <- ""
  }

  re_components <- names(re_flags)[re_flags]
  spatial_components <- names(spatial_flags)[spatial_flags]

  re_text_en <- if (length(re_components)) {
    paste0(
      " Ordinary Gaussian random effects were included in the ",
      join_labels(labels_en[re_components], "en"),
      if (length(re_components) == 1L) " predictor." else " predictors."
    )
  } else {
    ""
  }

  re_text_pt <- if (length(re_components)) {
    paste0(
      " Efeitos aleat\u00f3rios gaussianos foram inclu\u00eddos ",
      if (length(re_components) == 1L) "no preditor de " else "nos preditores de ",
      join_labels(labels_pt[re_components], "pt"),
      "."
    )
  } else {
    ""
  }

  spatial_text_en <- if (length(spatial_components)) {
    paste0(
      " Proper CAR spatial effects were included in the ",
      join_labels(labels_en[spatial_components], "en"),
      if (length(spatial_components) == 1L) {
        " predictor to account for areal dependence."
      } else {
        " predictors to account for areal dependence."
      }
    )
  } else {
    ""
  }

  spatial_text_pt <- if (length(spatial_components)) {
    paste0(
      " Efeitos espaciais CAR pr\u00f3prios foram inclu\u00eddos ",
      if (length(spatial_components) == 1L) {
        "no preditor de "
      } else {
        "nos preditores de "
      },
      join_labels(labels_pt[spatial_components], "pt"),
      " para acomodar depend\u00eancia espacial areal."
    )
  } else {
    ""
  }

  estimation_en <- if (has_latent) {
    paste0(
      "Parameters were estimated by maximum marginal likelihood using the ",
      "Laplace approximation implemented in Template Model Builder (TMB) with ",
      "automatic differentiation."
    )
  } else {
    paste0(
      "Parameters were estimated by maximum likelihood using Template Model ",
      "Builder (TMB) with automatic differentiation."
    )
  }

  estimation_pt <- if (has_latent) {
    paste0(
      "Os par\u00e2metros foram estimados por m\u00e1xima verossimilhan\u00e7a marginal ",
      "utilizando a aproxima\u00e7\u00e3o de Laplace implementada no Template Model ",
      "Builder (TMB) com diferencia\u00e7\u00e3o autom\u00e1tica."
    )
  } else {
    paste0(
      "Os par\u00e2metros foram estimados por m\u00e1xima verossimilhan\u00e7a utilizando o ",
      "Template Model Builder (TMB) com diferencia\u00e7\u00e3o autom\u00e1tica."
    )
  }

  if (lang == "en") {
    text <- paste0(
      "Data were analyzed using the 'unitregTMB' framework for bounded ",
      "responses. The continuous component was assumed to follow a ",
      family_name,
      " distribution. The ",
      parameter_text_en,
      " was related to covariates through a ",
      mu_link,
      " link. The second distributional parameter (phi) used a ",
      phi_link,
      " link.",
      augmentation_en,
      re_text_en,
      spatial_text_en,
      " ",
      estimation_en
    )
  } else {
    text <- paste0(
      "Os dados foram analisados utilizando o framework 'unitregTMB' para ",
      "respostas limitadas. O componente cont\u00ednuo foi modelado assumindo uma ",
      "distribui\u00e7\u00e3o ",
      family_name,
      ". A ",
      parameter_text_pt,
      " foi relacionada \u00e0s covari\u00e1veis por uma fun\u00e7\u00e3o de liga\u00e7\u00e3o ",
      mu_link,
      ". O segundo par\u00e2metro distribucional (phi) utilizou uma fun\u00e7\u00e3o de ",
      "liga\u00e7\u00e3o ",
      phi_link,
      ".",
      augmentation_pt,
      re_text_pt,
      spatial_text_pt,
      " ",
      estimation_pt
    )
  }

  version <- tryCatch(
    as.character(utils::packageVersion("unitregTMB")),
    error = function(e) "development"
  )

  year <- format(Sys.Date(), "%Y")

  citation_bib <- paste0(
    "@Manual{unitregTMB", year, ",\n",
    "  title = {unitregTMB: Regression Models for Correlated Continuous Bounded Data},\n",
    "  author = {Ricardo Rasmussen Petterle},\n",
    "  year = {", year, "},\n",
    "  note = {R package version ", version, "},\n",
    "  url = {https://github.com/rpetterle/unitregTMB}\n",
    "}"
  )

  wrapped_text <- paste(strwrap(text, width = 80), collapse = "\n")

  cat("\n------------------------------------------------------------\n",
      "SUGGESTED METHODS TEXT:\n",
      "------------------------------------------------------------\n",
      sep = "")
  cat(wrapped_text, "\n")

  cat("\n------------------------------------------------------------\n",
      "CITATION:\n",
      "------------------------------------------------------------\n",
      sep = "")
  cat(citation_bib, "\n")

  invisible(list(text = text, bibtex = citation_bib,
                 parameterization = model_for,
                 random_effect_components = re_components,
                 spatial_components = spatial_components))
}

