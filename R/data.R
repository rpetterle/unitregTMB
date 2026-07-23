#' Body Fat Percentage Data (Wide Format)
#'
#' @description 
#' A dataset containing body fat proportion measurements for 298 individuals, 
#' alongside demographic and biometric covariates. This dataset was originally 
#' analyzed in the Multivariate Quasi-Beta regression paper and features 5 distinct 
#' anatomical regions evaluated for each individual.
#' 
#' @format A \code{data.frame} with 298 observations and 10 variables:
#' \describe{
#'   \item{id}{Individual identifier.}
#'   \item{age}{Age of the individual in years.}
#'   \item{bmi}{Body Mass Index.}
#'   \item{gender}{Factor indicating biological sex (\code{Female}, \code{Male}).}
#'   \item{ipaq}{Factor indicating physical activity level (\code{Sedentary}, \code{Insufficiently active}, \code{Active}).}
#'   \item{arms}{Proportion of body fat in the arms (0 to 1).}
#'   \item{legs}{Proportion of body fat in the legs (0 to 1).}
#'   \item{trunk}{Proportion of trunk body fat (0 to 1).}
#'   \item{android}{Proportion of android body fat (0 to 1).}
#'   \item{gynoid}{Proportion of gynoid body fat (0 to 1).}
#' }
#' 
#' @references
#' Petterle, R. R., Bonat, W. H., Scarpin, C. T., Jonasson, T., & Borba, V. Z. C. (2020). 
#' Multivariate quasi-beta regression models for continuous bounded data. 
#' \emph{International Journal of Biostatistics}, 17(1), 39-53.
#' 
#' Petterle, R. R., Laureano, H. A., da Silva, G. P., & Bonat, W. H. (2021). 
#' Multivariate generalized linear mixed models for continuous bounded outcomes: 
#' Analyzing the body fat percentage data. 
#' \emph{Statistical Methods in Medical Research}, 30(12), 2619-2633.
#' 
#' Petterle, R. R., Taconeli, C. A., da Silva, J. L. P., da Silva, G. P., 
#' Laureano, H. A., & Bonat, W. H. (2021). 
#' Unit gamma mixed regression models for continuous bounded data. 
#' \emph{Journal of Statistical Computation and Simulation}, 92(1), 1-19.
#' 
#' @source \url{https://estatistica.c3sl.ufpr.br/doku.php/publications:papercompanions:multquasibeta}
"bodyfat"

#' Body Fat Percentage Data (Long Format)
#'
#' @description 
#' A longitudinal (stacked) version of the \code{bodyfat} dataset containing 
#' 1490 observations (5 regions x 298 individuals). This format is specifically 
#' designed to be used with mixed-effects models in \code{unitregTMB}, allowing 
#' the use of \code{(1 | id)} to account for intra-subject correlation across 
#' different anatomical regions.
#'
#' @format A \code{data.frame} with 1490 observations and 7 variables:
#' \describe{
#'   \item{id}{Individual identifier (used for random effects).}
#'   \item{age}{Age of the individual in years.}
#'   \item{gender}{Factor indicating biological sex (\code{Female}, \code{Male}).}
#'   \item{ipaq}{Factor indicating physical activity level.}
#'   \item{bmi}{Body Mass Index.}
#'   \item{regions}{Factor indicating the anatomical region of the measurement (\code{Arms}, \code{Legs}, \code{Trunk}, \code{Android}, \code{Gynoid}).}
#'   \item{y}{Proportion of body fat in the specified region (0 to 1).}
#' }
#' 
#' @references
#' Petterle, R. R., Laureano, H. A., da Silva, G. P., & Bonat, W. H. (2021). 
#' Multivariate generalized linear mixed models for continuous bounded outcomes: 
#' Analyzing the body fat percentage data. 
#' \emph{Statistical Methods in Medical Research}, 30(12), 2619-2633.
#' 
#' @source \url{https://estatistica.c3sl.ufpr.br/doku.php/publications:papercompanions:multquasibeta}
#'
#' @examples
#' data("bodyfat_long", package = "unitregTMB")
#' head(bodyfat_long)
#' \donttest{
#' # Example: Fitting a mixed-effects Vasicek model
#' # fit <- unitregTMB(y ~ regions + gender + (1 | id), 
#' #                   data = bodyfat_long, 
#' #                   family = vasicek())
#' }
"bodyfat_long"