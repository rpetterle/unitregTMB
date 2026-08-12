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

#' Periodontal Disease Data Set
#'
#' @description
#' The data come from a study conducted at the Medical University of South Carolina by Fernandes et al. (2006).
#' The data set involves 1160 clustered observations from 290 subjects.
#'
#' @format A data frame with 1160 rows and 8 variables:
#' \describe{
#'   \item{id}{Subject identifier.}
#'   \item{tooth}{Tooth types (molar; premolar; canine; incisor).}
#'   \item{y}{The response variable lies on the closed interval [0,1] and represents the proportion of diseased tooth sites for each of the four tooth types. Contains exact zeros and ones.}
#'   \item{gender}{Patient gender (Male; Female).}
#'   \item{age}{Patient age in years.}
#'   \item{hba1c}{Glycosylated hemoglobin status indicator (Controlled, <7\%; Uncontrolled, >=7\%).}
#'   \item{smoker}{Smoking status (Non-smoker; Smoker).}
#'   \item{bmi}{Body mass index.}
#' }
#'
#' @source Downloaded from Dr. Bandyopadhyay's software page: \url{https://sites.google.com/view/dbandyop/software}.
#' Original dataset archive: \url{https://drive.google.com/file/d/1-fyruHi5eMp6TvFegQazx5lios6igLGc/view?pli=1}.
#'
#' @references
#' Zhao, W., Lian, H., & Bandyopadhyay, D. (2018). A partially linear additive model for clustered proportion data. \emph{Statistics in Medicine}, 37(6), 1009-1030.
#'
#' Fernandes, J., et al. (2006). Prevalence of periodontal disease in Gullah African American diabetics. \emph{Journal of Dental Research}, 85(Special Issue A), 997.
#'
#' @examples
#' data("periodontal", package = "unitregTMB")
#' head(periodontal)
#' \donttest{
#' # Example: Fitting a zero-and-one augmented mixed-effects model
#' # fit <- unitregTMB(y ~ age + gender + tooth + (1 | id), 
#' #                   data = periodontal, 
#' #                   family = unitgamma())
#' }
"periodontal"

