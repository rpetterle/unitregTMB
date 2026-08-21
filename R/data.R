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

#' Water Quality Index (WQI) Data
#'
#' @description 
#' A longitudinal dataset containing the Water Quality Index (WQI) evaluated 
#' at different spatial locations and time periods (quarters). The response 
#' variable is strictly bounded in the open unit interval (0, 1). This dataset 
#' has been extensively used in the literature to illustrate mixed-effects 
#' models for continuous bounded data.
#' 
#' @format A \code{data.frame} containing the variables:
#' \describe{
#'   \item{id}{Factor. The sample/unit identifier.}
#'   \item{y}{Numeric. The calculated Water Quality Index (WQI), bounded between 0 and 1.}
#'   \item{quarter}{Factor. The time of measurement (1, 2, 3, or 4).}
#'   \item{location}{Factor. The spatial measurement site (\code{Upstream}, \code{Reservoir}, \code{Downstream}).}
#' }
#' 
#' @source Downloaded from \url{http://leg.ufpr.br/doku.php/publications:papercompanions:betamix}
#' 
#' @references
#' Bonat, W. H., Ribeiro Jr, P. J., & Zeviani, W. M. (2015). 
#' Likelihood analysis for a class of beta mixed models. 
#' \emph{Journal of Applied Statistics}, 42(2), 252-266. 
#' \doi{10.1080/02664763.2014.947248}
#' 
#' Bonat, W. H., Lopes, J. E., Shimakura, S. E. (2018). 
#' Likelihood analysis for a class of simplex mixed models. 
#' \emph{Chilean Journal of Statistics}, 9(1), 3-17.
#' 
#' Bonat, W. H., et al. (2019). 
#' Flexible quasi-beta regression models for continuous bounded data. 
#' \emph{Statistical Modelling}, 19(5), 525-543. 
#' \doi{10.1177/1471082X18790847}
#' 
#' Petterle, R. R., et al. (2019). 
#' Quasi-beta longitudinal regression model applied to water quality index data. 
#' \emph{Journal of Agricultural, Biological and Environmental Statistics}, 24(3), 469-487. 
#' \doi{10.1007/s13253-019-00360-8}
#' 
#' Petterle, R. R., Taconeli, C. A., da Silva, J. L. P., da Silva, G. P., 
#' Laureano, H. A., & Bonat, W. H. (2021). Unit gamma mixed regression models 
#' for continuous bounded data. \emph{Journal of Statistical Computation and Simulation}, 92(1), 1-19.
#'  
#' @examples
#' data("wqi", package = "unitregTMB")
#' head(wqi)
#' \donttest{
#' # Example: Fitting a mixed-effects Beta model for the mean
#' fit <- unitregTMB(y ~ location + quarter + (1 | id), 
#'                   family = unitgamma(model_for = "mean"),
#'                   data = wqi)
#' }
"wqi"