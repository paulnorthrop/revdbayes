#' revdbayes: Ratio-of-Uniforms Sampling for Bayesian Extreme Value Analysis
#'
#' Uses the multivariate generalized ratio-of-uniforms method to simulate
#' random samples from the posterior distributions commonly encountered in
#' Bayesian extreme value analyses.
#'
#' @details The main function in the revbayes package is \code{rpost}, which
#'   simulates random samples from the posterior distribution of extreme value
#'   model parameters using the function \code{ru} from the rust package. The
#'   user chooses the extreme value model, the prior density for the parameters
#'   and provides the data.  There are options to improve the probability of
#'   acceptance of the ratio-of-uniforms algorithm by working with
#'   transformation of the model parameters.
#'
#'   See \code{vignette("revdbayes-vignette", package = "revdbayes")} for an
#'   overview of the package.
#'
#' @references Northrop, P. J. (2016). rust: Ratio-of-Uniforms Simulation with
#'   Transformation. R package version 1.1.0.
#'   \url{https://cran.r-project.org/package=rust}.
#'
#' @seealso \code{\link{set_prior}} to set a prior density for extreme value
#'   parameters.
#' @seealso \code{\link{rpost}} to perform ratio-of-uniforms sampling from
#'   an extreme value posterior distribution.
#' @seealso The \code{\link[rust]{ru}} function in the \code{\link{rust}}
#'   package for details of the arguments that can be passed to \code{ru}
#'   via \code{rpost} and for the form of the object (of class "evprior")
#'   returned from \code{rpost}, which has the same structure as an object
#'   (of class "ru") returned by \code{ru}.
#' @docType package
#' @name revdbayes
#' @import methods
#' @importFrom stats runif
NULL

#' Annual Maximum Sea Levels at Port Pirie, South Australia
#'
#' A numeric vector of length 65 containing annual maximum sea levels,
#' in metres, from 1923 to 1987 at Port Pirie, South Australia.
#'
#' @format A numeric vector containing 65 observations.
#' @source Coles, S. G. (2001) \emph{An Introduction to Statistical Modelling
#'   of Extreme Values}. London: Springer.
#'   doi:\href{https://doi.org/10.1007/978-1-4471-3675-0}{10.1007/978-1-4471-3675-0}.
"portpirie"

#' Annual Maximum Temperatures at Oxford
#'
#' A numeric vector containing annual maximum temperatures, in degrees
#' Fahrenheit, from 1901 to 1980 at Oxford, England.
#'
#'@format A vector containing 80 observations.
#'@source Tabony, R. C. (1983) Extreme value analysis in meteorology.
#'  \emph{The Meteorological Magazine}, \strong{112}, 77-98.
"oxford"

#' Daily Aggregate Rainfall
#'
#' A numeric vector of length 20820 containing daily aggregate rainfall
#' observations, in millimetres, recorded at a rain gauge in England
#' over a 57 year period, beginning on a leap year. Three of these years
#' contain only missing values.
#'
#' @format A vector containing 20820 observations.
#' @source Unknown
"rainfall"

#' Storm peak significant wave heights from the Gulf of Mexico
#'
#' A numeric vector containing 315 hindcasts of storm peak significant wave
#' heights, metres, from 1900 to 2005 at an unnamed location in the Gulf
#' of Mexico.
#'
#'@format A vector containing 315 observations.
#'@source Oceanweather Inc. (2005) GOMOS -- Gulf of Mexico hindcast study.
#'@references Northrop, P. J., N. Attalides, and P. Jonathan. (2016).
#'  Cross-Validatory Extreme Value Threshold Selection and Uncertainty with
#'  Application to Ocean Storm Severity. \emph{Journal of the Royal
#'  Statistical Society: Series C (Applied Statistics)}.
#'  doi:\href{https://doi.org/10.1111/rssc.12159}{10.1111/rssc.12159}.
"gom"

#' Largest Sea Levels in Venice
#'
#' The \code{venice} data frame has 51 rows and 10 columns. The jth column
#' contains the jth largest sea levels in Venice, for the years 1931-1981. Only
#' the largest six measurements are available for the year 1935; the
#' corresponding row contains four missing values. The years for each set of
#' measurements are given as row names.
#'
#' @format A data frame with 51 rows and 10 columns.
#' @source Smith, R. L. (1986) Extreme value theory based on the \emph{r}
#'   largest annual events. \emph{Journal of Hydrology}, \strong{86}, 27-43.
#'   doi:\href{https://doi.org/10.1016/0022-1694(86)90004-1}{10.1016/0022-1694(86)90004-1}.
#'
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modelling of Extreme Values}. London: Springer.
#'   doi:\href{https://doi.org/10.1007/978-1-4471-3675-0}{10.1007/978-1-4471-3675-0}.
"venice"
