#' Vector Autoregressive Spatio-Temporal Model
#'
#' Analyze spatial variatoin in multiple variables representing data species, ages, size-classes, etc.
#' This can be used to implement many types of analysis including:
#' \itemize{
#' \item species distribution models with residual spatial and spatio-temporal variation, as well as covariates affecting density and catchability
#' \item joint species distribution models (both static and dynamic)
#' \item empirical orthogonal function (EOF) analysis,
#' \item vector autoregressive modelling,
#' \item spatial factor analysis,
#' \item standardization and biomass-expansion of spatially unbalanced samples
#' }
#' Features are built to be compatible among model types, e.g., by allowing catchability and density covariates to be included in EOF analysis.
#'
#' See \code{\link{fit_model}} for a simple example of high-level wrapper functions for using VAST.
#'
#' @seealso \code{\link[VAST]{VAST}} for general documentation, \code{\link{make_settings}} for generic settings, \code{\link{fit_model}} for model fitting, and \code{\link{plot_results}} for generic plots
#' @seealso VAST wiki \url{https://github.com/James-Thorson-NOAA/VAST/wiki} for examples documenting many different use-cases and features.
#' @seealso GitHub mainpage \url{https://github.com/James-Thorson-NOAA/VAST#description} for a list of user resources and publications documenting features
#' @name VAST
#' @importFrom utils packageVersion packageDescription
#' @importFrom stats gaussian model.matrix poisson update.formula
#' @importFrom grDevices jpeg rainbow rgb tiff
#' @importFrom graphics abline arrows contour legend lines matplot plot.new points polygon rect
#' @importFrom methods as formalArgs is slot
#' @importFrom stats .getXlevels cor cutree delete.response dist hclust lm median model.frame optimHess pchisq pgamma plnorm plogis pnorm qnorm quantile rbinom rgamma rlnorm rpois runif sd terms varimax
#' @importFrom utils data example head tail write.csv
NULL
