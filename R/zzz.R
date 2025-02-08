
# .onLoad <- function(libname, pkgname) {
# }

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("###########################################################################################")
  packageStartupMessage("Loading package VAST version ", packageVersion("VAST") )
  packageStartupMessage("For information and examples, please see http://github.com/james-thorson/VAST/")
  packageStartupMessage("###########################################################################################")
  #if( getOption("repos")["CRAN"] == "@CRAN@" ){
  #  options(repos = c("CRAN" = "http://cran.us.r-project.org"))
  #}
  #
  # URL for ZIP files to download and install locally: https://inla.r-inla-download.org/R/stable/src/contrib/
  #if( !"INLA" %in% utils::installed.packages()[,1] ){
  #  packageStartupMessage("Installing package: INLA...")
  #
  #  # Over-ride default install for R 3.5.0 through R 3.5.3
  #  Rvers = numeric_version(paste0(R.version[6:7],collapse="."))
  #  if( Rvers<numeric_version("3.6.0") & Rvers>numeric_version("3.5.0") ){
  #    utils::install.packages( "https://inla.r-inla-download.org/R/stable/bin/windows/contrib/3.5/INLA_18.07.12.zip" )
  #  }else{
  #    utils::install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
  #  }
  #}
  #
  ## Load `FishStatsUtils` via .onAttach because Remotes doessn't enforce branch properly
  #if( !"FishStatsUtils" %in% utils::installed.packages()[,1] || utils::packageVersion("FishStatsUtils") < numeric_version("2.11.0") ){
  #  packageStartupMessage("Updating package FishStatsUtils because previously using version < 2.11.0")
  #  remotes::install_github("james-thorson-NOAA/FishStatsUtils", ref="dev")
  #}
  #packageStartupMessage( "Loading package `FishStatsUtils` version ", packageVersion("FishStatsUtils") )
  #library(FishStatsUtils)
}

#' Copy of VAST::make_model
#'
#' Included for continuity when using old scripts
#'
#' Please use \code{?VAST::make_model} to see list of arguments
#' and usage
#' @param ... Arguments to be passed to \code{make_model}
#' @export
Build_TMB_Fn = function( ... ){
  .Deprecated( new="VAST::make_model" )
  VAST::make_model( ... )
}

#' Copy of VAST::make_data
#'
#' Included for continuity when using old scripts
#'
#' Please use \code{?VAST::make_data} to see list of arguments
#' and usage
#' @param ... Arguments to be passed to \code{make_data}
#' @export
Data_Fn = function( ... ){
  .Deprecated( new="VAST::make_data" )
  VAST::make_data( ... )
}

#' Copy of FishStatsUtils::plot_factors
#'
#' Included for continuity when using old scripts
#'
#' Please use \code{?FishStatsUtils::plot_factors} to see list of arguments and usage
#' @param ... Arguments to be passed to \code{plot_factors}
#' @export
Plot_factors = function( ... ){
  .Deprecated( new="FishStatsUtils::plot_factors" )
  FishStatsUtils::plot_factors( ... )
}

#' Copy of FishStatsUtils::plot_overdispersion
#'
#' Included for continuity when using old scripts
#'
#' Please use \code{?FishStatsUtils::plot_overdispersion} to see list of arguments and usage
#' @param ... Arguments to be passed to \code{plot_overdispersion}
#' @export
Plot_Overdispersion = function( ... ){
  .Deprecated( new="FishStatsUtils::plot_overdispersion" )
  FishStatsUtils::plot_overdispersion( ... )
}

#' Copy of FishStatsUtils::summarize_covariance
#'
#' Included for continuity when using old scripts
#'
#' Please use \code{?FishStatsUtils::summarize_covariace} to see list of arguments and usage
#' @param ... Arguments to be passed to \code{summarize_covariance}
#' @export
Summarize_Covariance = function( ... ){
  .Deprecated( new="FishStatsUtils::summarize_covariance" )
  FishStatsUtils::summarize_covariance( ... )
}

#' Copy of make_map
#'
#' Included for continuity when using old scripts
#'
#' Please use \code{?make_map} to see list of arguments and usage
#' @param ... Arguments to be passed to \code{make_map}
#' @noRd
## Note: noRD tag needed otherwise this overwrites documentation for
##  make_map b/c for some reason it's not case sensitive in the
##   documentation -Cole
Make_Map = function( ... ){
  .Deprecated( new="make_map" )
  make_map( ... )
}

#' Copy of make_parameters
#'
#' Included for continuity when using old scripts
#'
#' Please use \code{?make_parameters} to see list of arguments and usage
#' @param ... Arguments to be passed to \code{make_parameters}
#' @export
Param_Fn = function( ... ){
  .Deprecated( new="make_parameters" )
  make_parameters( ... )
}

PlotLoadings = function( ... ){
  .Deprecated( new="FishStatsUtils::plot_loadings" )
  FishStatsUtils::plot_loadings( ... )
}
Rotate_Fn = function( ... ){
  .Deprecated( new="FishStatsUtils::rotate_factors" )
  FishStatsUtils::rotate_factors( ... )
}
Calc_Anisotropic_Mesh = function( ... ){
  .Deprecated( new="FishStatsUtils::make_mesh" )
  FishStatsUtils::make_mesh( ... )
}
Calc_Kmeans = function( ... ){
  .Deprecated( new="FishStatsUtils::make_kmeans" )
  FishStatsUtils::make_kmeans( ... )
}


