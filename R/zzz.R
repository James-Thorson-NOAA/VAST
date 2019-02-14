
#.onLoad <- function(libname, pkgname) {
#}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("###########################################################################################")
  packageStartupMessage("Loading package VAST, developed by James Thorson for the Northwest Fisheries Science Center")
  packageStartupMessage("For details and citation guidance, please see http://github.com/james-thorson/VAST/")
  packageStartupMessage("###########################################################################################")
  if( getOption("repos")["CRAN"] == "@CRAN@" ){
    options(repos = c("CRAN" = "http://cran.us.r-project.org"))
  }
  if( !"INLA" %in% utils::installed.packages()[,1] ){
    packageStartupMessage("Installing package: INLA...")
    #utils::install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")
    utils::install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
  }
  #if( !"SpatialDeltaGLMM"%in%utils::installed.packages()[,1] || utils::packageVersion("SpatialDeltaGLMM")<3.40 ){
  #  packageStartupMessage("Installing package: SpatialDeltaGLMM...")
  #  devtools::install_github("nwfsc-assess/geostatistical_delta-GLMM")
  #}
  #if( !"SpatialDFA"%in%utils::installed.packages()[,1] ){
  #  packageStartupMessage("Installing package: SpatialDFA...")
  #  devtools::install_github("james-thorson/spatial_DFA")
  #}
  #if( !"ThorsonUtilities"%in%utils::installed.packages()[,1] ){
  #  packageStartupMessage("Installing package: ThorsonUtilities...")
  #  devtools::install_github("james-thorson/utilities")
  #}
  #if( !"TMBhelper"%in%utils::installed.packages()[,1] ){
  #  packageStartupMessage("Installing package: TMBhelper...")
  #  devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
  #}
}
