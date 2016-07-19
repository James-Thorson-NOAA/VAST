
#.onLoad <- function(libname, pkgname) {
#}

.onAttach <- function(libname, pkgname) {
  print("###########################################################################################")
  print("Loading package VAST, developed by James Thorson for the Northwest Fisheries Science Center")
  print("For details and citation guidance, please see http://github.com/james-thorson/VAST/")
  print("###########################################################################################")
  if( !"INLA" %in% installed.packages()[,1] ){
    print("Installing package: INLA...")
    install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")  
  }
  if( !"TMB" %in% installed.packages()[,1] ){
    print("Installing package: TMB...")
    devtools::install_github("kaskr/adcomp/TMB")
  }
  if( !"SpatialDeltaGLMM" %in% installed.packages()[,1] ){
    print("Installing package: SpatialDeltaGLMM...")
    devtools::install_github("nwfsc-assess/geostatistical_delta-GLMM")
  }
}
