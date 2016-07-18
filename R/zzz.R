
#.onLoad <- function(libname, pkgname) {
#}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Loading package VAST, developed by James Thorson for the Northwest Fisheries Science Center")
  packageStartupMessage("For details and citation guidance, please see http://github.com/james-thorson/VAST/")
  if( !"INLA" %in% installed.packages()[,1] ){
    message("Installing INLA...")
    install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")  
  }
  if( !"TMB" %in% installed.packages()[,1] ){
    message("Installing TMB...")
    devtools::install_github("kaskr/adcomp/TMB")
  }
}
