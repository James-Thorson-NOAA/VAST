library(testthat)
library(VAST)
library(TMB)

# Version
Version_VAST = "VAST_v1_8_0"

# Use "extdata" in "inst" because its loaded with R packages
example_path <- system.file("extdata", package="SpatialDeltaGLMM")

# setwd(system.file("tests", package="SpatialDeltaGLMM"))
testthat::test_check("VAST")

# Local testing
if(FALSE){
  # Use local path
  example_path <- "C:/Users/James.Thorson/Desktop/Project_git/geostatistical_delta-GLMM/inst/extdata/"
  # Run from local directory
  testthat::test_dir( "C:/Users/James.Thorson/Desktop/Project_git/VAST/tests/testthat/", reporter="check" )
}
