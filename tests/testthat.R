library(testthat)
library(VAST)
library(TMB)

# Version
Version_VAST = "VAST_v1_8_0"

# Use "extdata" in "inst" because its loaded with R packages
singlespecies_example_path <- system.file("extdata", package="SpatialDeltaGLMM")
multispecies_example_path <- system.file("extdata", package="VAST")

# Run tests for VAST
#setwd(system.file("tests", package="VAST"))
testthat::test_check("VAST")

# Local testing
  # Use local path
  #example_path <- "C:/Users/James.Thorson/Desktop/Project_git/VAST/inst/extdata/"

  # Run from local directory
  #testthat::test_dir( "C:/Users/James.Thorson/Desktop/Project_git/VAST/tests/testthat/", reporter="check" )
