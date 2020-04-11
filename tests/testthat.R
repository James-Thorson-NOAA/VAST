
# devtools::install_github("james-thorson/FishStatsUtils", ref="development")
# devtools::install_github("james-thorson/VAST", ref="development")
# devtools::install_local("C:/Users/James.Thorson/Desktop/Git/VAST", force=TRUE, dep=FALSE)
# devtools::install_local("C:/Users/James.Thorson/Desktop/Git/FishStatsUtils", force=TRUE, dep=FALSE)

library(testthat)
library(VAST)

# Version
Version_VAST = get_latest_version( package="VAST" )
TmbDir = system.file("executables", package="VAST")

# Use "extdata" in "inst" because its loaded with R packages
singlespecies_example_path <- system.file("extdata", package="SpatialDeltaGLMM")
multispecies_example_path <- system.file("extdata", package="VAST")

# Check that single-species results are available locally
test_path = file.path(singlespecies_example_path,"EBS_pollock")
file.exists( file.path(test_path,"opt.RData") )

# Run tests for VAST
#setwd(system.file("tests", package="VAST"))
testthat::test_check("VAST")

################
# Local testing
################

# Use local path
#singlespecies_example_path <- "C:/Users/James.Thorson/Desktop/Git/geostatistical_delta-GLMM/inst/extdata/"
#multispecies_example_path <- "C:/Users/James.Thorson/Desktop/Git/VAST/inst/extdata/"

# Run from local directory
#testthat::test_dir( "C:/Users/James.Thorson/Desktop/Git/VAST/tests/testthat/", reporter="check" )
#testthat::test_dir( "/media/sf_c/Users/jim/Desktop/Project_git/VAST/tests/testthat/", reporter="check" )
