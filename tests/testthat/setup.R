## Stuff to run before executing the tests
library(FishStatsUtils)
## skip the main tests? caught with skip_if(skip_local) in all
## tests except platform
skip_local <- FALSE
singlespecies_example_path <- system.file("extdata", package="FishStatsUtils")
multispecies_example_path <- system.file("extdata", package="VAST")
## working directory is VAST/tests/testthat when running test()
## Version_VAST <- FishStatsUtils::get_latest_version()
Version_VAST <- FishStatsUtils::get_latest_version(path='../../inst/executables/')
message("Using model version ", Version_VAST)
