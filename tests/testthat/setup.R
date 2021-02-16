## Stuff to run before executing the tests
library(FishStatsUtils)
## skip the main tests? caught with skip_if(skip_local) in all
## tests except platform
skip_local <- TRUE
singlespecies_example_path <- system.file("extdata", package="SpatialDeltaGLMM")
multispecies_example_path <- system.file("extdata", package="VAST")
Version <- FishStatsUtils::get_latest_version(path='../../inst/executables/')
