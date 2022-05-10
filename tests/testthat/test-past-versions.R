
############################
# Modifications from SpatialDeltaGLMM version
# 1. Change Version to Version_VAST
# 2. Change Data_Fn and Build_TMB_Fn to call VAST instead of SpatialDeltaGLMM
# 3. Change VesselConfig to OverdispersionConfig, and move to Data_Fn
# 4. Add c_i to Data_Fn
# 5. Change ObsModel to c(ObsModel,0)
# 6. Change expect_equal
# 7. Change v_i from Vessel to VesselxYear, because SpatialDeltaGLMM did this automatically for VesselConfig=c(0,1)
# 8. Change tolerance on Chatham Rise example (because SpatialDeltaGLMM hit bound and isn't *quite* converged
############################

# Tutorial: http://r-pkgs.had.co.nz/tests.html
# And example see: https://github.com/ss3sim/ss3sim/tree/master/tests/testthat
context("Testing examples")

# Eastern Bering Sea pollcok
test_that("Eastern Bering Sea pollock is working ", {
  # Previously worked with CI, but not anymore
  #skip_on_ci()
  skip_if(skip_local)

  # Prepping
  test_path = file.path(singlespecies_example_path,"EBS_pollock")

  # load data set
  example = load_example( data_set="EBS_pollock" )

  # Make settings
  settings = make_settings( n_x=50,
                            Region=example$Region,
                            purpose="index2" )
  settings$FieldConfig[c("Omega","Epsilon"),"Component_2"] = 0

  # Add
  version_set = list.files(system.file("executables", package = "VAST"))
    version_set = sapply( version_set, FUN=function(char){strsplit(char,".",fixed=TRUE)[[1]][1]} )
  #version_set = c(
  #  paste0("VAST_v8_",0:6,"_0"),
  #  paste0("VAST_v9_",0:4,"_0"),
  #  "VAST_v10_0_0",
  #  "VAST_v11_0_0",
  #  "VAST_v12_0_0",
  #  paste0("VAST_v13_",0:1,"_0"),
  #  "VAST_v14_0_0",
  #  "VAST_v14_0_1" )
  version_set = setdiff( version_set, get_latest_version() )

  cpp_exists = paste0(version_set,".cpp") %in% list.files(system.file("executables", package = "VAST"))
  expect_equal( all(cpp_exists), TRUE )

  for( vI in seq_along(version_set) ){
    # Run model
    settings$Version = version_set[vI]
    fit = fit_model( "settings"=settings,
                     "Lat_i"=example$sampling_data[,'Lat'],
                     "Lon_i"=example$sampling_data[,'Lon'],
                     "t_i"=example$sampling_data[,'Year'],
                     "b_i"=example$sampling_data[,'Catch_KG'],
                     "a_i"=example$sampling_data[,'AreaSwept_km2'],
                     "newtonsteps" = 0,
                     "test_fit" = FALSE,
                     "getsd" = FALSE,
                     "working_dir" = test_path )
    expect_equal( fit$parameter_estimates$objective, 58883.63, tol=1 )
  }
})

