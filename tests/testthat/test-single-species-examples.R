
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
  if(!require("INLA")){
    install.packages("INLA", dep=TRUE, repos=c( CRAN="https://cloud.r-project.org",
                                                INLA="https://inla.r-inla-download.org/R/stable") )
  }

  # Prepping
  test_path = file.path(singlespecies_example_path,"EBS_pollock")
  load( file.path(test_path,"opt.RData") )
  load( file.path(test_path,"Record.RData") )
  attach(Record)
  on.exit( detach(Record) )
  # Run model
  data( EBS_pollock_data, package="VAST" )
  EBS_pollock_data = EBS_pollock_data$sampling_data
  Data_Geostat = data.frame( "Catch_KG"=EBS_pollock_data[,'catch'], "Year"=EBS_pollock_data[,'year'], "Vessel"="missing", "AreaSwept_km2"=0.01, "Lat"=EBS_pollock_data[,'lat'], "Lon"=EBS_pollock_data[,'long'], "Pass"=0)
  Extrapolation_List = make_extrapolation_info( Region=Region, strata.limits=strata.limits )
  Spatial_List = make_spatial_info( backwards_compatible_kmeans=TRUE, grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=test_path )
  TmbData = make_data("Version"=Version_VAST, "OverdispersionConfig"=rep(VesselConfig[2],2), "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "ObsModel"=c(ObsModel,0), "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=as.numeric(factor(paste(Data_Geostat[,'Vessel'],Data_Geostat[,'Year'])))-1, "t_i"=Data_Geostat[,'Year'], "spatial_list"=Spatial_List )
  TmbList = make_model("TmbData"=TmbData, "build_model"=TRUE, "RunDir"=test_path, "Version"=Version_VAST, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x ) #, "Parameters"=Params)
  #on.exit( dyn.unload(paste0(system.file("executables", package = "VAST"),"/",TMB::dynlib(Version_VAST))), add=TRUE )
  Opt = fit_tmb( obj=TmbList[["Obj"]], getsd=FALSE, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]] )  # , rel.tol=1e-20
  # Comparisons
  Par1 = Opt$par[names(Opt$par)%in%c("ln_H_input","beta1_ct","beta1_ft","logkappa1","beta2_ct","beta2_ft","logkappa1","logSigmaM")]
  Par2 = opt$par[names(opt$par)%in%c("ln_H_input","beta1_t","logkappa1","beta2_t","logkappa1","logSigmaM")]
  expect_equal( as.vector(Par1), as.vector(Par2), tolerance=1e-3 )
})

