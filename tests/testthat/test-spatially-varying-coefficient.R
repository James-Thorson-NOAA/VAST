
context("Testing examples")

# Eastern Bering Sea pollcok
test_that("Spatially varying coefficient example is working ", {
  skip_on_ci()
  skip_if(skip_local)


  # Prepping
  test_path = file.path(multispecies_example_path,"Spatially_varying_coefficient")
  load( file=file.path(test_path,"Data.RData") )
  load( file.path(test_path,"saved_estimates.RData") )
  load( file.path(test_path,"settings.RData") )
  settings$Version = FishStatsUtils::get_latest_version()
  settings$Options = c( settings$Options, "report_additional_variables"=TRUE )
  attach(settings)
  on.exit( detach(settings) )

  # Renumber years to have same indexing as covariates
  t_i = match( Data$sampling_data[,'Year'], sort(unique(Data$sampling_data[,'Year'])) )
  covariate_data = data.frame( "Lat"=Data$sampling_data[,'Lat'], "Lon"=Data$sampling_data[,'Lon'], "Year"=t_i,
    "var1"=Data$X_itp[,1,1], "var2"=Data$X_itp[,1,2], "var3"=Data$X_itp[,1,3] )

  # Run model -- spatially varying density
  fitX = fit_model( "settings"=settings, "Lat_i"=Data$sampling_data[,'Lat'],
    "Lon_i"=Data$sampling_data[,'Lon'], "t_i"=t_i,
    "c_i"=rep(0,nrow(Data$sampling_data)), "b_i"=Data$sampling_data[,'Catch_KG'],
    "a_i"=Data$sampling_data[,'AreaSwept_km2'], "v_i"=Data$sampling_data[,'Vessel'],
    "X_gtp"=Data$X_gtp, "X_itp"=Data$X_itp, "Xconfig_zcp"=Data$Xconfig_zcp, "getsd"=FALSE,
    "working_dir"=test_path, "test_fit"=FALSE )

  # Run model -- spatially varying catchability
  fitQ = fit_model( "settings"=settings, "Lat_i"=Data$sampling_data[,'Lat'],
    "Lon_i"=Data$sampling_data[,'Lon'], "t_i"=t_i,
    "c_i"=rep(0,nrow(Data$sampling_data)), "b_i"=Data$sampling_data[,'Catch_KG'],
    "a_i"=Data$sampling_data[,'AreaSwept_km2'], "v_i"=Data$sampling_data[,'Vessel'],
    "Q_ik"=array(Data$X_itp[,1,],dim=dim(Data$X_itp)[c(1,3)]),
    "Q1config_k"=Data$Xconfig_zcp[1,1,], "Q2config_k"=Data$Xconfig_zcp[2,1,], "getsd"=FALSE,
    "working_dir"=test_path, "test_fit"=FALSE )

  # Run model -- spatially varying density -- FORMULA interface
  fitX_formula = fit_model( "settings"=settings, "Lat_i"=Data$sampling_data[,'Lat'],
    "Lon_i"=Data$sampling_data[,'Lon'], "t_i"=t_i,
    "c_i"=rep(0,nrow(Data$sampling_data)), "b_i"=Data$sampling_data[,'Catch_KG'],
    "a_i"=Data$sampling_data[,'AreaSwept_km2'], "v_i"=Data$sampling_data[,'Vessel'],
    "covariate_data"=covariate_data,
    "X2_formula"=~var1+var2+var3, "X2config_cp"=matrix(Data$Xconfig_zcp[2,,],nrow=1),
    "getsd"=FALSE, "working_dir"=test_path, "test_fit"=FALSE )

  # Run model -- spatially varying catchability -- FORMULA interface
  fitQ_formula = fit_model( "settings"=settings, "Lat_i"=Data$sampling_data[,'Lat'],
    "Lon_i"=Data$sampling_data[,'Lon'], "t_i"=t_i,
    "c_i"=rep(0,nrow(Data$sampling_data)), "b_i"=Data$sampling_data[,'Catch_KG'],
    "a_i"=Data$sampling_data[,'AreaSwept_km2'], "v_i"=Data$sampling_data[,'Vessel'],
    "catchability_data"=covariate_data,
    "Q2_formula"=~var1+var2+var3, "Q2config_k"=Data$Xconfig_zcp[2,1,],
    "getsd"=FALSE, "working_dir"=test_path, "test_fit"=FALSE )

  # Comparisons
  Par_orig = parameter_estimates$par[names(parameter_estimates$par)%in%c("ln_H_input","beta1_ft","logkappa1","gamma2_cp","beta2_ft","logkappa2","log_sigmaXi2_cp","logSigmaM")]
  ParX = fitX$parameter_estimates$par[names(fitX$parameter_estimates$par)%in%c("ln_H_input","beta1_ft","logkappa1","gamma2_ctp","beta2_ft","logkappa2","log_sigmaXi2_cp","logSigmaM")]
  ParQ = fitQ$parameter_estimates$par[names(fitQ$parameter_estimates$par)%in%c("ln_H_input","beta1_ft","logkappa1","gamma2_ctp","beta2_ft","logkappa2","log_sigmaPhi2_k","logSigmaM")]

  expect_equal( as.vector(Par_orig), as.vector(ParX), tolerance=1e-3 )
  expect_equal( as.vector(Par_orig), as.vector(ParQ), tolerance=1e-3 )
  expect_equal( fitX$parameter_estimates$objective, parameter_estimates$objective, tolerance=1e-3 )
  expect_equal( fitX$parameter_estimates$objective, fitQ$parameter_estimates$objective, tolerance=1e-3 )
  expect_equal( fitX$parameter_estimates$objective, fitX_formula$parameter_estimates$objective, tolerance=1e-3 )
  expect_equal( fitQ$parameter_estimates$objective, fitQ_formula$parameter_estimates$objective, tolerance=1e-3 )

})

