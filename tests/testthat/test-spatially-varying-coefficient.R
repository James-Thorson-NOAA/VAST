
context("Testing examples")

# Eastern Bering Sea pollcok
test_that("Spatially varying coefficient example is working ", {
  skip_on_travis()

  # Prepping
  test_path = file.path(multispecies_example_path,"Spatially_varying_coefficient")
  load( file=file.path(test_path,"Data.RData") )
  load( file.path(test_path,"saved_estimates.RData") )
  load( file.path(test_path,"settings.RData") )
  settings$Version = FishStatsUtils::get_latest_version()
  attach(settings)
  on.exit( detach(settings) )

  # Renumber years to have same indexing as covariates
  t_i = match( Data$sampling_data[,'Year'], sort(unique(Data$sampling_data[,'Year'])) )

  # Run model
  fit = fit_model( "settings"=settings, "Lat_i"=Data$sampling_data[,'Lat'],
    "Lon_i"=Data$sampling_data[,'Lon'], "t_i"=t_i,
    "c_i"=rep(0,nrow(Data$sampling_data)), "b_i"=Data$sampling_data[,'Catch_KG'],
    "a_i"=Data$sampling_data[,'AreaSwept_km2'], "v_i"=Data$sampling_data[,'Vessel'],
    "X_gtp"=Data$X_gtp, "X_itp"=Data$X_itp, "Xconfig_zcp"=Data$Xconfig_zcp, optimize_args=list("getsd"=FALSE),
    "working_dir"=test_path, "test_fit"=FALSE )

  # Comparisons
  Par1 = fit$parameter_estimates$par[names(fit$parameter_estimates$par)%in%c("ln_H_input","beta1_ft","logkappa1","gamma2_ctp","beta2_ft","logkappa2","log_sigmaXi2_cp","logSigmaM")]
  Par2 = parameter_estimates$par[names(parameter_estimates$par)%in%c("ln_H_input","beta1_ft","logkappa1","gamma2_cp","beta2_ft","logkappa2","log_sigmaXi2_cp","logSigmaM")]
  expect_equal( as.vector(Par1), as.vector(Par2), tolerance=1e-3 )
  expect_equal( parameter_estimates$objective, fit$parameter_estimates$objective, tolerance=1e-3 )

  # Get expanded composition estimates and input sample sizes
  #Opt = TMBhelper::fit_tmb( obj=fit$tmb_list$Obj, startpar=fit$parameter_estimates$par, getsd=TRUE, lower=fit$tmb_list$Lower, upper=fit$tmb_list$Upper, newtonsteps=0 )
  #Index = plot_biomass_index( TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, DirName=test_path )
  #calculate_proportion( TmbData=fit$data_list, Index=Index, DirName=test_path, Year_Set=fit$year_labels )
})

