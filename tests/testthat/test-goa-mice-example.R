
context("Testing examples")

# Eastern Bering Sea pollcok
test_that("Gulf of Alaska MICE-in-space example is working ", {
  #skip_on_ci()
  skip_if(skip_local)

  # Prepping
  test_path = file.path(multispecies_example_path,"goa_mice_example")
  example = load_example( "GOA_MICE_example" )
  load( file.path(test_path,"saved_estimates.RData") )
  load( file.path(test_path,"settings.RData") )
  #settings$Version = FishStatsUtils::get_latest_version()
  settings$Version = Version_VAST
  attach(settings)
  on.exit( detach(settings) )

  # Run model
  Fit = fit_model( "settings"=settings, "Lat_i"=example$sampling_data[,'Lat'],
    "Lon_i"=example$sampling_data[,'Lon'], "b_i"=example$sampling_data[,'Catch_KG'],
    "a_i"=example$sampling_data[,'AreaSwept_km2'], "v_i"=as.numeric(example$sampling_data[,'Vessel']),
    "t_i"=example$sampling_data[,'Year'], "c_i"=as.numeric(example$sampling_data[,'spp'])-1, "F_ct"=example$F_ct,
    "newtonsteps"=0, optimize_args=list("getsd"=FALSE),
    "working_dir"=test_path, "run_model"=FALSE )

  # Comparisons
  #Par1 = fit$parameter_estimates$par[names(fit$parameter_estimates$par)%in%c("ln_H_input","beta1_ft","logkappa1","gamma2_ctp","beta2_ft","logkappa2","log_sigmaXi2_cp","logSigmaM")]
  #Par2 = parameter_estimates$par[names(parameter_estimates$par)%in%c("ln_H_input","beta1_ft","logkappa1","gamma2_ctp","beta2_ft","logkappa2","log_sigmaXi2_cp","logSigmaM")]
  #expect_equal( as.vector(Par1), as.vector(Par2), tolerance=1e-3 )
  Obj = Fit$tmb_list$Obj$fn( parameter_estimates$par )
  expect_equal( parameter_estimates$objective, as.numeric(Obj), tolerance=1e-3 )
})

