
context("Testing examples")

# Eastern Bering Sea pollcok
test_that("Male lingcod compositional expansion is working ", {
  skip_on_ci()
  skip_if(skip_local)

  # Prepping
  test_path = file.path(multispecies_example_path,"Lingcod_comp_expansion")
  load( file=file.path(test_path,"Data_Geostat.RData") )
  load( file.path(test_path,"opt.RData") )
  load( file.path(test_path,"Record.RData") )
  attach(Record)
  on.exit( detach(Record) )
  file.copy( from=paste0(test_path,"/Kmeans-50.RData"), to=paste0(test_path,"/Kmeans_knots-50.RData") )

  # Make settings
  settings = make_settings( n_x = n_x,
           Region = "California_current",
           Version = Version_VAST )

  # Change settings from defaults
  settings$ObsModel = ObsModel_Est
  settings$use_anisotropy = Aniso
  settings$fine_scale = FALSE

  # Run model
  fit = fit_model( "settings" = settings,
      "Lat_i" = Data_Geostat[,'BEST_LAT_DD'],
      "Lon_i" = Data_Geostat[,'BEST_LON_DD'],
      "c_i" = as.numeric(Data_Geostat[,'Length_bin'])-1,
      "b_i" = Data_Geostat[,'First_stage_expanded_numbers'],
      "a_i" = Data_Geostat[,'AreaSwept_km2'],
      "v_i" = rep(0,nrow(Data_Geostat)),
      "t_i" = Data_Geostat[,'Year'],
      "Npool" = Npool,
      "getsd" = FALSE,
      "savedir" = NULL,
      #run_model  =  FALSE,
      "newtonsteps" = 1,
      "test_fit" = FALSE,
      "working_dir" = test_path )

  # Comparisons
  Par1 = fit$parameter_estimates$par[names(fit$parameter_estimates$par)%in%c("ln_H_input","beta1_ct","beta1_ft","logkappa1","beta2_ct","beta2_ft","logkappa1","logSigmaM")]
  Par2 = opt$par[names(opt$par)%in%c("ln_H_input","beta1_ct","logkappa1","beta2_ct","logkappa1","logSigmaM")]
  expect_equal( as.vector(Par1), as.vector(Par2), tolerance=1e-3 )
  expect_equal( as.numeric(opt$objective), as.numeric(fit$parameter_estimates$objective)+1, tolerance=1e-3 )

  # Get expanded composition estimates and input sample sizes
  #Opt = TMBhelper::fit_tmb( obj=fit$tmb_list$Obj, startpar=fit$parameter_estimates$par, getsd=TRUE, lower=fit$tmb_list$Lower, upper=fit$tmb_list$Upper, newtonsteps=0 )
  #Index = plot_biomass_index( TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, DirName=test_path )
  #calculate_proportion( TmbData=fit$data_list, Index=Index, DirName=test_path, Year_Set=fit$year_labels )
})

