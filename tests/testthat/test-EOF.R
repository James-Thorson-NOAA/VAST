

# Tutorial: http://r-pkgs.had.co.nz/tests.html
# And example see: https://github.com/ss3sim/ss3sim/tree/master/tests/testthat
context("Testing examples")

# Eastern Bering Sea pollcok
test_that("EOF is working ", {
  # Previously worked with CI, but not anymore
  #skip_on_ci()
  skip_if(skip_local)

  test_path = file.path(multispecies_example_path,"EOF")
  load( file.path(test_path,"parameter_estimates.RData") )
  #file.copy( from=paste0(test_path,"/Kmeans_extrapolation-2000.RData"), to=paste0(test_path,"/Kmeans_extrapolation-2000.RData") )
  #file.copy( from=paste0(test_path,"/Kmeans_knots-50.RData"), to=paste0(test_path,"/Kmeans_knots-50.RData") )

  # load data set
  example = load_example( data_set="five_species_ordination" )
  which_rows = which( example$sampling_data[,'species_number'] %in% c(1,2) &
                      example$sampling_data[,'Year'] %in% 2006:2015 )

  # Make settings:
  # including modifications from default settings to match
  # analysis in original paper
  settings = make_settings( n_x=50,
    Region=example$Region,
    purpose="EOF3",
    n_categories=2,
    ObsModel=c(1,1),
    RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0) )
  #settings$Version = "VAST_v14_0_0"

  # Run model (including settings to speed up run)
  #dyn.unload( paste0(TmbDir = system.file("executables", package = "VAST"), "/", settings$Version) )
  fit = fit_model( settings = settings,
    Lat_i = example$sampling_data[which_rows,'Lat'],
    Lon_i = example$sampling_data[which_rows,'Lon'],
    t_i = example$sampling_data[which_rows,'Year'],
    c_i = example$sampling_data[which_rows,'species_number']-1,
    b_i = example$sampling_data[which_rows,'Catch_KG'],
    a_i = example$sampling_data[which_rows,'AreaSwept_km2'],
    #Parameters  =  fit$ParHat,
    #run_model  =  TRUE,
    newtonsteps = 0,
    getsd = FALSE,
    Use_REML = TRUE,
    working_dir = test_path )

  # Comparisons -- Use abs(.) to avoid label switching
  Par1 = fit$parameter_estimates$par[names(fit$parameter_estimates$par)%in%c("ln_H_input","L_omega1_z","Ltime_epsilon1_z","logkappa1","logSigmaM")]
  Par2 = parameter_estimates$par[names(parameter_estimates$par)%in%c("ln_H_input","L_omega1_z","Ltime_epsilon1_z","logkappa1","logSigmaM")]
  expect_equal( abs(as.vector(Par1)), abs(as.vector(Par2)), tolerance=1e-3 )
  expect_equal( as.numeric(parameter_estimates$objective), as.numeric(fit$parameter_estimates$objective), tolerance=1e-3 )

})

