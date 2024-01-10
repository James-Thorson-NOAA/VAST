
context("Testing examples")

# Eastern Bering Sea pollcok
test_that("Condition-and-density example is working ", {
  skip_on_ci()
  skip_if(skip_local)
  # Prepping
  test_path = file.path(multispecies_example_path,"Condition_and_density")
  load( file.path(test_path,"saved_estimates.RData") )
  setwd(test_path)

  # START CODE BLOCK FROM WIKI
  example = load_example( data_set = "goa_arrowtooth_condition_and_density" )

  # Format data
  b_i = ifelse( !is.na(example$sampling_data[,'cpue_kg_km2']),
    example$sampling_data[,'cpue_kg_km2'],
    example$sampling_data[,'weight_g'] )
  c_i = ifelse( !is.na(example$sampling_data[,'cpue_kg_km2']), 0, 1 )
  Q_i = ifelse(!is.na(example$sampling_data[,'cpue_kg_km2']),
    0, log(example$sampling_data[,'length_mm']) )

  # Make settings
  settings = make_settings( n_x = 250,
    Region = example$Region,
    purpose = "condition_and_density",
    bias.correct = FALSE,
    knot_method = "grid",
    Version=Version_VAST )
  settings$FieldConfig[c("Omega","Epsilon"),"Component_1"] = "IID"
  Expansion_cz = matrix( c( 0,0, 2,0 ), ncol=2, byrow=TRUE )
  settings$ObsModel = matrix( c(2,4, 1,4), ncol=2, byrow=TRUE )

  # Run model
  fit = fit_model( settings = settings,
    Lat_i = example$sampling_data[,'latitude'],
    Lon_i = example$sampling_data[,'longitude'],
    t_i = example$sampling_data[,'year'],
    c_i = c_i,
    b_i = b_i,
    a_i = rep(1, nrow(example$sampling_data)),
    Q_ik = matrix(Q_i, ncol=1),
    Expansion_cz = Expansion_cz,
    build_model = FALSE,
    backwards_compatible_kmeans=TRUE )

  #
  if( FALSE ){
    tapply( b_i, INDEX=list(example$sampling_data[,'year'],c_i), FUN=function(x){mean(x>0)} )
    DataList = fit$data_list
    TmbParams = fit$tmb_list$Parameters
    RhoConfig = c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0)
    Npool = 0
  }

  # Modify Map
  Map = fit$tmb_list$Map
    Map$lambda2_k = factor(NA)

  # Run model
  fit = fit_model( settings = settings,
    Lat_i = example$sampling_data[,'latitude'],
    Lon_i = example$sampling_data[,'longitude'],
    t_i = example$sampling_data[,'year'],
    c_i = c_i,
    b_i = b_i,
    a_i = rep(1, nrow(example$sampling_data)),
    Q_ik = matrix(Q_i, ncol=1),
    Expansion_cz = Expansion_cz,
    Map = Map,
  # END CODE BLOCK FROM WIKI
    #framework = "CppAD",
    getsd=FALSE,
    #run_model=FALSE,
    #newtonsteps=0,
    backwards_compatible_kmeans=TRUE )

  # Comparisons
  Par1 = fit$parameter_estimates$par[names(fit$parameter_estimates$par)%in%c("ln_H_input","logkappa1","logSigmaM")]
  Par2 = parameter_estimates$par[names(parameter_estimates$par)%in%c("ln_H_input","logkappa1","logSigmaM")]
  expect_equal( as.vector(Par1), as.vector(Par2), tolerance=1e-3 )
  expect_equal( parameter_estimates$objective, fit$parameter_estimates$objective, tolerance=1e-3 )
})

