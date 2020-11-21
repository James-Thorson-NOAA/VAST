
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
test_that("Catchability covariates give identical results to glm(.) ", {
  skip_on_travis()

  # load data set
  example = load_example( data_set="covariate_example" )

  # Make new factor
  example$covariate_data = cbind( example$covariate_data, "Depth_bin"=factor(ifelse(example$covariate_data[,'BOT_DEPTH']>200,'Deep','Shallow')) )

  # Scramble order of example$covariate_data and example$sampling_data for testing purposes
  Reorder = sample(1:nrow(example$covariate_data),replace=FALSE)
  example$covariate_data = example$covariate_data[ Reorder, ]
  example$sampling_data = example$sampling_data[ Reorder, ]

  # Make settings (turning off bias.correct to save time for example)
  settings1 = make_settings( n_x=100, Region=example$Region, purpose="index2",
    use_anisotropy=FALSE, strata.limits=example$strata.limits, bias.correct=FALSE, fine_scale=TRUE,
    FieldConfig=c(0,0,0,0), ObsModel=c(1,0) )

  # Define formula
  formula = ~ factor(Depth_bin)
  Q_ik = matrix( ifelse(example$covariate_data[,'Depth_bin']=="Deep",0,1) )

  # Run model -- Lognormal
  fit1 = fit_model( settings=settings1, Lat_i=example$sampling_data[,'Lat'],
    Lon_i=example$sampling_data[,'Lon'], t_i=example$sampling_data[,'Year'],
    b_i=example$sampling_data[,'Catch_KG'], a_i=example$sampling_data[,'AreaSwept_km2'],
    Q1_formula=formula, Q2_formula=formula, catchability_data=example$covariate_data,
    working_dir=multispecies_example_path )

  # Define via Q_ik (attempt at eliminating Warning by unloading previous run fails)
  # suppressWarnings(dyn.unload( file.path(multispecies_example_path,TMB::dynlib(settings1$Version)) ))
  fit1B = fit_model( settings=settings1, Lat_i=example$sampling_data[,'Lat'],
    Lon_i=example$sampling_data[,'Lon'], t_i=example$sampling_data[,'Year'],
    b_i=example$sampling_data[,'Catch_KG'], a_i=example$sampling_data[,'AreaSwept_km2'],
    Q_ik=Q_ik, working_dir=multispecies_example_path )

  # Glm fits
  Data1 = Data2 = cbind( example$sampling_data, example$covariate_data )
  Data1[,'Catch_KG'] = ifelse( Data1[,'Catch_KG']>0, 1, 0 )
  Data2 = Data2[ which(Data2[,'Catch_KG']>0), ]

  # Function to facilitate comparison
  extract = function(vec, name, remove=FALSE){
    if(remove==FALSE) return( as.vector(vec[grep(name,names(vec))]) )
    if(remove==TRUE) return( as.vector(vec[-grep(name,names(vec))]) )
  }

  # Glm0 -- Bernoulli
  Glm0 = stats::glm( formula=update.formula(formula, Catch_KG~0+factor(Year)+.),
    data=Data1, family="binomial" )

  # Glm1 -- Lognormal
  Glm1 = stats::glm( formula=update.formula(formula, log(Catch_KG)~0+factor(Year)+.),
    data=Data2, offset=log(AreaSwept_km2) )

  # Comparison with Glm0
  expect_equal( extract(fit1$parameter_estimates$par,"beta1_ft"), extract(Glm0$coef,"Depth_bin",remove=TRUE), tolerance=0.001 )
  expect_equal( extract(fit1$parameter_estimates$par,"lambda1_k"), extract(Glm0$coef,"Depth_bin",remove=FALSE), tolerance=0.001 )

  # Comparison with Glm1
  expect_equal( extract(fit1$parameter_estimates$par,"beta2_ft") - exp(2*fit1$parameter_estimates$par['logSigmaM'])/2, extract(Glm1$coef,"Depth_bin",remove=TRUE), tolerance=0.001 )
  expect_equal( extract(fit1$parameter_estimates$par,"lambda2_k"), extract(Glm1$coef,"Depth_bin",remove=FALSE), tolerance=0.001 )

  # Comparison between fit1 and fit1B
  expect_equal( fit1$parameter_estimates$par, fit1B$parameter_estimates$par, tolerance=0.001 )

})

