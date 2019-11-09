
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
test_that("Covariates give identical results to glm(.) ", {
  skip_on_travis()

  # load data set
  example = load_example( data_set="covariate_example" )

  # Scramble order of example$covariate_data and example$sampling_data for testing purposes
  Reorder = sample(1:nrow(example$covariate_data),replace=FALSE)
  example$covariate_data = example$covariate_data[ Reorder, ]
  example$sampling_data = example$sampling_data[ Reorder, ]

  # Make settings (turning off bias.correct to save time for example)
  settings = make_settings( n_x=100, Region=example$Region, purpose="index", use_anisotropy=FALSE,
    strata.limits=example$strata.limits, bias.correct=FALSE, fine_scale=TRUE,
    FieldConfig=c(0,0,0,0), ObsModel=c(1,0) )

  # Define formula
  formula = ~ BOT_DEPTH:factor(Year) + I(BOT_DEPTH^2)

  # set Year = NA to treat all covariates as "static" (not changing among years)
  # If using a mix of static and dynamic covariates, please email package author to add easy capability
  example$covariate_data[,'Year'] = NA

  # Rescale covariates being used to have an SD >0.1 and <10 (for numerical stability)
  example$covariate_data[,c("BOT_DEPTH","BOT_TEMP","SURF_TEMP","TEMP_STRAT")] =
    example$covariate_data[,c("BOT_DEPTH","BOT_TEMP","SURF_TEMP","TEMP_STRAT")] / 100

  # Run model
  fit = fit_model( "settings"=settings, "Lat_i"=example$sampling_data[,'Lat'],
    "Lon_i"=example$sampling_data[,'Lon'], "t_i"=example$sampling_data[,'Year'],
    "b_i"=example$sampling_data[,'Catch_KG'], "a_i"=example$sampling_data[,'AreaSwept_km2'],
    "formula"=formula, "covariate_data"=example$covariate_data, "working_dir"=multispecies_example_path )

  # Glm fits
  Data1 = Data2 = cbind( example$sampling_data, example$covariate_data )
  Data1[,'Catch_KG'] = ifelse( Data1[,'Catch_KG']>0, 1, 0 )
  Data2 = Data2[ which(Data2[,'Catch_KG']>0), ]

  # Function to facilitate comparison
  extract = function(vec, name, remove=FALSE){
    if(remove==FALSE) return( as.vector(vec[grep(name,names(vec))]) )
    if(remove==TRUE) return( as.vector(vec[-grep(name,names(vec))]) )
  }

  # Glm1
  Glm1 = stats::glm( formula=update.formula(formula, Catch_KG~0+factor(Year)+.),
    data=Data1, family="binomial" )

  # Glm2
  Glm2 = stats::glm( formula=update.formula(formula, log(Catch_KG/AreaSwept_km2)~0+factor(Year)+.),
    data=Data2 )

  # Comparison with Glm1
  expect_equal( extract(fit$parameter_estimates$par,"beta1_ft"), extract(Glm1$coef,"BOT_DEPTH",remove=TRUE), tolerance=0.0001 )
  expect_equal( extract(fit$parameter_estimates$par,"gamma1_ctp"), extract(Glm1$coef,"BOT_DEPTH",remove=FALSE), tolerance=0.0001 )

  # Comparison with Glm2
  expect_equal( extract(fit$parameter_estimates$par,"beta2_ft") - exp(2*fit$parameter_estimates$par['logSigmaM'])/2, extract(Glm2$coef,"BOT_DEPTH",remove=TRUE), tolerance=0.0001 )
  expect_equal( extract(fit$parameter_estimates$par,"gamma2_ctp"), extract(Glm2$coef,"BOT_DEPTH",remove=FALSE), tolerance=0.0001 )
})

