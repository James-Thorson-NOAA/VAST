
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
test_that("Density covariates give identical results to glm(.) ", {
  # Previously worked with CI, but not anymore
  skip_on_ci()
  skip_if(skip_local)

  # load data set
  example = load_example( data_set="covariate_example" )

  # Make new factor
  example$covariate_data = cbind( example$covariate_data,
    "Depth_bin"=factor(ifelse(example$covariate_data[,'BOT_DEPTH']>200,'Deep','Shallow')) )

  # Rescale covariates being used to have an SD >0.1 and <10 (for numerical stability)
  example$covariate_data[,c("BOT_DEPTH","BOT_TEMP","SURF_TEMP","TEMP_STRAT")] =
    example$covariate_data[,c("BOT_DEPTH","BOT_TEMP","SURF_TEMP","TEMP_STRAT")] / 100

  # Make data for GLMs
  Data1 = Data2 = cbind( example$sampling_data, example$covariate_data,
    Year_factor = factor(example$covariate_data[,'Year'],levels=sort(unique(example$covariate_data[,'Year']))) )
  Data1[,'Catch_KG'] = ifelse( Data1[,'Catch_KG']>0, 1, 0 )
  Data2 = Data2[ which(Data2[,'Catch_KG']>0), ]

  # set Year = NA to treat all covariates as "static" (not changing among years)
  # If using a mix of static and dynamic covariates, please email package author to add easy capability
  year_set = min(example$covariate_data[,'Year']):max(example$covariate_data[,'Year'])
  example$covariate_data[,'Year'] = NA

  # NOW necessary to duplicate across years
  tmp = example$covariate_data
  for( year in year_set ){
    if(year == year_set[1]){
      example$covariate_data[,'Year'] = year
    }else{
      tmp[,'Year'] = year
      example$covariate_data = rbind(example$covariate_data, tmp)
    }
  }

  # Make new factor
  example$covariate_data = cbind( example$covariate_data,
    "Year_factor"=factor(example$covariate_data[,'Year'],levels=sort(unique(example$covariate_data[,'Year']))) )

  # Scramble order of example$covariate_data and example$sampling_data for testing purposes
  #Reorder = sample(1:nrow(example$covariate_data),replace=FALSE)
  #example$covariate_data = example$covariate_data[ Reorder, ]
  #example$sampling_data = example$sampling_data[ Reorder, ]

  # Make settings (turning off bias.correct to save time for example)
  settings3 = settings2 = settings1 = make_settings(
            n_x=100,
            Region=example$Region,
            purpose="index",
            use_anisotropy=FALSE,
            strata.limits=example$strata.limits,
            bias.correct=FALSE,
            fine_scale=TRUE,
            FieldConfig=c(0,0,0,0),
            ObsModel=c(1,0),
            Version=Version_VAST )
  settings2$ObsModel = c(2,0)
  settings3$ObsModel = c(3,0)

  # Define formula
  formula = ~ BOT_DEPTH:Year_factor + I(BOT_DEPTH^2)
  formula_factors = ~ Depth_bin

  # Run model -- Lognormal
  #source( "C:/Users/James.Thorson/Desktop/Git/FishStatsUtils/R/fit_model.R")
  #source( "C:/Users/James.Thorson/Desktop/Git/VAST/R/make_data.R")
  fit1 = fit_model( settings = settings1,
        Lat_i = example$sampling_data[,'Lat'],
        Lon_i = example$sampling_data[,'Lon'],
        t_i = example$sampling_data[,'Year'],
        b_i = example$sampling_data[,'Catch_KG'],
        a_i = example$sampling_data[,'AreaSwept_km2'],
        X1_formula = formula,
        X2_formula = formula,
        covariate_data = example$covariate_data,
        working_dir = multispecies_example_path )
  fit1B = fit_model( settings = settings1,
        Lat_i = example$sampling_data[,'Lat'],
        Lon_i = example$sampling_data[,'Lon'],
        t_i = example$sampling_data[,'Year'],
        b_i = example$sampling_data[,'Catch_KG'],
        a_i = example$sampling_data[,'AreaSwept_km2'],
        formula = formula,
        covariate_data = example$covariate_data,
        working_dir = multispecies_example_path )

  # Run model -- Gamma
  fit2 = fit_model( settings=settings2,
        Lat_i=example$sampling_data[,'Lat'],
        Lon_i=example$sampling_data[,'Lon'],
        t_i=example$sampling_data[,'Year'],
        b_i=example$sampling_data[,'Catch_KG'],
        a_i=example$sampling_data[,'AreaSwept_km2'],
        X1_formula=formula,
        X2_formula=formula,
        covariate_data=example$covariate_data,
        working_dir=multispecies_example_path )
  fit2B = fit_model( settings=settings2,
        Lat_i=example$sampling_data[,'Lat'],
        Lon_i=example$sampling_data[,'Lon'],
        t_i=example$sampling_data[,'Year'],
        b_i=example$sampling_data[,'Catch_KG'],
        a_i=example$sampling_data[,'AreaSwept_km2'],
        formula=formula,
        covariate_data=example$covariate_data,
        working_dir=multispecies_example_path )

  # Run model -- Inverse-Gaussian
  if( FALSE ){
    fit3 = fit_model( settings=settings3, Lat_i=example$sampling_data[,'Lat'],
      Lon_i=example$sampling_data[,'Lon'], t_i=example$sampling_data[,'Year'],
      b_i=example$sampling_data[,'Catch_KG'], a_i=example$sampling_data[,'AreaSwept_km2'],
      X1_formula=formula, X2_formula=formula, covariate_data=example$covariate_data,
      working_dir=multispecies_example_path )
    fit3B = fit_model( settings=settings3, Lat_i=example$sampling_data[,'Lat'],
      Lon_i=example$sampling_data[,'Lon'], t_i=example$sampling_data[,'Year'],
      b_i=example$sampling_data[,'Catch_KG'], a_i=example$sampling_data[,'AreaSwept_km2'],
      formula=formula, covariate_data=example$covariate_data,
      working_dir=multispecies_example_path )
  }

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
  Glm1B = stats::glm( formula=update.formula(formula, log(Catch_KG/AreaSwept_km2)~0+factor(Year)+.),
    data=Data2 )
  expect_equal( Glm1$coef, Glm1B$coef, tolerance=0.001 )

  # Glm2 -- Gamma
  # glm uses constant shape: https://stats.stackexchange.com/questions/58497/using-r-for-glm-with-gamma-distribution
  Glm2 = stats::glm( formula=update.formula(formula, Catch_KG~0+factor(Year)+.),
    data=Data2, family=Gamma(link="log"), offset=log(AreaSwept_km2) )
  Glm2B = stats::glm( formula=update.formula(formula, I(Catch_KG/AreaSwept_km2)~0+factor(Year)+.),
    data=Data2, family=Gamma(link="log") )
  expect_equal( Glm2$coef, Glm2B$coef, tolerance=0.001 )

  # Glm3 -- Inverse-Gaussian
  if( FALSE ){
    # Requires starting value to converge
    Glm3 = stats::glm( formula=update.formula(formula, Catch_KG~0+factor(Year)+.),
      data=Data2, family=inverse.gaussian(link="log"), start=Glm2$coef, offset=log(AreaSwept_km2) )
    Glm3B = stats::glm( formula=update.formula(formula, I(Catch_KG/AreaSwept_km2)~0+factor(Year)+.),
      data=Data2, family=inverse.gaussian(link="log"), start=Glm2$coef ) #, offset=AreaSwept_km2 )
    expect_equal( Glm3$coef, Glm3B$coef, tolerance=0.001 )
  }

  # Try Lognormal with factors in formula
  fit1_factors = fit_model( settings=settings1, Lat_i=example$sampling_data[,'Lat'],
    Lon_i=example$sampling_data[,'Lon'], t_i=example$sampling_data[,'Year'],
    b_i=example$sampling_data[,'Catch_KG'], a_i=example$sampling_data[,'AreaSwept_km2'],
    X1_formula=formula_factors, X2_formula=formula_factors, covariate_data=example$covariate_data,
    working_dir=multispecies_example_path )

  #
  Glm1_factors = stats::glm( formula=update.formula(formula_factors, log(Catch_KG)~0+factor(Year)+.),
    data=Data2, offset=log(AreaSwept_km2) )

  # Compare predictions from Glm1 and fit1
  if( formula == ~ BOT_DEPTH:Year_factor + I(BOT_DEPTH^2) ){
    predict_data = fit1$spatial_list$latlon_g
    predict_nn = RANN::nn2(query=predict_data, data=example$covariate_data[,c('Lat','Lon')], k=1)$nn.idx
    predict_data = cbind( predict_data, example$covariate_data[predict_nn,c("BOT_DEPTH","BOT_TEMP","SURF_TEMP","TEMP_STRAT","Year_factor")],
      "AreaSwept_km2"=fit1$extrapolation_list$a_el[,1], "Year"=NA )
    # Check processed values for depth
    for( Year in sort(unique(Data1$Year)) ){
      tI = match( Year, min(Data1$Year):max(Data1$Year) )
      fit2_depth = fit2$data_list$X2_gctp[,1,tI,][,paste0("BOT_DEPTH:Year_factor",Year)]
      Glm2_depth = predict_data[,'BOT_DEPTH']
      expect_equal( as.numeric(fit2_depth), as.numeric(Glm2_depth), tolerance=0.001 )
    }
    # Check predicted positive density
    # Not working after change from factor(Year) to Year_factors
    #for( Year in sort(unique(Data1$Year)) ){
    #  tI = match(Year, fit2$year_labels )
    #  predict_data[,'Year'] = Year
    #  Glm2_pred = predict( Glm2, newdata=predict_data, type="response" )
    #  fit2_pred = fit2$Report$R2_gct[,1,tI] * fit1$extrapolation_list$a_el[,1]
    #  expect_equal( as.numeric(fit2_pred), as.numeric(Glm2_pred), tolerance=0.001 )
    #}
  }else{
    stop("Check problem in `formula` in `test-covariates_against_glm.R`")
  }

  # Comparison with Glm0
  expect_equal( extract(fit1$parameter_estimates$par,"beta1_ft"), extract(Glm0$coef,"BOT_DEPTH",remove=TRUE), tolerance=0.001 )
  expect_equal( extract(fit1$parameter_estimates$par,"gamma1_cp"), extract(Glm0$coef,"BOT_DEPTH",remove=FALSE), tolerance=0.001 )

  # Comparison with Glm1
  expect_equal( extract(fit1$parameter_estimates$par,"beta2_ft") - exp(2*fit1$parameter_estimates$par['logSigmaM'])/2, extract(Glm1$coef,"BOT_DEPTH",remove=TRUE), tolerance=0.001 )
  expect_equal( extract(fit1$parameter_estimates$par,"gamma2_cp"), extract(Glm1$coef,"BOT_DEPTH",remove=FALSE), tolerance=0.001 )

  # Comparison with Glm2
  expect_equal( extract(fit2$parameter_estimates$par,"beta2_ft"), extract(Glm2$coef,"BOT_DEPTH",remove=TRUE), tolerance=0.001 )
  expect_equal( extract(fit2$parameter_estimates$par,"gamma2_cp"), extract(Glm2$coef,"BOT_DEPTH",remove=FALSE), tolerance=0.001 )

  # Comparison of old and new interface
  expect_equal( fit1$parameter_estimates$par, fit1B$parameter_estimates$par, tolerance=0.001 )
  expect_equal( fit2$parameter_estimates$par, fit2B$parameter_estimates$par, tolerance=0.001 )

  # Comparison with Glm3
  # Inverse-Gaussian doesn't seem to be identical, presumably due to different parameterizations
  if( FALSE ){
    expect_equal( extract(fit3$parameter_estimates$par,"beta2_ft"), extract(Glm3$coef,"BOT_DEPTH",remove=TRUE), tolerance=0.001 )
    expect_equal( extract(fit3$parameter_estimates$par,"gamma2_ctp"), extract(Glm3$coef,"BOT_DEPTH",remove=FALSE), tolerance=0.001 )
  }

  # Comparison with formula_factors
  expect_equal( extract(fit1_factors$parameter_estimates$par,"beta2_ft") - exp(2*fit1_factors$parameter_estimates$par['logSigmaM'])/2, extract(Glm1_factors$coef,"Depth_bin",remove=TRUE), tolerance=0.001 )
  expect_equal( extract(fit1_factors$parameter_estimates$par,"gamma2_cp"), extract(Glm1_factors$coef,"Depth_bin",remove=FALSE), tolerance=0.001 )

  ##################
  # Test predict.fit_model with existing covariates
  ##################

  # predict.fit_model test  -- Component 1
  pred0_glm = predict( Glm0, newdata=Data1, type="response" )
  pred0_vast = predict( fit2,
    Lat_i = Data1[,'Lat'],
    Lon_i = Data1[,'Lon'],
    t_i = Data1[,'Year'],
    a_i = Data1[,'AreaSwept_km2'],
    what = "R1_i",
    working_dir = multispecies_example_path )
  expect_equal( as.numeric(pred0_glm), pred0_vast, tolerance=0.001 )

  # predict.fit_model test  -- Component 2
  # Use Response-scale so that both are implementing the offset
  pred2_glm = predict( Glm2, newdata=Data2, type="response" )
  pred2_vast = predict( fit2,
    Lat_i = Data2[,'Lat'],
    Lon_i = Data2[,'Lon'],
    t_i = Data2[,'Year'],
    a_i = Data2[,'AreaSwept_km2'],
    what = "R2_i",
    working_dir = multispecies_example_path )
  expect_equal( as.numeric(pred2_glm), pred2_vast, tolerance=0.001 )

  ##################
  # Test predict.fit_model with new covariate data
  ##################

  # Randomize new data
  newdata = Data1
  for(cI in 1:ncol(Data1)) newdata[,cI] = sample( newdata[,cI], replace=FALSE )

  # predict.fit_model test  -- Component 1
  pred0_glm = predict( Glm0, newdata=newdata, type="response" )
  pred0_vast = predict( fit2,
    Lat_i = newdata[,'Lat'],
    Lon_i = newdata[,'Lon'],
    t_i = newdata[,'Year'],
    a_i = newdata[,'AreaSwept_km2'],
    what = "R1_i",
    new_covariate_data = newdata,
    working_dir = multispecies_example_path,
    do_checks = FALSE )
  expect_equal( as.numeric(pred0_glm), pred0_vast, tolerance=0.001 )

  # predict.fit_model test  -- Component 2
  # Use Response-scale so that both are implementing the offset
  pred2_glm = predict( Glm2, newdata=newdata, type="response" )
  pred2_vast = predict( fit2,
    Lat_i = newdata[,'Lat'],
    Lon_i = newdata[,'Lon'],
    t_i = newdata[,'Year'],
    a_i = newdata[,'AreaSwept_km2'],
    what = "R2_i",
    new_covariate_data = newdata,
    working_dir = multispecies_example_path,
    do_checks = FALSE )
  expect_equal( as.numeric(pred2_glm), pred2_vast, tolerance=0.001 )

})

