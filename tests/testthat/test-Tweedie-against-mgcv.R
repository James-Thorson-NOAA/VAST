
# Tutorial: http://r-pkgs.had.co.nz/tests.html
# And example see: https://github.com/ss3sim/ss3sim/tree/master/tests/testthat
context("Testing examples")

# Tweedie distribution
test_that("Tweedie gives identical results to mgcv::gam(.) ", {
  # Previously worked with CI, but not anymore
  #skip_on_ci()
  skip_if(skip_local)
  library(mgcv)

  #library(tweedie)  # Installed from locally from tweedie_2.3.2.tar.gz here: https://cran.r-project.org/web/packages/tweedie/index.html
  # Simulate
  n_obs = 50
  mu = 2
  phi = 1
  p = 1.5

  # Simulate data
  set.seed(101)
  d_i = tweedie::rtweedie( n=n_obs, mu=mu,  phi=phi, p=p )
  Lon_i = runif( n=n_obs, min=-12, max=-5)
  Lat_i = runif( n=n_obs, min=48, max=52)
  t_i = rep(1, n_obs)

  # Fit using GAM
  mgcv_gam = gam( d_i ~ 1, family=tw )
  p = 1 + plogis(mgcv_gam$family$getTheta())

  # load data set
  example = load_example( data_set="EBS_pollock" )

  # Make settings (turning off bias.correct to save time for example)
  settings = make_settings( n_x=100,
    Region="other",
    purpose="index2",
    ObsModel = c(10,2),
    max_cells = Inf,
    Version = Version_VAST )
  settings$FieldConfig[c('Omega','Epsilon'),c('Component_1','Component_2')] = 0

  # Run model
  fit = fit_model( settings = settings,
    Lat_i = Lat_i,
    Lon_i = Lon_i,
    t_i = t_i,
    b_i = d_i,
    a_i = rep(1, n_obs),
    observations_LL = cbind("Lat"=Lat_i, "Lon"=Lon_i),
    working_dir = multispecies_example_path,
    getsd = FALSE,
    Use_REML = TRUE )
  p_hat = 1 + plogis( fit$ParHat$logSigmaM[1,1] )

  # Comparison -- MEAN
  expect_equal( as.numeric(fit$ParHat$beta1_ft+fit$ParHat$beta2_ft), as.numeric(mgcv_gam$coefficients['(Intercept)']), tolerance=0.001 )
  # Comparison -- theta a.k.a. p
  expect_equal( as.numeric(p_hat), as.numeric(p), tolerance=0.01 )
  # Comparison -- scale
  expect_equal( as.numeric(exp(fit$ParHat$beta1_ft)), as.numeric(mgcv_gam$scale), tolerance=0.05 )
})

# Eastern Bering Sea pollcok
test_that("Covariate effects when using a smoother gives identical results to mgcv::gam(.) ", {
  #skip_on_ci()
  skip_if(skip_local)
  library(mgcv)

  # Simulate
  data( EBS_pollock_data, package="FishStatsUtils" )
  EBS_pollock_data = EBS_pollock_data$sampling_data
  pollock_data = EBS_pollock_data[ which(EBS_pollock_data$waterTmpC > -999), ]

  # load data set
  example = load_example( data_set="EBS_pollock" )
  covariate_data = data.frame( "Lat"=pollock_data$lat, "Lon"=pollock_data$long, "Year"=pollock_data$year, "Temp"=pollock_data$waterTmpC)

  ##################
  # Tweedie
  ##################

  # Make settings (turning off bias.correct to save time for example)
  settings = make_settings( n_x = 100,
    Region = "Eastern_Bering_Sea",
    purpose = "index2",
    ObsModel = c(10,2),
    max_cells = Inf,
    fine_scale = FALSE,
    use_anisotropy = FALSE,
    Version = Version_VAST )
  settings$FieldConfig[c('Omega','Epsilon'),c('Component_1','Component_2')] = 0
  settings$FieldConfig['Omega','Component_2'] = 1
  settings$RhoConfig['Beta1'] = 3

  # Run model
  fit_tweedie = fit_model( settings = settings,
    Lat_i = pollock_data[,'lat'],
    Lon_i = pollock_data[,'long'],
    t_i = pollock_data[,'year'],
    b_i = pollock_data[,'catch'],
    a_i = rep(0.01, nrow(pollock_data) ),
    covariate_data = covariate_data,
    X2_formula = ~ Temp + I(Temp^2),
    working_dir = multispecies_example_path,
    getsd = FALSE,
    Use_REML = TRUE )
  p_hat = 1 + plogis( fit_tweedie$ParHat$logSigmaM[1,1] )

  # Fit using GAM
  x_i = ( fit_tweedie$spatial_list$A_is%*%(1:ncol(fit_tweedie$spatial_list$A_is)) )[,1]
  Data = cbind(pollock_data, fit_tweedie$spatial_list$loc_i, "Area_km2"=0.01, "pres"=ifelse(pollock_data$catch==0,0,1),
    Temp_i=fit_tweedie$data_list$X2_ip[,1], "east_i"=fit_tweedie$spatial_list$loc_x[x_i,"E_km"],
    "north_i"=fit_tweedie$spatial_list$loc_x[x_i,"N_km"] )
  mgcv_gam = gam( catch ~ 0 + factor(year) + Temp_i + I(Temp_i^2) + s(east_i,north_i,bs="gp"), family=tw,
    offset=log(Area_km2), data=Data, method="ML" )
  p = 1 + plogis(mgcv_gam$family$getTheta())
  #mgcv_gam$scale

  # slopes
  expect_equal( as.numeric(fit_tweedie$ParHat$gamma2_cp), as.numeric(summary(mgcv_gam)$p.coeff[c('Temp_i','I(Temp_i^2)')]), tolerance=0.01 )
  expect_equal( as.numeric(p_hat), as.numeric(p), tolerance=0.01 )

  ##################
  # Conventional delta-lognormal
  #  1.  Doesn't work well for 1st linear predictor for reasons I don't understand
  ##################

  # Make settings (turning off bias.correct to save time for example)
  settings = make_settings( n_x = 100,
    Region = "Eastern_Bering_Sea",
    purpose = "index2",
    ObsModel = c(1,0),
    max_cells = Inf,
    fine_scale = FALSE,
    use_anisotropy = FALSE,
    Version = Version_VAST )
  settings$FieldConfig['Epsilon',c('Component_1','Component_2')] = 0
  #settings$FieldConfig['Omega',c('Component_1','Component_2')] = 0

  # Run model
  fit = fit_model( settings = settings,
    Lat_i = pollock_data[,'lat'],
    Lon_i = pollock_data[,'long'],
    t_i = pollock_data[,'year'],
    b_i = pollock_data[,'catch'],
    a_i = rep(0.01, nrow(pollock_data) ),
    covariate_data = covariate_data,
    X1_formula = ~ Temp + I(Temp^2),
    X2_formula = ~ Temp + I(Temp^2),
    working_dir = multispecies_example_path,
    getsd = FALSE,
    Use_REML = TRUE )

  # Fit using GAM
  x_i = ( fit$spatial_list$A_is%*%(1:ncol(fit$spatial_list$A_is)) )[,1]
  Data = cbind(pollock_data, fit$spatial_list$loc_i, "Area_km2"=0.01, "pres"=ifelse(pollock_data$catch==0,0,1),
    Temp_i=fit$data_list$X1_ip[,1], "east_i"=fit$spatial_list$loc_x[x_i,"E_km"], "north_i"=fit$spatial_list$loc_x[x_i,"N_km"] )
  gam1 = gam( pres ~ 0 + factor(year) + Temp_i + I(Temp_i^2) + s(east_i,north_i,bs='gp'),
    family=binomial(link="logit"), data=Data, method="ML" )  # + te(lat,long,bs="gp") + te(lat,long,by=Year,bs="gp")
  gam2 = gam( log(catch) ~ 0 + factor(year)+ Temp_i + I(Temp_i^2) + s(east_i,north_i,bs='gp'),
    family=gaussian, offset=log(Area_km2), data=Data[which(Data$catch>0),], method="ML" )

  # slopes
  expect_equal( as.numeric(fit$ParHat$gamma1_cp), as.numeric(summary(gam1)$p.coeff[c('Temp_i','I(Temp_i^2)')]), tolerance=0.1 )
  expect_equal( as.numeric(fit$ParHat$gamma2_cp), as.numeric(summary(gam2)$p.coeff[c('Temp_i','I(Temp_i^2)')]), tolerance=0.1 )
})

