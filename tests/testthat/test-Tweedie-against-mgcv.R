
# Tutorial: http://r-pkgs.had.co.nz/tests.html
# And example see: https://github.com/ss3sim/ss3sim/tree/master/tests/testthat
context("Testing examples")

# Eastern Bering Sea pollcok
test_that("Tweedie gives identical results to mgcv::gam(.) ", {
  skip_on_travis()

  #library(mgcv)
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
  mgcv_gam = mgcv::gam( d_i ~ 1, family=tw )
  p = 1 + plogis(mgcv_gam$family$getTheta())
  #mgcv_gam$scale

  # Load package
  library(VAST)

  # load data set
  # see `?load_example` for list of stocks with example data
  # that are installed automatically with `FishStatsUtils`.
  example = load_example( data_set="EBS_pollock" )

  # Make settings (turning off bias.correct to save time for example)
  settings = make_settings( n_x=100,
    Region="other",
    purpose="index2",
    ObsModel = c(10,2),
    max_cells = Inf )
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
