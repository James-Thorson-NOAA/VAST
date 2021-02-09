
# Tutorial: http://r-pkgs.had.co.nz/tests.html
# And example see: https://github.com/ss3sim/ss3sim/tree/master/tests/testthat
context("Testing examples")

# Eastern Bering Sea pollcok
test_that("Zero-inflated Poisson gives identical results to pscl::zeroinfl(.) ", {
  skip_on_travis()

  ## Simulate
  set.seed(101)
  n <- 50
  lambda <- 2
  zeroprob <- 0.2
  xyz <- data.frame(year = 1, Lon = runif(n, -12, -5), Lat = runif(n, 48, 52))
  m <- nrow(xyz)
  xyz$z <- rpois(n, lambda = 2)
  xyz$z <- xyz$z * rbinom(n, size=1, prob=1-zeroprob)

  ## Make settings
  settings <- make_settings(n_x = 10,
                            Region = "other",
                            ObsModel = c(7, 0),
                            purpose = "index")

  ## Run zero-inflated Poisson
  fit0 <- fit_model(settings = settings,
                   Lat_i= xyz$Lat,
                   Lon_i = xyz$Lon,
                   t_i = xyz$year,
                   c_i = rep(0,nrow(xyz)),
                   b_i = xyz$z,
                   a_i = rep(1, nrow(xyz)),
                   observations_LL = xyz[,c('Lat','Lon')],
                   ObsModel = c(7, 0),
                   Aniso = FALSE,
                   working_dir = multispecies_example_path,
                   FieldConfig = c(Omega1 = 0, Epsilon1 = 0,
                                   Omega2 = 0, Epsilon2 = 0))

  #
  fit0_pscl <- pscl::zeroinfl(z ~ 1, data = xyz, dist = "poisson")
  pscl_parhat = rbind( c(-1,1)*summary(fit0_pscl)$coef$zero[1,1:2], summary(fit0_pscl)$coef$count[1,1:2] )

  # Comparison with formula_factors
  expect_equal( as.numeric(summary(fit0$parameter_estimates$SD,"fixed")), as.numeric(pscl_parhat), tolerance=0.001 )
})

