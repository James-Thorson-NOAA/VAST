context("Testing combining acoustic and bottom trawl data")

test_that("Combining categories example is working ", {
  ## This is a stripped down version of the example at
  ##  https://github.com/James-Thorson-NOAA/VAST/wiki/Combine-acoustic-and-bottom-trawl-data
  ## to work with CI in github actions

  # Previously worked with CI, but not anymore
  skip_on_ci()
  skip_if(skip_local)

  # Disabled because it's using hte old EBS grid
  skip_if(TRUE)

  ## Read in prepared data which comes from a simulated example
  ## loosely conditioned on EBS pollock.

  ## !! THIS IS NOT REAL DATA !! ###
  data(acoustic_and_trawl, package = "VAST" )
  dat <- subset(acoustic_and_trawl, Year<2012)
  ## Make default settings for index standardization
  settings <- make_settings(
    n_x = 100,
    Version=Version_VAST,
    Region = "eastern_bering_sea",
    purpose = "index2",
    bias.correct = FALSE,
    fine_scale = TRUE,
    max_cells=Inf)
  settings$FieldConfig[2,2] <- 0 # turn off spatiotemporal LP2
  ## Setup the "combined" part of the model
  ## Associate depth strata (0 = <0.5m, 1 = 0.5-16m, 2 = >16m) with each
  ## observation (row). This is two columns b/c there are
  ## observations from at most 2 strata.
  c_iz <- matrix( c(1,2, 2,NA, 3,NA),
                 byrow = TRUE,
                 nrow = 3,
                 ncol = 2)[as.numeric(dat$Gear)] - 1

  ## Pulled from prevoius run and used dput
  par.mle <- c(0.0248120108120263, -0.626485433756577, 1.89822973227332, -1.42246167112512,
               -3.1797375137418, 1.48377620922032, -1.11812118230679, -2.90645503396005,
               1.83959007633375, -1.50781345548461, -4.39058347401738, 1.36484647915788,
               -1.69239396837296, -3.57146409043051, 2.00057539077216, 1.32373753739063,
               2.01084490099894, 1.91215326325085, 0.91380567507344, 0.714918993846101,
               0.509157821168581, -4.51106915786979, 4.83111649420249, 4.70704612274629,
               7.17579951218695, 4.84730143896306, 4.32995280146583, 6.2720447083131,
               3.61336547280719, 3.96945635517732, 6.51369500367214, 5.49121847155114,
               4.86820559346004, 6.95211186071455, 5.17633364166951, 2.10557108901664,
               2.44688675934284, 1.89110172145742, -4.41252477454976, 0.305806625171559,
               0.538372354396114, 0.496393515603893)

  ## Test model at MLE
  wd <- tempdir()
  fit <- fit_model(settings=settings, Lat_i=dat$Lat, Lon_i=dat$Lon,
                   t_i=dat$Year, c_i=c_iz, b_i=dat$Catch_KG,
                   a_i=dat$AreaSwept_km2, run_model=FALSE,
                   working_dir=paste0(wd, '/'))

  Obj <- fit$tmb_list$Obj
  Obj$env$beSilent()
  nll <- as.numeric(Obj$fn(par.mle))
  grads <- Obj$gr(par.mle)
  expect_equal(max(abs(grads)), 5.279877e-07)
  expect_equal(nll,32294.2448)
})
