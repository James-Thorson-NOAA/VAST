

context("Testing cross platform and R version compatibility")

# Eastern Bering Sea pollcok
test_that("Eastern Bering Sea pollock is working ", {
  ## Prep really simple example using built-in data set, adapted
  ## from simple example on wiki
  example <- load_example( data_set="EBS_pollock" )
  dat <- example$sampling_data[example$sampling_data$Year==2012,]
  settings <- make_settings(n_x=100, Region=example$Region,
                            purpose="index2",
                            strata.limits=example$strata.limits,
                            Version=Version,
                            bias.correct=FALSE,
                            fine_scale=FALSE, max_cells=Inf)
  settings$FieldConfig[1:2, 1:2] <- 0
  ## Run model
  wd <- tempdir()
  fit <- fit_model(settings=settings,
                   working_dir=paste0(wd, '/'),
                   Lat_i=dat$Lat,
                   Lon_i=dat$Lon,
                   t_i=dat$Year,
                   c_i=rep(0,nrow(dat)),
                   b_i=dat$Catch_KG,
                   a_i=dat$AreaSwept_km2,
                   v_i=dat$Vessel)
})

