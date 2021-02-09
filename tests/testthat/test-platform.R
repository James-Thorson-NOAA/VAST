

context("Testing cross platform and R version compatibility")

# Eastern Bering Sea pollcok
test_that("Eastern Bering Sea pollock is working ", {
  ## Prep really simple example using built-in data set, adapted
  ## from simple example on wiki
  example <- load_example( data_set="EBS_pollock" )
  dat <- example$sampling_data[example$sampling_data$Year==2012,]
  ## Make settings (turning off bias.correct to save time for example)
  settings <- make_settings(n_x=100, Region=example$Region,
                            purpose="index2",
                            strata.limits=example$strata.limits,
                            bias.correct=FALSE,
                            finescale=FALSE)

# Run model
fit = fit_model( settings = settings,
  Lat_i = example$sampling_data[,'Lat'],
  Lon_i = example$sampling_data[,'Lon'],
  t_i = example$sampling_data[,'Year'],
  c_i = rep(0,nrow(example$sampling_data)),
  b_i = example$sampling_data[,'Catch_KG'],
  a_i = example$sampling_data[,'AreaSwept_km2'],
  v_i = example$sampling_data[,'Vessel'] )
})

