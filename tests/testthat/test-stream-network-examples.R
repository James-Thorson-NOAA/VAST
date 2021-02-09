context("Testing stream network example")

# Eastern Bering Sea pollcok
test_that("Stream network example is working ", {
  skip_on_ci()
  skip_if(skip_local)

  # Prepping
  test_path = file.path(multispecies_example_path,"Stream_network")
  load( file.path(test_path,"parameter_estimates.RData") )
  load( file.path(test_path, "settings.RData"))
  settings$max_cells = Inf
  settings$Version = FishStatsUtils::get_latest_version()
  attach(settings)
  on.exit( detach(settings) )

  # Run model
  data( "stream_network_eel_example", package="FishStatsUtils" )

  ## load each dataset object
  network <- stream_network_eel_example[["network"]]   ## network
  obs <- stream_network_eel_example[["observations"]]  ## sampling data
  hab <- stream_network_eel_example[["habitat"]]       ## habitat information 

  ## name stream network info
  Network_sz_LL = network
  Network_sz = network[,c('parent_s','child_s','dist_s')] 

  ## add a small value to presence observations
  present <- obs$data_value
  devs <- rnorm(length(present), 0, 0.01)
  present_new <- sapply(1:length(present), function(x) ifelse(present[x]==1, present[x]+devs[x], present[x]))
  obs$data_value <- present_new 

  ## setup dataset
  Data_Geostat <- data.frame( "Catch_KG" = obs$data_value,
                #"Year" = obs$year,  # CHANGING TO MAINTAIN COMPATIBILITY WITH HOW data.frame worked prior to R version 4.0.0
                "Year" = match(obs$year,sort(unique(obs$year))),
                 "Vessel" = obs$fishmethod,
                 "AreaSwept_km2" = obs$dist_i,
                 "Lat" = obs$lat,
                 "Lon" = obs$long,
                 "Pass" = 0,
                 "Knot" = obs$child_i,
                 "Category" = "Longfin_eels")

  fit = fit_model( settings=settings, 
                  Lat_i=Data_Geostat[,"Lat"], 
                  Lon_i=Data_Geostat[,"Lon"], 
                  t_i=as.numeric(Data_Geostat[,'Year']),
                  c_iz=rep(0,nrow(Data_Geostat)), 
                  b_i=Data_Geostat[,'Catch_KG'], 
                  a_i=Data_Geostat[,'AreaSwept_km2'], 
                  v_i=Data_Geostat[,'Vessel'], 
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"],
                                     "Lon"=Data_Geostat[,"Lon"],
                                     "Area_km2"=Data_Geostat[,"AreaSwept_km2"],
                                     "child_i"=Data_Geostat[,"Knot"]),
                  spatial_args=list(Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=1),
                  working_dir=multispecies_example_path,
                  test_fit=FALSE,
                  run_model=TRUE )

  # Comparisons
  Par1 = parameter_estimates$par[c("logkappa1","Beta_mean1_c")] # Not logSigmaM or beta2_ft, which depends on jittered values
  Par2 = fit$parameter_estimates$par[c("logkappa1","Beta_mean1_c")]
  expect_equal( as.vector(Par1), as.vector(Par2), tolerance=1e-3 )
  Par1 = abs(parameter_estimates$par[c("L_omega1_z","L_epsilon1_z")]) # Abs(L) to deal with sign-switching
  Par2 = abs(fit$parameter_estimates$par[c("L_omega1_z","L_epsilon1_z")])
  expect_equal( as.vector(Par1), as.vector(Par2), tolerance=1e-3 )
})
