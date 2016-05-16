

# Install TMB
# Must be installed from: https://github.com/kaskr/adcomp

# Install INLA
# Must be installed from: http://www.r-inla.org/download

# Install geostatistical delta-GLMM package
devtools::install_github("james-thorson/VAST", auth_token="550d63896f37bfef240ba2803f38b4a1387d5a86")
devtools::install_github("nwfsc-assess/geostatistical_delta-GLMM") # This is the developement version.  Please check GitHub for the latest release number.
devtools::install_github("james-thorson/utilities")

# setwd("C:/Users/James.Thorson/Desktop/Project_git/VAST/examples/")

# Load libraries
library(TMB)
library(INLA)
library(ThorsonUtilities)
library(VAST)

# This is where all runs will be located
DateFile = paste(getwd(),'/',Sys.Date(),'/',sep='')
  dir.create(DateFile)

###############
# Settings
###############

  Data_Set = c("Iceland_cod", "WCGBTS_canary", "GSL_american_plaice", "BC_pacific_cod", "EBS_pollock", "GOA_Pcod", "GOA_pollock", "GB_spring_haddock", "GB_fall_haddock", "SAWC_jacopever", "Sim")[5]
  Sim_Settings = list("Species_Set"=1:100, "Nyears"=10, "Nsamp_per_year"=600, "Depth_km"=-1, "Depth_km2"=-1, "Dist_sqrtkm"=0, "SigmaO1"=0.5, "SigmaO2"=0.5, "SigmaE1"=0.5, "SigmaE2"=0.5, "SigmaVY1"=0.05, "Sigma_VY2"=0.05, "Range1"=1000, "Range2"=500, "SigmaM"=1)
  Version = "geo_index_v3m"
  n_x = c(100, 250, 500, 1000, 2000)[3] # Number of stations
  FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1) # 1=Presence-absence; 2=Density given presence; #Epsilon=Spatio-temporal; #Omega=Spatial
  RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) # Structure for beta or epsilon over time: 0=None (default); 1=WhiteNoise; 2=RandomWalk; 3=Constant
  VesselConfig = c("Vessel"=0, "VesselYear"=0)
  ObsModel = c(1,0)  # 0=normal (log-link); 1=lognormal; 2=gamma; 4=ZANB; 5=ZINB; 11=lognormal-mixture; 12=gamma-mixture
  Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )     # Samples: Do K-means on trawl locs; Domain: Do K-means on extrapolation grid

  # Determine region
  Region = switch( Data_Set, "Iceland_cod"="Iceland", "WCGBTS_canary"="California_current", "GSL_american_plaice"="Gulf_of_St_Lawrence", "BC_pacific_cod"="British_Columbia", "EBS_pollock"="Eastern_Bering_Sea", "GOA_Pcod"="Gulf_of_Alaska", "GOA_pollock"="Gulf_of_Alaska", "GB_spring_haddock"="Northwest_Atlantic", "GB_fall_haddock"="Northwest_Atlantic", "SAWC_jacopever"="South_Africa", "Sim"="California_current")

# Decide on strata for use when calculating indices
  if( Data_Set %in% c("WCGBTS_canary","Sim")){
    # In this case, it will calculate a coastwide index, and also a separate index for each state (although the state lines are approximate)
    strata.limits <- data.frame(
      'STRATA' = c("Coastwide","CA","OR","WA"),
      'north_border' = c(49.0, 42.0, 46.0, 49.0),
      'south_border' = c(32.0, 32.0, 42.0, 46.0),
      'shallow_border' = c(55, 55, 55, 55),
      'deep_border' = c(1280, 1280, 1280, 1280)
    )
    # Override default settings for vessels
    VesselConfig = c("Vessel"=0, "VesselYear"=1)
  }
  if( Data_Set %in% c("GSL_american_plaice")){
    strata.limits <- data.frame('STRATA'="All_areas")
  }
  if( Data_Set %in% c("BC_pacific_cod")){
    # In this case, will not restrict the extrapolation domain at all while calculating an index
    strata.limits <- data.frame( 'STRATA'="All_areas")
  }
  if( Data_Set %in% c("EBS_pollock")){
    # In this case, will not restrict the extrapolation domain at all while calculating an index
    strata.limits <- data.frame( 'STRATA'="All_areas")
  }
  if( Data_Set %in% c("GOA_Pcod","GOA_pollock")){
    # In this case, will calculating an unrestricted index and a separate index restricted to west of -140W
    strata.limits <- data.frame(
      'STRATA' = c("All_areas", "west_of_140W"),
      'west_border' = c(-Inf, -Inf),
      'east_border' = c(Inf, -140)
    )
  }
  if( Data_Set %in% c("GB_spring_haddock","GB_fall_haddock")){
    # For NEFSC indices, strata must be specified as a named list of area codes
    strata.limits = list( 'Georges_Bank'=c(1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200, 1210, 1220, 1230, 1240, 1250, 1290, 1300) )
  }
  if( Data_Set %in% c("SAWC_jacopever")){
    strata.limits = data.frame( 'STRATA'="All_areas" )
  }
  if( Data_Set %in% c("Iceland_cod")){
    strata.limits = data.frame( 'STRATA'="All_areas" )
    # Turn off all spatial, temporal, and spatio-temporal variation in probability of occurrence, because they occur almost everywhere
    FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=1, "Epsilon2"=1)
    RhoConfig = c("Beta1"=3, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) # 0=Off; 1=WhiteNoise; 2=RandomWalk; 3=Constant
  }

  # Save options for future records
  Record = bundlelist( c("Data_Set","Sim_Settings","Version","n_x","FieldConfig","RhoConfig","VesselConfig","ObsModel","Kmeans_Config") )
  capture.output( Record, file=paste0(DateFile,"Record.txt"))

################
# Prepare data
# (THIS WILL VARY FOR DIFFERENT DATA SETS) 
################

# Read or simulate trawl data
  if(Data_Set=="WCGBTS_canary"){
    data( WCGBTS_Canary_example, package="SpatialDeltaGLMM" )
    Data_Geostat = data.frame( "Catch_KG"=WCGBTS_Canary_example[,'HAUL_WT_KG'], "Year"=as.numeric(sapply(WCGBTS_Canary_example[,'PROJECT_CYCLE'],FUN=function(Char){strsplit(as.character(Char)," ")[[1]][2]})), "Vessel"=WCGBTS_Canary_example[,"VESSEL"], "AreaSwept_km2"=WCGBTS_Canary_example[,"AREA_SWEPT_HA"]/1e2, "Lat"=WCGBTS_Canary_example[,'BEST_LAT_DD'], "Lon"=WCGBTS_Canary_example[,'BEST_LON_DD'], "Pass"=WCGBTS_Canary_example[,'PASS']-1.5)
  }
  if( Data_Set %in% c("BC_pacific_cod")){
    data( BC_pacific_cod_example, package="SpatialDeltaGLMM" )
    Data_Geostat = data.frame( "Catch_KG"=BC_pacific_cod_example[,'PCOD_WEIGHT'], "Year"=BC_pacific_cod_example[,'Year'], "Vessel"="missing", "AreaSwept_km2"=BC_pacific_cod_example[,'TOW.LENGTH..KM.']/100, "Lat"=BC_pacific_cod_example[,'LAT'], "Lon"=BC_pacific_cod_example[,'LON'], "Pass"=0)
    Data_Geostat = na.omit( Data_Geostat )
    Data_Geostat$Year = as.numeric( factor(Data_Geostat$Year))
  }
  if( Data_Set %in% c("GSL_american_plaice")){
    data( GSL_american_plaice, package="SpatialDeltaGLMM" )
    Print_Message( "GSL_american_plaice" )
    Data_Geostat = data.frame( "Year"=GSL_american_plaice[,'year'], "Lat"=GSL_american_plaice[,'latitude'], "Lon"=GSL_american_plaice[,'longitude'], "Vessel"="missing", "AreaSwept_km2"=GSL_american_plaice[,'swept'], "Catch_KG"=GSL_american_plaice[,'biomass']*GSL_american_plaice[,'vstd'] )
  }
  if(Data_Set=="EBS_pollock"){
    data( EBS_pollock_data, package="SpatialDeltaGLMM" )
    Data_Geostat = data.frame( "Catch_KG"=EBS_pollock_data[,'catch'], "Year"=EBS_pollock_data[,'year'], "Vessel"="missing", "AreaSwept_km2"=0.01, "Lat"=EBS_pollock_data[,'lat'], "Lon"=EBS_pollock_data[,'long'], "Pass"=0)
  }
  if(Data_Set=="GOA_Pcod"){
    data( GOA_pacific_cod, package="SpatialDeltaGLMM" )
    Data_Geostat = data.frame( "Catch_KG"=GOA_pacific_cod[,'catch'], "Year"=GOA_pacific_cod[,'year'], "Vessel"="missing", "AreaSwept_km2"=0.01, "Lat"=GOA_pacific_cod[,'lat'], "Lon"=GOA_pacific_cod[,'lon'], "Pass"=0)
    # Rename years and keep track of correspondance (for computational speed, given that there's missing years)
    Data_Geostat$Year = as.numeric( factor(Data_Geostat$Year))
  }
  if(Data_Set=="GOA_pollock"){
    data( GOA_walleye_pollock, package="SpatialDeltaGLMM" )
    Data_Geostat = data.frame( "Catch_KG"=GOA_walleye_pollock[,'catch'], "Year"=GOA_walleye_pollock[,'year'], "Vessel"="missing", "AreaSwept_km2"=0.01, "Lat"=GOA_walleye_pollock[,'lat'], "Lon"=GOA_walleye_pollock[,'lon'], "Pass"=0)
    # Rename years and keep track of correspondance (for computational speed, given that there's missing years)
    Data_Geostat$Year = as.numeric( factor(Data_Geostat$Year))
  }
  if( Data_Set=="GB_spring_haddock"){
    data( georges_bank_haddock_spring, package="SpatialDeltaGLMM" )         # standardized area swept = 0.0112 nm^2 = 0.0112*1.852^2 km^2
    Print_Message( "GB_haddock" )
    Data_Geostat = data.frame( "Catch_KG"=georges_bank_haddock_spring[,'CATCH_WT_CAL'], "Year"=georges_bank_haddock_spring[,'YEAR'], "Vessel"="missing", "AreaSwept_km2"=0.0112*1.852^2, "Lat"=georges_bank_haddock_spring[,'LATITUDE'], "Lon"=georges_bank_haddock_spring[,'LONGITUDE'])
  }
  if( Data_Set=="GB_fall_haddock"){
    data( georges_bank_haddock_fall, package="SpatialDeltaGLMM" )         # standardized area swept = 0.0112 nm^2 = 0.0112*1.852^2 km^2
    Print_Message( "GB_haddock" )
    Data_Geostat = data.frame( "Catch_KG"=georges_bank_haddock_fall[,'CATCH_WT_CAL'], "Year"=georges_bank_haddock_fall[,'YEAR'], "Vessel"="missing", "AreaSwept_km2"=0.0112*1.852^2, "Lat"=georges_bank_haddock_fall[,'LATITUDE'], "Lon"=georges_bank_haddock_fall[,'LONGITUDE'])
  }
  if( Data_Set=="SAWC_jacopever"){
    data( south_africa_westcoast_jacopever, package="SpatialDeltaGLMM" )         # standardized area swept = 0.0112 nm^2 = 0.0112*1.852^2 km^2
    #Data = read.csv( paste0(getwd(),"/../../examples/archive of data inputs for creation of grid files/South Africa/SAWC_geodata.csv") )
    Data_Geostat = data.frame( "Catch_KG"=south_africa_westcoast_jacopever[,'HELDAC'], "Year"=south_africa_westcoast_jacopever[,'Year'], "Vessel"="missing", "AreaSwept_km2"=south_africa_westcoast_jacopever[,'area_swept_nm2']*1.852^2, "Lat"=south_africa_westcoast_jacopever[,'cen_lat'], "Lon"=south_africa_westcoast_jacopever[,'cen_long'])
    Data_Geostat$Year = as.numeric( factor(Data_Geostat$Year))
  }
  if(Data_Set=="Sim"){
    Sim_DataSet = Geostat_Sim(Sim_Settings=Sim_Settings, Extrapolation_List=Extrapolation_List, MakePlot=TRUE)
    Data_Geostat = Sim_DataSet[['Data_Geostat']]
    True_Index = Sim_DataSet[['True_Index']]
  }
  if( Data_Set %in% c("Iceland_cod")){
    # WARNING:  This data set has not undergone much evaluation for spatio-temporal analysis
    data( iceland_cod, package="SpatialDeltaGLMM" )
    Data_Geostat = data.frame( "Catch_KG"=iceland_cod[,'Catch_b'], "Year"=iceland_cod[,'year'], "Vessel"=1, "AreaSwept_km2"=iceland_cod[,'towlength'], "Lat"=iceland_cod[,'lat1'], "Lon"=iceland_cod[,'lon1'])
    Data_Geostat = na.omit( Data_Geostat )
  }

# Get extrapolation data
  if( Region == "California_current" ){
    Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits )
  }
  if( Region == "British_Columbia" ){
    Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits, strata_to_use=c("HS","QCS") )
  }
  if( Region == "Eastern_Bering_Sea" ){
    Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits )
  }
  if( Region == "Gulf_of_Alaska" ){
    Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits )
  }
  if( Region == "Northwest_Atlantic" ){
    Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits )
  }
  if( Region == "South_Africa" ){
    Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits, region="west_coast" )
  }
  if( Region == "Iceland" ){
    Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits, observations_LL=Data_Geostat[,c('Lat','Lon')], maximum_distance_from_sample=15 )
  }
  if( Region == "Gulf_of_St_Lawrence" ){
    Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits )
  }

# Convert to an Eastings-Northings coordinate system
  tmpUTM = SpatialDeltaGLMM::Convert_LL_to_UTM_Fn( Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], zone=Extrapolation_List[["zone"]] )                                                         #$
  Data_Geostat = cbind( Data_Geostat, 'E_km'=tmpUTM[,'X'], 'N_km'=tmpUTM[,'Y'])

# Calculate k-means centroids (but only once for all species)
  Kmeans = SpatialDeltaGLMM::Calc_Kmeans(n_x=n_x, loc_orig=Data_Geostat[,c("E_km", "N_km")], randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile)
  loc_x = Kmeans$centers
  Data_Geostat = cbind( Data_Geostat, "knot_i"=Kmeans$cluster )

# Calc design matrix and areas
  PolygonList = SpatialDeltaGLMM::Calc_Polygon_Areas_and_Polygons_Fn( loc_x=loc_x, Data_Extrap=Extrapolation_List[["Data_Extrap"]], Covariates=c("none"), a_el=Extrapolation_List[["a_el"]])
  a_xl = PolygonList[["a_xl"]]

# Make mesh and info for anisotropy
  MeshList = SpatialDeltaGLMM::Calc_Anisotropic_Mesh(loc_x=loc_x)

################
# Make and Run TMB model
# (THIS WILL BE SIMILAR FOR EVERY DATA SET) 
################

  # Make TMB data list
  TmbData = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel, "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, "s_i"=Data_Geostat[,'knot_i']-1, "t_i"=Data_Geostat[,'Year'], "a_xl"=a_xl, "MeshList"=MeshList )

  # Make TMB object
  TmbList = Build_TMB_Fn("TmbData"=TmbData, "RunDir"=DateFile, "Version"=Version, "RhoConfig"=RhoConfig, "VesselConfig"=VesselConfig, "loc_x"=loc_x)
  Obj = TmbList[["Obj"]]

  # Run first time -- marginal likelihood
  Start_time = Sys.time()
  Obj$fn(Obj$par)
  # Run first time -- gradient with respect to fixed effects
  Obj$gr(Obj$par)

  # Run model
  for(i in 1:2) Opt = nlminb(start=Obj$env$last.par.best[-Obj$env$random], objective=Obj$fn, gradient=Obj$gr, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], control=list(eval.max=1e4, iter.max=1e4, trace=1))  # , rel.tol=1e-20
  Opt[["final_diagnostics"]] = data.frame( "Name"=names(Opt$par), "Lwr"=TmbList[["Lower"]], "Est"=Opt$par, "Upr"=TmbList[["Upper"]], "Gradient"=Obj$gr(Opt$par) )
  Opt[["total_time_to_run"]] = Sys.time() - Start_time
  Opt[["number_of_coefficients"]] = c("Total"=length(unlist(Obj$env$parameters)), "Fixed"=length(Obj$par), "Random"=length(unlist(Obj$env$parameters))-length(Obj$par) )
  capture.output( Opt, file=paste0(DateFile,"Opt.txt"))
    
  # Reports
  Report = Obj$report()                                      
  Sdreport = sdreport(Obj, bias.correct=TRUE)
  
  # Save stuff
  Save = list("Opt"=Opt, "Report"=Report, "Sdreport"=Sdreport, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
  save(Save, file=paste0(DateFile,"Save.RData"))
  capture.output( Opt, file=paste0(DateFile,"Opt.txt"))
  capture.output( Sdreport, file=paste0(DateFile,"Sdreport.txt"))

################
# Make diagnostic plots
################

  # Plot Anisotropy  
  if( TmbData$Options_vec['Aniso']==1 ){
    PlotAniso_Fn( FileName=paste0(DateFile,"Aniso.png"), Report=Report )
  }

  # Plot surface
  Dim = c( "Nrow"=ceiling(sqrt(TmbData$n_t)), "Ncol"=ceiling(TmbData$n_t/ceiling(sqrt(TmbData$n_t))) )
  par( mfrow=Dim )
  MapDetails_List = MapDetails_Fn( "Region"=Region, "NN_Extrap"=PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
  PlotResultsOnMap_Fn(plot_set=1:3, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=paste0(DateFile,"Field_"), Year_Set=sort(unique(Data_Geostat[,'Year'])), Rotate=MapDetails_List[["Rotate"]], mfrow=Dim, mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), Cex=MapDetails_List[["Cex"]])
                                                                                                                           
  # Plot index
  PlotIndex_Fn( DirName=DateFile, TmbData=TmbData, Sdreport=Sdreport, Year_Set=sort(unique(Data_Geostat[,'Year'])), strata_names=strata.limits[,1], use_biascorr=TRUE )

  # Positive catch rate Q-Q plot
  Q = QQ_Fn( TmbData=TmbData, Report=Report, FileName_PP=paste0(DateFile,"Posterior_Predictive.jpg"), FileName_Phist=paste0(DateFile,"Posterior_Predictive-Histogram.jpg"), FileName_QQ=paste0(DateFile,"Q-Q_plot.jpg"), FileName_Qhist=paste0(DateFile,"Q-Q_hist.jpg"))

  # Plot center of gravity
  Plot_range_shifts(Sdreport=Sdreport, Report=Report, TmbData=TmbData, Znames=colnames(TmbData$Z_xm), FileName_COG=paste0(DateFile,"center_of_gravity.png"))

  # Vessel effects
  #Return = Vessel_Fn(TmbData=TmbData, Sdreport=Sdreport, FileName_VYplot=paste0(DateFile,"VY-effect.jpg"))


  