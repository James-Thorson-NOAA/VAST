# VAST will often leave you in the subdirectory of the current run. Using HomeDir helps get you back where you started.
# Only do this once per R session, after you are in the your main working directory:

HomeDir <- getwd()

# =============================================

#Test run of single species spatial delta glmm
#Test, canary data; implimentation, Lingcod groundfish survey data
# Based on single-species example
# Revised by M. Haltuch, Feb 2017
# Revised by J. Wallace Mar 2017
# Revised by James Thorson April 2017

if (!any(installed.packages()[, 1] %in% "JRWToolBox"))
     devtools::install_github("John-R-Wallace/R-ToolBox")

if (!any(installed.packages()[, 1] %in% "VAST"))
    devtools::install_github("james-thorson/VAST")

if (!any(installed.packages()[, 1] %in% "pander"))
     install.packages("pander")

if (!any(installed.packages()[, 1] %in% "nwfscDeltaGLM"))
     devtools::install_github("nwfsc-assess/nwfscDeltaGLM")

require(TMB)
require(VAST)

# Yelloweye
data( WC_triennial_yelloweye, package="SpatialDeltaGLMM" )
Data_Set = WC_triennial_yelloweye

# Restrict to triennial
Data_Set = Data_Set[ which(Data_Set$SURVEY=="Tri.Shelf"), ]

# Exclude foreign hauls
data( AFSCforeign_hauls, package="nwfscDeltaGLM" )
Exclude = which(Data_Set$HAULJOIN %in% AFSCforeign_hauls$HAULJOIN)
if(length(Exclude)>0) Data_Set = Data_Set[-Exclude,]

# Reformat to same headers at Shelf-Slope (so that code is otherwise similar)
Data_Set = data.frame( "Year"=Data_Set[,'YEAR'], "Pass"=1.5, "Wt_kg"=Data_Set[,'WEIGHT'], "Latitude_dd"=rowMeans(Data_Set[,c('START_LATITUDE','END_LATITUDE')]), "Area_Swept_ha"=Data_Set[,'NET_WIDTH']*Data_Set[,'DISTANCE_FISHED']/10, "Longitude_dd"=rowMeans(Data_Set[,c('START_LONGITUDE','END_LONGITUDE')]), "Vessel"="none" )

# Look at the data by year and pass - showing 'NA's if any via JRWToolBox::Table function.
JRWToolBox::Table(Data_Set$Year, Data_Set$Pass)

# Versions of VAST you can use:
list.files(R.home(file.path("library", "VAST", "executables")))
# This gives the latest version available. (Up to v10_0_0 - then broken.)
(Version <- substr(list.files(R.home(file.path("library", "VAST", "executables")))[length(list.files(R.home(file.path("library", "VAST", "executables"))))], 1, 11))

#define the spatial resolution for the model, and whether to use a grid or mesh approximation
#mesh is default recommendation, number of knots need to be specified
#do not modify Kmeans set up
Method = c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km = 25     # Value only matters if Method="Grid"
n_x = 250  # Number of "knots" used when Method="Mesh"
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )   # Controls K-means algorithm to define location of knots when Method="Mesh"

# Model settings

# define whether to include spatial and spatio-temporal variation, whether its autocorrelated, and whether there's overdispersion
# field config - for both model components
# Omega- spatial variation
# Epsilon - temporal spatial variation
# review these settings
# if all field config settings are zero it is a fixed effects model
# RhoConfig - autocorrelation across time: defaults to zero, both annual intercepts (beta) and spatio-temporal (epsilon)
# OverdispersionConfig, vessel effects for both components of the model?
# settings can be on or off; 0,1
# obs model - distribution for errors and which model to run (e.g. default is delta model with standard link functions)
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1)
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Delta1"=0, "Delta2"=0)  # Turn on vessel-year effects for both components if using WCGBTS
ObsModel = c(2,0)

#outputs calculated after model runs, essentially reports to create
Options =  c("SD_site_density"=0, "SD_site_logdensity"=0, "Calculate_Range"=0, "Calculate_evenness"=0, "Calculate_effective_area"=0, "Calculate_Cov_SE"=0,
             'Calculate_Synchrony'=0, 'Calculate_Coherence'=0)

#strata limits, run model but then calculate area specific indices
  (strata.limits <- data.frame(
    'STRATA' = c("Coastwide","CA","OR","WA"),
    'north_border' = c(49.0, 42.0, 46.0, 49.0),
    'south_border' = c(32.0, 32.0, 42.0, 46.0),
    'shallow_border' = c(55, 55, 55, 55),
    'deep_border' = c(1280, 1280, 1280, 1280)
    ))

setwd(HomeDir)  # Make sure that the working directory is back where it started

#region that tells software which grid to use
Region = "California_current"
Domain = c("WCGBTS", "Triennial")[2]

#save files setting

# DateFile = paste0(getwd(),'/VAST_output/')  # Simple, but requires manually changing the directory to save different runs
(DateFile <- paste0(getwd(),'/VAST_output_', Sys.Date(), '_Yelloweye-Triennial_nx=', n_x, '_Domain=',Domain,'/')) # Change '_LingCod_nx=' for different runs, e.g. '_LingCod_Pass_nx=' for including pass
dir.create(DateFile)

#save all settings
# Record = ThorsonUtilities::bundlelist( c("Data_Set","Version","Method","grid_size_km","n_x","FieldConfig","RhoConfig","OverdispersionConfig","ObsModel","Kmeans_Config") )
# save( Record, file=file.path(DateFile,"Record.RData"))
# capture.output( Record, file=paste0(DateFile,"Record.txt"))

#set up data frame from data set
#creates data geostat...need this data format
# Vessel has a unique value for each boat-licence and calendar year (i.e., its a "Vessel-Year" effect)
Data_Geostat = data.frame(Catch_KG = Data_Set$Wt_kg, Year = Data_Set$Year, Vessel = paste(Data_Set$Vessel,Data_Set$Year,sep="_"),
             AreaSwept_km2 = Data_Set$Area_Swept_ha/100, Lat =Data_Set$Latitude_dd,
             Lon = Data_Set$Longitude_dd, Pass = Data_Set$Pass - 1.5)

#see data format
head(Data_Geostat)

# shows data being used, read this document
pander::pandoc.table( Data_Geostat[1:6,], digits=3 )

#extrapolation grid
if( Domain=="WCGBTS" ){
  Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits, surveyname='propInWCGBTS' )
}
if( Domain=="Triennial" ){
  Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits, surveyname='propInTriennial' )
}
#derived objects for spatio-temporal estiamtion
Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn( grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'],
                         Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]],
                         DirPath=DateFile, Save_Results=FALSE )

# Add knots to Data_Geostat
Data_Geostat = cbind( Data_Geostat, knot_i = Spatial_List$knot_i )
head(Data_Geostat)

#build model, this is where you could specify new covariates using Data_Fn...read more on this
# No Pass included
TmbData = VAST::Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "OverdispersionConfig"=OverdispersionConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel,
                   "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1,
                   "s_i"=Data_Geostat[,'knot_i']-1, "t_i"=Data_Geostat[,'Year'], "a_xl"=Spatial_List$a_xl, "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList,
                   "Method"=Spatial_List$Method, "Options"=Options )

################
# Do the estimation
################

# Build tmb object
TmbList = VAST::Build_TMB_Fn("TmbData"=TmbData, "RunDir"=DateFile, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]

# Run optimizer with Newton steps to improve convergence
Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE,
    newtonsteps=2, savedir=DateFile, bias.correct=TRUE )

# Create the report
Report = Obj$report()

# Save everything in object "Save" so that if you load it again, you can attach Save or not,
  # and know you haven't polluted your workspace
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file=paste0(DateFile,"Save.RData"))

# Check convergence via gradient (should be TRUE)
all( abs(Opt$diagnostics[,'final_gradient'])<1e-6 )
# Check convergence via Hessian (should be TRUE)
all( eigen(Opt$SD$cov.fixed)$values>0 )

setwd(HomeDir)

################
# Model output (some of the diagnostic plots are slow, so do the model ouptut first)
################

# Decide which years to plot
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

# Get region-specific settings for plots
MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( Region = Region, NN_Extrap = Spatial_List$PolygonList$NN_Extrap,  Extrapolation_List = Extrapolation_List )

#Plot Anisotropy
SpatialDeltaGLMM::PlotAniso_Fn( FileName=paste0(DateFile,"Aniso.png"), Report=Report, TmbData=TmbData )

# Annual density surface, use plot_set = 3 to start and then do plot_set=c(3:10) to see more output on a good working model
SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=3, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]],
        MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set,
        Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]],
        mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)

#index of abundance
Index = SpatialDeltaGLMM::PlotIndex_Fn( DirName=DateFile, TmbData=TmbData, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include,
  strata_names=strata.limits[,1], use_biascorr=TRUE )
pander::pandoc.table( Index$Table[,c("Year","Fleet","Estimate_metric_tons","SD_log","SD_mt")] )

#center of gravity / range exapnsion
SpatialDeltaGLMM::Plot_range_shifts(Report=Report, TmbData=TmbData, Sdreport=Opt[["SD"]], Znames=colnames(TmbData$Z_xm), PlotDir=DateFile, Year_Set=Year_Set)


################
# Make diagnostic plots
################

SpatialDeltaGLMM::Plot_data_and_knots(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=DateFile )

#convergence
pander::pandoc.table( Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')] )

# Plot encounter probability diagnostics p/a
Enc_prob = SpatialDeltaGLMM::Check_encounter_prob( Report=Report, Data_Geostat=Data_Geostat, DirName=DateFile)

# QQ plot
Q = SpatialDeltaGLMM::QQ_Fn( TmbData=TmbData, Report=Report, FileName_PP=paste0(DateFile,"Posterior_Predictive.jpg"), FileName_Phist=paste0(DateFile,"Posterior_Predictive-Histogram.jpg"),
                  FileName_QQ=paste0(DateFile,"Q-Q_plot.jpg"), FileName_Qhist=paste0(DateFile,"Q-Q_hist.jpg"))


# Residuals
SpatialDeltaGLMM::plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=TmbData, Report=Report, Q=Q, savedir=DateFile, MappingDetails=MapDetails_List[["MappingDetails"]],
           PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile,
           Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]],
           zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)

# Histogram of quantiles...should be a flat line for well behaved model; also can use the Q-Q plot

# Model selection, see example code, just run one model for now

setwd(HomeDir)

