
context("Testing examples")

# Eastern Bering Sea pollcok
test_that("Male lingcod compositional expansion is working ", {
  skip_on_travis()
  #test_path = "C:/Users/James.Thorson/Desktop/UW Hideaway/AFSC/2018-08 -- Lingcod testthat/"
  test_path = file.path(multispecies_example_path,"Lingcod_comp_expansion")
  load( file=file.path(test_path,"Data_Geostat.RData") )
  load( file.path(test_path,"opt.RData") )
  load( file.path(test_path,"Record.RData") )
  attach(Record)
  on.exit( detach(Record) )
  # Run model
  #data( EBS_pollock_data, package="FishStatsUtils" )
  Extrapolation_List = make_extrapolation_info( Region="California_current", strata.limits=strata.limits )
  Spatial_List = make_spatial_info( grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon=Data_Geostat[,'BEST_LON_DD'], Lat=Data_Geostat[,'BEST_LAT_DD'], Extrapolation_List=Extrapolation_List, DirPath=test_path )
  Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )
  TmbData = VAST::Data_Fn("CheckForErrors"=FALSE, "Version"=Version_VAST, "OverdispersionConfig"==OverdispersionConfig, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel_Est, "c_i"=as.numeric(Data_Geostat[,'Length_bin'])-1, "b_i"=Data_Geostat[,'First_stage_expanded_numbers'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=rep(0,nrow(Data_Geostat)), "s_i"=Data_Geostat[,'knot_i']-1, "t_i"=Data_Geostat[,'Year'], "a_xl"=Spatial_List$a_xl, "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList, "Method"=Spatial_List$Method, "Aniso"=Aniso )
  # Switch missing categories to NA
  EncProb_tc = tapply( TmbData$b_i, INDEX=list(TmbData$t_i,TmbData$c_i), FUN=function(vec){mean(vec>0)} )
  TmbData$b_i = ifelse( EncProb_tc[1+cbind(TmbData$t_i,TmbData$c_i)]==0, NA, TmbData$b_i )
  # Build object
  TmbList = Build_TMB_Fn("build_model"=TRUE, "TmbData"=TmbData, "RunDir"=test_path, "Version"=Version_VAST, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Npool"=Npool)
  # Optimize
  on.exit( dyn.unload(paste0(test_path,"/",TMB::dynlib(Version_VAST))), add=TRUE )
  Opt = TMBhelper::Optimize( obj=TmbList[["Obj"]], getsd=TRUE, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], newtonsteps=1 )  # , rel.tol=1e-20
  # Comparisons
  Par1 = Opt$par[names(Opt$par)%in%c("ln_H_input","beta1_ct","logkappa1","beta2_ct","logkappa1","logSigmaM")]
  Par2 = opt$par[names(opt$par)%in%c("ln_H_input","beta1_ct","logkappa1","beta2_ct","logkappa1","logSigmaM")]
  expect_equal( as.vector(Par1), as.vector(Par2), tolerance=1e-3 )
  #
  #Index = plot_biomass_index( TmbData=TmbData, Sdreport=Opt$SD, DirName=test_path )
  #calculate_proportion( TmbData=TmbData, Index=Index, DirName=test_path, Year_Set=unique(Data_Geostat$Year) )
})

