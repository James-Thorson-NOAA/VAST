
#' @export
Prepare_AI_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, zone=NA, flip_around_dateline=TRUE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( aleutian_islands_grid, package="VAST" )
  Data_Extrap <- aleutian_islands_grid

  # Survey areas
  Area_km2_x = Data_Extrap[,'Area_km2']

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LAT_DD"=Data_Extrap[,'Lat'], "propInSurvey"=1)
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }
             #
  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_BC_Coast_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, strata_to_use=c('SOG','WCVI','QCS','HS','WCHG'), zone=NA, flip_around_dateline=FALSE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( bc_coast_grid, package="VAST" )
  bc_coast_grid <- cbind( bc_coast_grid, "ALL"=13.74)
  Data_Extrap <- bc_coast_grid

  # Survey areas
  Area_km2_x = rowSums( Data_Extrap[,strata_to_use,drop=FALSE] )

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LAT_DD"=Data_Extrap[,'Lat'], "BEST_LON_DD"=Data_Extrap[,'Lon'])
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=length(strata.limits), dimnames=list(NULL,names(strata.limits))))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }

  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_BSslope_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, zone=NA, flip_around_dateline=FALSE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( bering_sea_slope_grid, package="VAST" )
  Data_Extrap <- bering_sea_slope_grid

  # Survey areas
  Area_km2_x = Data_Extrap[,'Area_in_survey_km2']

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=Data_Extrap[,'Depth_km']*1000, "BEST_LAT_DD"=Data_Extrap[,'Lat'], "propInSurvey"=1)
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }
             #
  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_EBS_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, zone=NA, flip_around_dateline=TRUE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( eastern_bering_sea_grid, package="VAST" )
  Data_Extrap <- eastern_bering_sea_grid

  # Survey areas
  #Area_km2_x = 4 * 1.852^2 * ifelse( Data_Extrap[,'EBS_STRATUM']!=0, 1, 0 )
  Area_km2_x = Data_Extrap[,'Area_in_survey_km2']

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LON_DD"=Data_Extrap[,'Lon'], "BEST_LAT_DD"=Data_Extrap[,'Lat'], "propInSurvey"=1)
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }
             #
  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap = cbind( Data_Extrap, "E_km"=tmpUTM[,'X'], "N_km"=tmpUTM[,'Y'] )

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_GOA_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, zone=NA, flip_around_dateline=FALSE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( gulf_of_alaska_grid, package="VAST" )
  Data_Extrap <- gulf_of_alaska_grid

  # Survey areas
  #Area_Tmp = 4 * 1.852^2 * ifelse( Data_Extrap[,'EBS_STRATUM']!=0, 1, 0 )
  Area_km2_x = Data_Extrap[,'Area_in_survey_km2']

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LAT_DD"=Data_Extrap[,'Lat'], "BEST_LON_DD"=Data_Extrap[,'Lon'], "propInSurvey"=1)
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }

  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_GOM_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, zone=NA, ... ){

  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( gulf_of_mexico_grid, package="VAST" )
  Data_Extrap <- gulf_of_mexico_grid

  # Survey areas
  Area_km2_x = Data_Extrap[,'Area_in_survey_km2']

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=Data_Extrap[,'Depth_m'], "BEST_LAT_DD"=Data_Extrap[,'Lat'], "BEST_LON_DD"=Data_Extrap[,'Lon'], "propInSurvey"=1)
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }

  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=FALSE)                                                         #$
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=FALSE, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_GSL_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, zone=NA, flip_around_dateline=FALSE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( gulf_of_st_lawrence_grid, package="VAST" )
  Data_Extrap <- gulf_of_st_lawrence_grid

  # Survey areas
  Area_km2_x = Data_Extrap[,'area'] # Convert from nm^2 to km^2

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LAT_DD"=Data_Extrap[,'latitude'], "BEST_LON_DD"=Data_Extrap[,'longitude'])
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=length(strata.limits), dimnames=list(NULL,names(strata.limits))))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }

  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'longitude'], Lat=Data_Extrap[,'latitude'], zone=zone, flip_around_dateline=flip_around_dateline)
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]
  colnames(Data_Extrap)[1:2] = c("Lat","Lon")

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_HabCam_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, zone=NA, flip_around_dateline=TRUE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( habcam_grid, package="VAST" )
  Data_Extrap <- habcam_grid

  # Survey areas
  Area_km2_x = Data_Extrap[,'Cell_area_km2']

  # Augment with strata for each extrapolation cell
  Tmp = cbind("Longitude"=Data_Extrap[,'Lon'], "BEST_LAT_DD"=Data_Extrap[,'Lat'], "propInSurvey"=1)
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }
             #
  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_NBS_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, zone=NA, flip_around_dateline=FALSE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( northern_bering_sea_grid, package="VAST" )
  Data_Extrap <- northern_bering_sea_grid

  # Survey areas
  Area_km2_x = Data_Extrap[,'Area_in_survey_km2']

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LAT_DD"=Data_Extrap[,'Lat'], "BEST_LON_DD"=Data_Extrap[,'Lon'], "propInSurvey"=1)
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }
             #
  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap = cbind( Data_Extrap, "E_km"=tmpUTM[,'X'], "N_km"=tmpUTM[,'Y'] )

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_Chukchi_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, zone=NA, flip_around_dateline=FALSE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( chukchi_sea_grid, package="VAST" )
  Data_Extrap <- chukchi_sea_grid

  # Survey areas
  Area_km2_x = Data_Extrap[,'Area_in_survey_km2']

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LAT_DD"=Data_Extrap[,'Lat'], "BEST_LON_DD"=Data_Extrap[,'Lon'], "propInSurvey"=1)
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }
             #
  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap = cbind( Data_Extrap, "E_km"=tmpUTM[,'X'], "N_km"=tmpUTM[,'Y'] )

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_NWA_Extrapolation_Data_Fn <-
function( strata.limits=NULL, epu_to_use = c('All', 'Georges_Bank','Mid_Atlantic_Bight','Scotian_Shelf','Gulf_of_Maine','Other')[1], projargs=NA, zone=NA, flip_around_dateline=FALSE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = list('All_areas'=1:1e5)
  }
  message("Using strata ", strata.limits)

  if(any(tolower(epu_to_use) %in% "all")) {
    epu_to_use <- c('Georges_Bank','Mid_Atlantic_Bight','Scotian_Shelf','Gulf_of_Maine','Other')
  }

  # Read extrapolation data
  utils::data( northwest_atlantic_grid, package="VAST" )
  Data_Extrap <- northwest_atlantic_grid

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LAT_DD"=Data_Extrap[,'Lat'], "BEST_LON_DD"=Data_Extrap[,'Lon'])
  if( length(strata.limits)==1 && strata.limits[1]=="EPU" ){
    # Specify epu by 'epu_to_use'
    Data_Extrap <- Data_Extrap[Data_Extrap$EPU %in% epu_to_use, ]
    Data_Extrap$EPU <- droplevels(Data_Extrap$EPU)

    a_el = matrix(NA, nrow=nrow(Data_Extrap), ncol=length(epu_to_use), dimnames=list(NULL, epu_to_use) )
    Area_km2_x = Data_Extrap[, "Area_in_survey_km2"]
    for(l in 1:ncol(a_el)){
      a_el[,l] = ifelse( Data_Extrap[,'EPU'] %in% epu_to_use[[l]], Area_km2_x, 0 )
    }
  }else{
    # Specify strata by 'stratum_number'
    a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=length(strata.limits), dimnames=list(NULL,names(strata.limits))))

        # Survey areas
    Area_km2_x = Data_Extrap[,'Area_in_survey_km2']
    for(l in 1:ncol(a_el)){
      a_el[,l] = ifelse( Data_Extrap[,'stratum_number'] %in% strata.limits[[l]], Area_km2_x, 0 )
    }
  }

  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_NZ_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, zone=NA, survey="Chatham_rise", flip_around_dateline=FALSE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( chatham_rise_grid, package="VAST" )
  Data_Extrap <- chatham_rise_grid

  # Survey areas
  #Area_km2_x = 4 * 1.852^2 * ifelse( Data_Extrap[,'EBS_STRATUM']!=0, 1, 0 )
  Area_km2_x = Data_Extrap[,'area_km2']

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LAT_DD"=Data_Extrap[,'Lat'], "propInSurvey"=1)
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }

  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #if( is.na(zone) || is.numeric(zone) ){
  #  tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)
  #  colnames(tmpUTM)[3:4] = c('E_km','N_km')
  #}else{
  #  tmpUTM = Convert_LL_to_EastNorth_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], crs=zone )
  #}
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_Other_Extrapolation_Data_Fn <-
function( strata.limits, observations_LL, projargs=NA, grid_dim_km=c(2,2), maximum_distance_from_sample=NULL,
          zone=NA, grid_in_UTM=TRUE, grid_dim_LL=c(0.1,0.1), flip_around_dateline=FALSE, ... ){

  # Local function
  rename_columns = function( DF, origname=colnames(DF), newname ){
    DF_new = DF
    for(i in 1:length(origname)){
      Match = match( origname[i], colnames(DF_new) )
      if(length(Match)==1) colnames(DF_new)[Match] = newname[i]
    }
    return(DF_new)
  }

  if( grid_in_UTM==TRUE ){
    # Fill in missing inputs
    if( is.null(maximum_distance_from_sample) ) maximum_distance_from_sample = sqrt((grid_dim_km[1]/2)^2+(grid_dim_km[2]/2)^2)

    # Get range
    #TmpUTM = Convert_LL_to_UTM_Fn( Lon=observations_LL[,'Lon'], Lat=observations_LL[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$
    TmpUTM = project_coordinates( X=observations_LL[,'Lon'], Y=observations_LL[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$
    E_lim = mean(range(TmpUTM[,'X'])) + c(-0.6,0.6)*diff(range(TmpUTM[,'X']))
    N_lim = mean(range(TmpUTM[,'Y'])) + c(-0.6,0.6)*diff(range(TmpUTM[,'Y']))

    # Detect Northern or Southern hemisphere
    #NorthernTF = all( observations_LL[,'Lat']>0 )
    #if( any(observations_LL[,'Lat']>0) & any(observations_LL[,'Lat']<0) ) warning( "PBSmapping doesn't work with observations in both Northern and Southern hemisphere" )

    # Make grid
    E_seq = seq(E_lim[1], E_lim[2], by=grid_dim_km[1])
    N_seq = seq(N_lim[1], N_lim[2], by=grid_dim_km[2])
    if( length(E_seq)*length(N_seq) > 1e5 ){
      warning("Please increase `grid_dim_km` to avoid memory issues caused by creating a large extrapolation-grid")
    }
    if( length(E_seq)*length(N_seq) > 1e7 ){
      stop("Please increase `grid_dim_km` to avoid memory-associated crashes caused by creating a huge extrapolation-grid")
    }
    Data_Extrap = expand.grid( "E_km"=E_seq, "N_km"=N_seq, "Area_km2"=prod(grid_dim_km) )

    # Add LL
    # NOTE:  inputs to project_coordinates are reversed due to projecting from projargs back to origargs
    #TmpUTM = rename_columns( Data_Extrap[,c("E_km","N_km")], newname=c("X","Y"))
    #attr(TmpUTM, "projection") = "UTM"
    #attr(TmpUTM, "zone") = attr(TmpUTM,"zone")
    #TmpLL = PBSmapping::convUL(TmpUTM, southern=!NorthernTF )
    #if( flip_around_dateline==TRUE ){
    #  TmpLL[,1] = TmpLL[,1] - 180
    #}
    TmpLL = project_coordinates( X=Data_Extrap[,"E_km"], Y=Data_Extrap[,"N_km"], projargs=attr(TmpUTM,"origargs"), origargs=attr(TmpUTM,"projargs"))
    Data_Extrap = cbind( Data_Extrap, rename_columns(TmpLL,newname=c("Lon","Lat")) )

    # Restrict to grid locations near samples
    NN_Extrap = RANN::nn2( query=Data_Extrap[,c("E_km","N_km")], data=TmpUTM[,c("X","Y")], k=1)
    Data_Extrap = cbind( Data_Extrap, "Include"=ifelse(NN_Extrap$nn.dists<maximum_distance_from_sample,1,0))

    # Survey areas
    Area_km2_x = Data_Extrap[,'Area_km2'] * Data_Extrap[,'Include']
  }else{
    # Fill in missing inputs
    if( is.null(maximum_distance_from_sample) ) maximum_distance_from_sample = sqrt((grid_dim_LL[1]/2)^2+(grid_dim_LL[2]/2)^2)

    # Get range
    Lat_lim = mean(range(observations_LL[,'Lat'])) + c(-0.6,0.6)*diff(range(observations_LL[,'Lat']))
    if( flip_around_dateline==FALSE ) Lon_lim = mean(range(observations_LL[,'Lon'])) + c(-0.6,0.6)*diff(range(observations_LL[,'Lon']))
    if( flip_around_dateline==TRUE ){
      Tmp_Lon = 180 + ifelse( observations_LL[,'Lon']>0, observations_LL[,'Lon']-360, observations_LL[,'Lon'])
      Lon_lim = mean(range(Tmp_Lon)) + c(-0.6,0.6)*diff(range(Tmp_Lon))
    }

    # Make grid
    Data_Extrap = expand.grid( "Lon"=seq(Lon_lim[1],Lon_lim[2],by=grid_dim_LL[1]), "Lat"=seq(Lat_lim[1],Lat_lim[2],by=grid_dim_LL[2]), "Area_km2"=110^2*prod(grid_dim_LL) )
    if( flip_around_dateline==TRUE ){
      Data_Extrap[,'Lon'] = ifelse( Data_Extrap[,'Lon']<0, Data_Extrap[,'Lon']+360, Data_Extrap[,'Lon'] ) - 180
    }

    # Add empty UTM
    Data_Extrap = cbind( Data_Extrap, "E_km"=NA, "N_km"=NA, "Area_km2"=prod(grid_dim_km) )
    TmpUTM = NA
    attr(TmpUTM, "zone") = NA

    # Restrict to grid locations near samples
    NN_Extrap = RANN::nn2( query=Data_Extrap[,c("Lat","Lon")], data=observations_LL[,c('Lat','Lon')], k=1)
    Data_Extrap = cbind( Data_Extrap, "Include"=ifelse(NN_Extrap$nn.dists<maximum_distance_from_sample,1,0))

    # Survey areas
    Area_km2_x = Data_Extrap[,'Area_km2'] * Data_Extrap[,'Include']
  }

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LAT_DD"=Data_Extrap[,'Lat'], "BEST_LON_DD"=Data_Extrap[,'Lon'])
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(TmpUTM,"zone"), "projargs"=attr(TmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_SA_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, region=c("south_coast","west_coast"), zone=NA, flip_around_dateline=FALSE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( south_africa_grid, package="VAST" )
  Data_Extrap <- south_africa_grid

  # Survey areas
  Area_km2_x = Data_Extrap[,'nm2'] * 1.852^2 * ifelse( Data_Extrap[,'stratum']%in%region, 1, 0 ) # Convert from nm^2 to km^2

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LAT_DD"=Data_Extrap[,'cen_lat'], "BEST_LON_DD"=Data_Extrap[,'cen_long'])
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=length(strata.limits), dimnames=list(NULL,names(strata.limits))))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }

  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'cen_long'], Lat=Data_Extrap[,'cen_lat'], zone=zone, flip_around_dateline=flip_around_dateline)
  tmpUTM = project_coordinates( X=Data_Extrap[,'cen_long'], Y=Data_Extrap[,'cen_lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]
  colnames(Data_Extrap)[1:2] = c("Lon","Lat")

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_SMI_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, zone=NA, flip_around_dateline=TRUE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( st_matthews_island_grid, package="VAST" )
  Data_Extrap <- st_matthews_island_grid

  # Survey areas
  #Area_km2_x = 4 * 1.852^2 * ifelse( Data_Extrap[,'EBS_STRATUM']!=0, 1, 0 )
  Area_km2_x = Data_Extrap[,'Area_in_survey_km2']

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LAT_DD"=Data_Extrap[,'Lat'], "propInSurvey"=1)
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }
             #
  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_User_Extrapolation_Data_Fn <-
function( input_grid, strata.limits=NULL, projargs=NA, zone=NA, flip_around_dateline=TRUE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  Data_Extrap <- input_grid

  # Survey areas
  Area_km2_x = Data_Extrap[,'Area_km2']

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_LAT_DD"=Data_Extrap[,'Lat'], "BEST_LON_DD"=Data_Extrap[,'Lon'])
  if( "Depth" %in% colnames(Data_Extrap) ){
    Tmp = cbind( Tmp, "BEST_DEPTH_M" = Data_Extrap[,'Depth'] )
  }
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }

  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  if( all(c("E_km","N_km") %in% colnames(Data_Extrap)) ){
    Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]
  }else{
    Data_Extrap = cbind( Data_Extrap, 'E_km'=tmpUTM[,'X'], 'N_km'=tmpUTM[,'Y'] )
  }

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_WCGBTS_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, surveyname='propInWCGBTS', zone=NA, flip_around_dateline=FALSE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( california_current_grid, package="VAST" )
  Data_Extrap <- california_current_grid

  # Survey areas
  Area_km2_x = 4*apply(Data_Extrap[,surveyname,drop=FALSE], MARGIN=1, FUN=min)

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=(-1000)*Data_Extrap[,'Depth_km'], "BEST_LAT_DD"=Data_Extrap[,'Lat'])
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    #a_el[,l] = apply(Tmp, MARGIN=1, FUN=nwfscDeltaGLM::strata.fn, Strata.df=strata.limits[l,,drop=FALSE])
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }

  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$
  tmpUTM = project_coordinates( X=Data_Extrap[,'Lon'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)
  Data_Extrap = cbind( Data_Extrap, 'Include'=(Data_Extrap[,'Cowcod']==0 & Data_Extrap[,'Ngdc_m']<(-35)))
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_WCGHL_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, zone=NA, flip_around_dateline=FALSE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( WCGHL_grid, package="VAST" )
  Data_Extrap <- WCGHL_grid

  # Survey areas
  Area_km2_x = Data_Extrap[,'Area_km2'] # Convert from nm^2 to km^2

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=0, "BEST_LAT_DD"=Data_Extrap[,'Lat'], "BEST_LON_DD"=Data_Extrap[,'Long'])
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=length(strata.limits), dimnames=list(NULL,names(strata.limits))))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }

  # Convert extrapolation-data to an Eastings-Northings coordinate system
  #tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Long'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)
  tmpUTM = project_coordinates( X=Data_Extrap[,'Long'], Y=Data_Extrap[,'Lat'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]
  colnames(Data_Extrap)[1:2] = c("Lat","Lon")

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

#' @export
Prepare_BFISH_MHI_Extrapolation_Data_Fn <-
function( strata.limits=NULL, projargs=NA, zone=NA, flip_around_dateline=FALSE, ... ){
  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  utils::data( BFISH_MHI, package="VAST" )
  Data_Extrap <- BFISH_MHI

  # Survey areas
  Area_km2_x = Data_Extrap[,'Area_km2']

  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=Data_Extrap[,'Depth_MEAN_m'], "BEST_LAT_DD"=Data_Extrap[,'lat_deg'], "BEST_LON_DD"=Data_Extrap[,'lon_deg'])
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=length(strata.limits), dimnames=list(NULL,names(strata.limits))))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }

  # Convert extrapolation-data to an Eastings-Northings coordinate system
  tmpUTM = project_coordinates( X=Data_Extrap[,'lon_deg'], Y=Data_Extrap[,'lat_deg'], projargs=projargs, zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]
  colnames(Data_Extrap)[1:2] = c("Lon","Lat")

  # Return
  Return = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "projargs"=attr(tmpUTM,"projargs"),
    "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  return( Return )
}

