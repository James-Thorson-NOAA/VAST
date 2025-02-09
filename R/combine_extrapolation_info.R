
#' Combines multiple extrapolation grids
#'
#' \code{combine_extrapolation_info} combines multiple extrapolation grids when combining data from multiple surveysn
#'
#' @param ... a sequence of outputs from \code{make_extrapolation_info}

#' @return Identical output from \code{make_extrapolation_info}, but combined from each input

#' @export
combine_extrapolation_info = function( ..., create_strata_per_region=FALSE ){

  input_list = list( ... )

  for( lI in 1:length(input_list)){
    if( !all( c("a_el","Data_Extrap","zone","flip_around_dateline","Area_km2_x") %in% names(input_list[[lI]]) )){
      stop( "Check inputs to `combine_extrapolation_grids`" )
    }
  }
  Zone = sapply( input_list, FUN=function(List){List[["zone"]]} )
  Flip = sapply( input_list, FUN=function(List){List[["flip_around_dateline"]]} )
  Projargs = sapply( input_list, FUN=function(List){List[["projargs"]]} )
  if( all(!is.na(Projargs)) ){
    if( length(unique(Projargs)) > 1 ){
      stop( "Must use same `projargs` for projection of all extrapolation grids" )
    }
  }else{
    if( sd(Zone)>0 | sd(Flip)>0 ){
      stop( "Must use same Zone for UTM conversion for all extrapolation grids" )
    }
  }

  # Combine stuff
  Data_Extrap = Area_km2_x = NULL
  a1_el = matrix(0, nrow=0, ncol=1)
  a2_el = matrix(0, nrow=0, ncol=0)
  #assign( x="input_list", value=input_list, envir = .GlobalEnv )

  for( lI in 1:length(input_list) ){
    #Tmp = input_list[[lI]]$Data_Extrap
    #colnames(Tmp) = ifelse( colnames(Tmp)=="Area_in_survey_km2", "Area_km2", colnames(Tmp) )
    #Data_Extrap = rbind( Data_Extrap, Tmp[,c('E_km','N_km','Lon','Lat','Include','Area_km2')] )

    # Warnings
    if( ncol(input_list[[lI]]$a_el)>1 ){
      if( !(create_strata_per_region==TRUE & lI==1) ){
        stop("`combine_extrapolation_info` isn't designed to combine regions with multiple identified strata, except when `create_strata_per_region=TRUE`")
      }
    }

    # Combine Data_Extrap
    Data_Extrap = rbind( Data_Extrap, input_list[[lI]]$Data_Extrap[,c('E_km','N_km','Lon','Lat','Include')] )

    # Combine area vector
    Area_km2_x = c( Area_km2_x, input_list[[lI]]$Area_km2_x )

    # Combine strata definitions
    a1_el = rbind( as.matrix(a1_el), as.matrix(input_list[[lI]]$a_el[,1,drop=FALSE]) )

    # Make one stratum per region
    a2_el = rbind(
      as.matrix(cbind(a2_el, matrix(0,nrow=nrow(a2_el),ncol=ncol(input_list[[lI]]$a_el)))),
      as.matrix(cbind(matrix(0,nrow=nrow(input_list[[lI]]$a_el),ncol=ncol(a2_el)), input_list[[lI]]$a_el))
    )
  }

  # Only pass back the
  if( create_strata_per_region==TRUE ){
    a_el = a2_el
  }else{
    a_el = a1_el
  }

  # Return
  Return = list( "a_el"=a_el,
                 "Data_Extrap"=Data_Extrap,
                 "zone"=Zone[1],
                 "flip_around_dateline"=Flip[1],
                 "projargs"=Projargs[1],
                 "Area_km2_x"=Area_km2_x)

  # Fix rowname duplicates
  rownames(Return$Data_Extrap) = seq_len(nrow(Return$Data_Extrap))

  # Add units
  units(Return$a_el) = "km^2"
  units(Return$Area_km2_x) = "km^2"

  # Return
  class(Return) = "make_extrapolation_info"
  return( Return )
}

