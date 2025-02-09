#' @export
match_strata_fn <- function(x, strata_dataframe) {
  # Default all strata
  match_latitude_TF = match_longitude_TF = match_depth_TF = rep( TRUE, nrow(strata_dataframe))
  if( all(c("south_border","north_border") %in% names(strata_dataframe)) ){
    match_latitude_TF = as.numeric(x["BEST_LAT_DD"])>strata_dataframe[,'south_border'] & as.numeric(x["BEST_LAT_DD"])<=strata_dataframe[,'north_border']
  }
  if( all(c("west_border","east_border") %in% names(strata_dataframe)) ){
    match_longitude_TF = as.numeric(x["BEST_LON_DD"])>strata_dataframe[,'west_border'] & as.numeric(x["BEST_LON_DD"])<=strata_dataframe[,'east_border']
  }
  if( all(c("shallow_border","deep_border") %in% names(strata_dataframe)) ){
    match_depth_TF = as.numeric(x["BEST_DEPTH_M"])>strata_dataframe[,'shallow_border'] & as.numeric(x["BEST_DEPTH_M"])<=strata_dataframe[,'deep_border']
  }
  # Return stuff
  Char = as.character(strata_dataframe[match_latitude_TF & match_longitude_TF & match_depth_TF,"STRATA"]) 
  return(ifelse(length(Char)==0,NA,Char))
}
