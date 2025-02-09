
#' Convert from Lat-Long to UTM
#'
#' \code{project_coordinates} converts from Latitude-Longitude to Universal Transverse Mercator projections for a given location
#'
#' @param projargs A character string of projection arguments; the arguments must be entered exactly as in the PROJ.4 documentation; if the projection is unknown, use \code{as.character(NA)}, it may be missing or an empty string of zero length and will then set to the missing value. With \pkg{rgdal} built with PROJ >= 6 and GDAL >= 3, the \code{+init=} key may only be used with value \code{epsg:<code>}. From \pkg{sp} version 1.4-4, the string associated with the SRS_string argument may be entered as-is and will be set as SRS_string if the projargs argument does not begin with a \code{+} (suggested by Mikko Vihtakari).
#' @param X vector of horizontal-coordinates (e.g., longitude by default)
#' @param Y vector of vertical-coordinates (e.g., latitude by default)
#' @param zone DEPRECATED INPUT; UTM zone (integer between 1 and 60) or alphanumeric CRS code used by package rgdal to convert latitude-longitude coordinates to projection in kilometers; \code{zone=NA} uses UTM and automatically detects the appropriate zone
#' @param flip_around_dateline DEPRECATED INPUT; boolean specifying whether to flip Lat-Lon locations around the dateline, and then retransform back (only useful if Lat-Lon straddle the dateline)

#' @return A data frame with the following columns
#' \describe{
#'   \item{X}{The horizonal coordinates after projection (e.g., Eastings in UTM kilometers by default)}
#'   \item{Y}{The vertical coordinates after projection (e.g., Northings in UTM kilometers by default)}
#' }

#' @export
project_coordinates <-
function( X,
          Y,
          projargs = NA,
          origargs = "+proj=longlat +datum=WGS84",
          zone = NA,
          flip_around_dateline = FALSE,
          silent = FALSE ){

  # Original projection
  origCRS = sp::CRS(origargs)

  # Default -- Use zone value if supplied for backwards compatibility
  if( !is.na(zone) ){
    if( !is.na(projargs) ){
      stop("Please do not specify both `zone` and `projargs`")
    }
    # Transform around dateline if requested
    if( flip_around_dateline==TRUE ){
      zone = zone + 30
      X = X + 180
    }
    # Project to UTM
    xy_sp = sp::SpatialPoints( coords=cbind(X,Y), proj4string=origCRS ) # expects in long-lat format
    projargs = paste0("+proj=utm +datum=WGS84 +units=km",ifelse(mean(Y)<0," +south","")," +zone=",zone)
    if(silent==FALSE) message("For the UTM conversion, used zone ",zone," as specified")
  }

  # New option
  if( is.na(zone) ){
    # Auto-detect the UTM zone to match PBSmapping
    if( is.na(projargs) && length(grep("+proj=longlat",origargs))==1 ){
      lon1 = X
      lon2 = ifelse( lon1 < 0, lon1 + 360, lon1 )
      # Determine whether to transform around dateline or not
      if( diff(range(lon1)) < diff(range(lon2)) ){
        mean_lon = mean(lon1)
      }else{
        mean_lon = mean(lon2)
      }
      # Detect UTM based on mean of longitudes to match PBSmapping behavior
      zone = 1 + floor( (mean_lon-180)/6 )
      zone = (zone-1) %% 60 + 1
      projargs = paste0("+proj=utm +datum=WGS84 +units=km",ifelse(mean(Y)<0," +south","")," +zone=",zone)
      if(silent==FALSE) message("For the UTM conversion, automatically detected zone ",zone,".")
    }
    xy_sp = sp::SpatialPoints( coords=cbind(X,Y), proj4string=origCRS ) # expects in long-lat format
  }

  # Check for issues
  if( !is.na(zone) && (zone<1 | zone>60) ) stop("Check zone in `project_coordinates(.)`")
  projCRS = sp::CRS( projargs )

  # Conduct projection
  xynew = sp::spTransform( xy_sp, projCRS )
  xynew_i = xynew@coords
  colnames(xynew_i) = c("X", "Y")

  # Add attributes such that it includes everything in `PBSmapping::convUL`, i.e., attribute "zone"
  attr(xynew_i,"zone") = zone
  attr(xynew_i,"origargs") = origargs
  attr(xynew_i,"projargs") = projargs
  attr(xynew_i,"origCRS") = origCRS
  attr(xynew_i,"projCRS") = projCRS

  # Return results
  return( xynew_i )
}
