

#' Convert shapefile to extrapolation-grid
#'
#' \code{convert_shapefile} reads in a shapefile and creates an extrapolation with a chosen resolution
#'
#' @inheritParams sp::CRS
#' @inheritParams make_extrapolation_info
#'
#' @param file_path path for shapefile on harddrive
#' @param make_plots Boolean indicating whether to visualize inputs and outputs as maps
#' @param projargs_for_shapefile projection-arguments (e.g., as parsed by \code{sp::CRS}), that are used when reading shapefile and overriding any projection arguments already saved there;  Default \code{projargs_for_shapefile=NULL} uses the projection arguments available in the shapefile
#' @param projargs A character string of projection arguments used to project the shapefile prior to construction an extrapolation-grid; the arguments must be entered exactly as in the PROJ.4 documentation; Default \code{projargs=NULL} uses UTM and infers the UTM zone based on the centroid of the shapefile.

#' @return extrapolation-grid

#' @examples
#' \dontrun{
#'  convert_shapefile( file_path="C:/Users/James.Thorson/Desktop/Work files/AFSC/2020-03 -- Add ICES grids/IBTS grids/BITS/Shapefile.shp", make_plots=TRUE )
#' }
#'
#' @author Cecilia O'Leary, James Thorson
#' @export
convert_shapefile = function( file_path,
    projargs = NULL,
    grid_dim_km = c(2,2),
    projargs_for_shapefile = NULL,
    make_plots = FALSE,
    quiet = TRUE,
    area_tolerance = 0.05,
    ... ){

  #shapefile_input = rgdal::readOGR( file_path, verbose=FALSE, p4s=projargs_for_shapefile )
  shapefile_input = sf::st_read( file_path )

  # Need to run sf::st_make_valid to ensure sf::st_area works later
  shapefile_input = sf::st_make_valid( shapefile_input )
  #  sf::st_crs(shapefile_input) = sf::st_crs(projargs_for_shapefile)
  # rgdal::writeOGR
  message("Reading shapefile with projargs: ", print(sf::st_crs(shapefile_input)))
  #
  if( is.null(sf::st_crs(shapefile_input)) ){
    sf::st_crs(shapefile_input) = sf::st_crs(projargs_for_shapefile)
  }
  # raster::shapefile(.) has simplified read-write interface for future reference
  #if( !(shapefile_input %in% c("SpatialPolygonsDataFrame","SpatialPolygons")) ){
  #  warning( "object at `file_path` doesn't appear to be a shapefile")
  #}
  proj_orig = "+proj=longlat +ellps=WGS84 +no_defs"
  shapefile_orig = sf::st_transform(shapefile_input, crs=sf::st_crs(proj_orig) )
  #shapefile_orig@proj4string = sp::CRS(proj_orig)

  # Infer projargs if missing, and project
  utm_zone = floor((mean(sf::st_bbox(shapefile_orig)[c('xmin','xmax')]) + 180) / 6) + 1
  projargs_utm = paste0("+proj=utm +zone=",utm_zone," +ellps=WGS84 +datum=WGS84 +units=km +no_defs ")
  if( is.null(projargs) || is.na(projargs) ){
    projargs = projargs_utm
  }
  shapefile_proj = sf::st_transform(shapefile_orig, crs=sf::st_crs(projargs) )

  # Check for total area
  Total_area = as.numeric(sum(sf::st_area(shapefile_proj)))
  if( (Total_area) < 100 ){
    warning("Total area for polygons in shapefile is ", Total_area, " km^2 and this might indicate an issue with identifying or specifying the projection" )
  }

  # Determine bounds for box
  bbox = sf::st_bbox(shapefile_proj)
  bbox[c('xmin','xmax')] = floor(bbox[c('xmin','xmax')] / grid_dim_km[1]) * grid_dim_km[1]
  bbox[c('ymin','ymax')] = ceiling(bbox[c('ymin','ymax')] / grid_dim_km[2]) * grid_dim_km[2]

  # Make grid
  grid_proj = expand.grid(X=seq(bbox['xmin'],bbox['xmax'],by=grid_dim_km[1]), Y=seq(bbox['ymin'],bbox['ymax'],by=grid_dim_km[2]))
  grid_proj = sf::st_as_sf( grid_proj[,c("X","Y")], coords=1:2, crs=sf::st_crs(projargs) )

  # Restrict points to those within a polygon
  grid_proj = sf::st_intersection(grid_proj, shapefile_proj)
  grid_proj_coordinates = sf::st_coordinates(grid_proj)

  # Return to original coordinates
  grid_orig = sf::st_transform( grid_proj, crs=sf::st_crs(proj_orig) )
  grid_orig = sf::st_coordinates(grid_orig)

  # Combine
  extrapolation_grid = data.frame( "Lat"=grid_orig[,'Y'], "Lon"=grid_orig[,'X'], "Area_km2"=prod(grid_dim_km), # "Stratum"=grid_orig[,'AreaName'],
    "Include"=1, "E_km"=grid_proj_coordinates[,'X'], "N_km"=grid_proj_coordinates[,'Y'] )

  # make plots
  if( make_plots==TRUE ){
    par( mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )

    # Plot projected shapefile
    sp::plot( shapefile_proj, main="shapefile in original projection" )
    axis(1); axis(2); box()

    # Plot projected grid
    sp::plot( grid_proj, main="Grid in new projection" )
    axis(1); axis(2); box()

    # plot original shapefile
    sp::plot( shapefile_orig, main="shapefile in new projection" )
    axis(1); axis(2); box()

    # plot original shapefile
    sp::plot( grid_orig[,c('X','Y')], main="Grid in original coordinates" )
    axis(1); axis(2); box()
  }

  # Compare areas
  area_shapefile_input = as.numeric(sum( sf::st_area(shapefile_input) ))
  area_shapefile_orig = as.numeric(sum( sf::st_area(shapefile_orig) )) / 1000^2  # Convert to square-kiometers
  area_shapefile_proj = as.numeric(sum( sf::st_area(shapefile_proj) ))
  area_grid_proj = sum( extrapolation_grid[,'Area_km2'] )

  # Messages
  Area_ratio = area_grid_proj / area_shapefile_orig
  if( quiet==FALSE ){
    message( "Area of projected extrapolation-grid is ", formatC(Area_ratio*100,format="f",digits=2), "% of the original shapefile area" )
  }
  if( Area_ratio>(1+area_tolerance) | Area_ratio<(1-area_tolerance) ){
    stop( "Please use a different projection to decrease this conversion error; perhaps `projargs=", projargs_utm, "`" )
  }

  # Return output
  Return = list( "extrapolation_grid"=extrapolation_grid, "projargs"=projargs )
  return( Return )
}

