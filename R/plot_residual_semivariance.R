#' @title
#' Plot semivariance for quantile residuals
#'
#' @description
#' \code{plot_residual_semivariance} calculates the spatial and temporal semivariance from standard-normal quantile residuals
#'
#' This function plots takes quantile residuals, converts them to a standard-normal distribution,
#'      fits the two-dimensional semi-variance in space and time, and plots this.
#'
#' @export
plot_residual_semivariance <-
function( fit,
          dharma_raster,
          dharmaRes,
          file_name = "quantile_residuals_semivariance",
          working_dir = getwd() ){

  # Transform
  tmp = dharma_raster$Raster_proj
  for(i in seq_along(tmp)) raster::values(tmp[[i]]) = qnorm((1-1/dharmaRes$nSim)*raster::values(tmp[[i]]) + 1/dharmaRes$nSim/2)
  tmp = raster::stack(tmp)

  # Set up
  sp = raster::rasterToPoints(tmp, spatial=TRUE)
  time = as.POSIXct( paste0(fit$year_labels[fit$years_to_plot],"-01-01") )
  mydata = data.frame(values = unlist(sp@data) ) # , ID=IDs)
  if(length(time)==1){
    endTime = time+1
  }else{
    endTime = spacetime::delta(time)
  }
  stfdf = spacetime::STFDF(sp=sp, time=time, data=mydata, endTime=endTime)
  residual_semivariance = gstat::variogram(values~1, stfdf, width=20, cutoff = 500, tlags=0:min(7,length(time)-1) )

  # Plot
  png( file=file.path(working_dir,paste0(file_name,".png")), width=6, height=5, units="in", res=200 )
    # gstat:::plot.gstatVariogram
    # gstat:::plot.StVariogram
    semivariance_plot = plot( residual_semivariance,
          at = seq(0, max(residual_semivariance$gamma,na.rm=TRUE), length=23) )#, wireframe=TRUE), plot.numbers=TRUE
    plot(semivariance_plot)
  dev.off()

  #
  return( invisible(residual_semivariance) )
}
