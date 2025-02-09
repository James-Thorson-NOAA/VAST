
#' Plot timeseries estimates
#'
#' \code{plot_timeseries} plots a timeseries with standard errors shown as whiskers or a shaded polygon
#'
#' @param x numeric-vector of x-coordinates to plot
#' @param y numeric-vector of y-coordinates to plot
#' @param ysd numeric-vector of standard errors in y-coordinates, used in normal approximation to interval
#' @param ybounds alternate specification of y-coordinate whiskers; useful to avoid normal approximation
#' @param fn what plotting function to use; default \code{fn=lines} does not create a new plotting window
#' @param bounds_type which type either \code{"whiskers"} or \code{"shading"}
#' @param bounds_args tagged list of additional arguments to pass to \code{\link[graphics]{lines}} or \code{\link[graphics]{polygon}}
#' @param interval_width width of interval in normal approximation; only used when specifying \code{ysd} without \code{ybounds}
#'
#' @export
plot_timeseries <-
function( x,
          y,
          y_sd,
          ybounds,
          fn = lines,
          ylim  =  NULL,
          bounds_type = "whiskers",
          bounds_args = list(),
          interval_width = 1,
          connect_ascending = TRUE,
          ... ){

  # fill in missing
  if(missing(ybounds)) ybounds = cbind(y-interval_width*y_sd, y+interval_width*y_sd)
  if(is.null(ylim)) ylim = range(ybounds,na.rm=TRUE)

  # Decide which to break up
  if( connect_ascending==TRUE ){
    if(length(x)==1){
      obs_order = 1
    }else{
      obs_order = 1
      for( xI in 2:length(x) ){
        if( x[xI]>x[xI-1] ){
          obs_order = c(obs_order,xI)
        }else{
          obs_order = c(obs_order,NA,xI)
        }
      }
    }
    x = x[obs_order]
    y = y[obs_order]
    ybounds = ybounds[obs_order,]
  }

  # Plot lines
  fn( y=y, x=x, ylim=ylim, ... )

  # Plot shading
  if( bounds_type=="whiskers" ){
    for(t in 1:length(y)){
      lines_args = combine_lists( input=bounds_args,
                   default=list(col="black", lty="solid", lwd=1, x=rep(x[t],2), y=ybounds[t,]) )
      do.call( what=lines, args=lines_args )
    }
  }
  if( bounds_type=="shading" ){
    # Combine entries
    polygon_args = combine_lists( input=bounds_args,
                   default=list(col="black", border=NA, lty="solid", lwd=1, x=c(x,rev(x)), y=c(ybounds[,1],rev(ybounds[,2]))) )
    #polygon(, col=col_bounds, border=border, lty=border_lty)
    do.call( what=polygon, args=polygon_args )
  }
}
