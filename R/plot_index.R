
#' Plot an index
#'
#' \code{plot_index} takes output from a VAST run and plots a panel figure of time-series estimates
#'
#' @inheritParams plot_biomass_index
#' @param Index_ctl A matrix or array of time-series estimates for multiple categories \code{c}, years \code{t}, and strata \code{l}
#' @param sd_Index_ctl A matrix or array of variances for each estimate
#' @param Yrange lower and upper bound for left-hand y-axis, corresponding to input \code{Index_ctl}
#'        (use \code{Yrange[1]=NA} and/or \code{Yrange[2]=NA} for using the lower and upper bound of estimate intervals)
#' @param Y2range lower and upper bound for right-hand y-axis, corresponding to input \code{SampleSize_ctz} (see Yrange for more info)
#' @param SampleSize_ctz optional array of sample sizes for each category and year to be plotted on each panel
#' @param plot_args additional arguments to pass to \code{plot}
#' @param plot_lines_args additional arguments to pass to \code{plot_lines}
#' @param ... list of settings to pass to \code{par} when making plot
#'
#' @return NULL
#'
#' @export
plot_index <-
function( Index_ctl,
          sd_Index_ctl = array(0,dim(Index_ctl)),
          year_labels = NULL,
          years_to_plot = NULL,
          strata_names = NULL,
          category_names = NULL,
          scale = "uniform",
          plot_legend = NULL,
          DirName = getwd(),
          PlotName = "Index.png",
          interval_width = 1,
          width = NULL,
          height = NULL,
          xlab = "Year",
          ylab = "Index",
          bounds_type = "whiskers",
          col = NULL,
          col_bounds = NULL,
          Yrange = c(0,NA),
          type = "b",
          plot_lines_args = list(),
          plot_args = list(),
          SampleSize_ctz = NULL,
          Y2range = c(0,NA),
          y2lab = "",
          add = FALSE,
          ... ){

  # Local function
  plot_lines <-
  function( x,
            y,
            ybounds,
            fn = lines,
            col_bounds = "black",
            bounds_type = "whiskers",
            border = NA,
            border_lty = "solid",
            lwd_bounds = 1,
            ... ){

    # Function still used in plot_index
    #warning( "`plot_lines` is soft-deprecated" )

    fn( y=y, x=x, ... )
    if( bounds_type=="whiskers" ){
      for(t in 1:length(y)){
        lines( x=rep(x[t],2), y=ybounds[t,], col=col_bounds, lty=border_lty, lwd=lwd_bounds)
      }
    }
    if( bounds_type=="shading" ){
      polygon( x=c(x,rev(x)), y=c(ybounds[,1],rev(ybounds[,2])), col=col_bounds, border=border, lty=border_lty)
    }
  }

  # Change inputs
  if( length(dim(Index_ctl))==length(dim(sd_Index_ctl)) ){
    if( length(dim(Index_ctl))==2 ){
      Index_ctl = Index_ctl %o% rep(1,1)
      sd_Index_ctl = sd_Index_ctl %o% rep(1,1)
    }
  }else{
    stop("Mismatch in dimensions for `Index_ctl` and `sd_Index_ctl` in `plot_index`")
  }
  n_categories = dim(Index_ctl)[1]
  n_years = dim(Index_ctl)[2]
  n_strata = dim(Index_ctl)[3]
  mfrow = c( ceiling(sqrt(n_categories)), ceiling(n_categories/ceiling(sqrt(n_categories))) )
  if( !is.null(SampleSize_ctz) ){
    if( !all( dim(SampleSize_ctz)[1:2] == dim(Index_ctl)[1:2] ) ){
      stop("Check input `SampleSize_ctz`")
    }
  }

  if(any(is.na(as.numeric(year_labels)))){
    x_Years = 1:length(year_labels)
  }else{
    x_Years = as.numeric(year_labels)
  }
  #if( all(is.numeric(Year_Set)) ){
  #  year_names = Year_Set
  #}else{
  #  year_names = Year_Set
  #  Year_Set = 1:length(Year_Set)
  #}

  Pretty = function(vec){
    Return = pretty(vec)
    Return = Return[which(Return %in% vec)]
    return(Return)
  }

  # Defaults
  if( is.null(col)) col = rainbow(n_strata)
  if( is.null(col_bounds)) col_bounds = rainbow(n_strata)
  if( is.null(width)) width = mfrow[2] * 3
  if( is.null(height)) height = mfrow[1] * 3

  # Fill in missing
  if( is.null(year_labels) ) year_labels = 1:n_years
  if( is.null(years_to_plot) ) years_to_plot = 1:n_years
  if( is.null(strata_names) ) strata_names = 1:n_strata
  if( is.null(category_names) ) category_names = 1:n_categories
  if( is.null(plot_legend)) plot_legend = ifelse(n_strata>1, TRUE, FALSE)

  # Plot
  Par = combine_lists( default=list(mar=c(2,2,1,0),mgp=c(2,0.5,0),tck=-0.02,yaxs="i",oma=c(2,2,0,0),mfrow=mfrow), input=list(...) )
  if(!is.na(PlotName)){
    png( file=file.path(DirName,PlotName), width=width, height=height, res=200, units="in")  # paste0(DirName,ifelse(DirName=="","","/"),PlotName)
    on.exit( dev.off() )
  }
  if(add==FALSE) par( Par )
  for( z1 in 1:n_categories ){
    # Calculate y-axis limits
    if(scale=="uniform") Ylim = range(Index_ctl[z1,years_to_plot,]%o%c(1,1) + sd_Index_ctl[z1,years_to_plot,]%o%c(-interval_width,interval_width)*1.05, na.rm=TRUE)
    if(scale=="log") Ylim = range(Index_ctl[z1,years_to_plot,]%o%c(1,1) * exp(sd_Index_ctl[z1,years_to_plot,]%o%c(-interval_width,interval_width)*1.05), na.rm=TRUE)
    Ylim = ifelse( is.na(Yrange), Ylim, Yrange )
    Xlim = range(x_Years[years_to_plot]) + c(-1,1) * diff(range(x_Years[years_to_plot]))/20
    # Plot stuff
    plot_inputs = combine_lists( default=list(1, type="n", xlim=Xlim, ylim=Ylim, xlab="", ylab="", main=ifelse(n_categories>1,category_names[z1],""), xaxt="n"),
      input=plot_args )
    do.call( what=plot, args=plot_inputs )
    for(z3 in 1:n_strata){
      if(scale=="uniform") ybounds = Index_ctl[z1,years_to_plot,z3]%o%c(1,1) + sd_Index_ctl[z1,years_to_plot,z3]%o%c(-interval_width,interval_width)
      if(scale=="log") ybounds = Index_ctl[z1,years_to_plot,z3]%o%c(1,1) * exp(sd_Index_ctl[z1,years_to_plot,z3]%o%c(-interval_width,interval_width))
      if( n_strata==1 ) x_offset = 0
      if( n_strata>=2 ) x_offset = seq(-0.1, 0.1, length=n_strata)[z3]
      plot_lines_defaults = list( y=Index_ctl[z1,years_to_plot,z3], x=x_Years[years_to_plot]+x_offset, ybounds=ybounds, type=type, col=col[z3], col_bounds=col_bounds[z3],
        ylim=Ylim, bounds_type=bounds_type )
      plot_lines_inputs = combine_lists( default=plot_lines_defaults, input=plot_lines_args )
      do.call( what=plot_lines, args=plot_lines_inputs )
    }
    # Plot lines for sample size
    if( !is.null(SampleSize_ctz) ){
      Y2lim = c(1, 1.2) * range(SampleSize_ctz[z1,years_to_plot,], na.rm=TRUE)
      Y2lim = ifelse( is.na(Y2range), Y2lim, Y2range )
      Labels = pretty(Y2lim)
      At = Labels / diff(range(Y2lim,na.rm=TRUE)) * diff(Ylim) + Ylim[1]
      axis( side=4, at=At, labels=Labels )
      for( z3 in 1:dim(SampleSize_ctz)[3] ){
        Y = SampleSize_ctz[z1,years_to_plot,z3] / diff(range(Y2lim,na.rm=TRUE)) * diff(Ylim) + Ylim[1]
        lines( x=x_Years, y=Y, col=col[z3], lwd=3, lty="dotted" )
        #points( x=year_labels, y=Y, col=col[z3], cex=1.5 )
      }
    }
    if(plot_legend==TRUE){
      legend( "top", bty="n", fill=rainbow(n_strata), legend=as.character(strata_names), ncol=2 )
    }
    axis( 1, at=Pretty(x_Years[years_to_plot]), labels=year_labels[match(Pretty(x_Years[years_to_plot]),x_Years)] )
  }
  mtext( side=c(1,2,4), text=c(xlab,ylab,y2lab), outer=TRUE, line=c(0,0) )

  return(invisible(plot_lines_inputs))
}
