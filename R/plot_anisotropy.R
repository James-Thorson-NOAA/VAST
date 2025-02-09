
#' Plot geometric anisotropy
#'
#' \code{plot_anisotropy} plots the degree of geometric anisotropy
#'
#' VAST estimates a matrix \code{H}, representing a linear transformation of coordinates (typically eastings-northings)
#' when computing a distance measure used to calculate correlations between sites. Geometric anisotropy represents the
#' tendancy for correlations to decline faster in some direction than others, as represented by \code{H}.
#' \code{plot_anistropy} visualizes the distance needed to achieve a correlation of approximately 10% from a location
#' centered at coordinates (0,0).  Therefore, an ellipse with a major (long) axis pointed northwest-southeast will
#' have correlations that decline slower along this axis then movement northeast-southwest. This ellipse is shown for
#' both linear predictors, which share the same \code{H} estimate but differ in their overall estimated decorrelation rate.
#'
#' @export
plot_anisotropy <-
function( Obj,
          FileName,
          ControlList = list("Width"=4, "Height"=5, "Res"=200, "Units"='in'),
          type = "ellipse",
          Report = Obj$report(),
          TmbData = Obj$env$data ){

  # Extract Options and Options_vec (depends upon version)
  if( all(c("Options","Options_vec") %in% names(TmbData)) ){
    Options_vec = TmbData$Options_vec
    Options = TmbData$Options
  }
  if( "Options_list" %in% names(TmbData) ){
    Options_vec = TmbData$Options_list$Options_vec
    Options = TmbData$Options_list$Options
  }

  # extract map
  Map = Obj$env$map
  Params = Obj$env$last.par.best

  if( Options_vec['Aniso']!=1 | all(c("logkappa1","logkappa2") %in% names(Map)) ){
    message("Skipping plot of geometric anisotropy because it has been turned off")
  }else{
    # Decomposition
    Eigen = eigen(Report$H)

    # Arrows
    if( type=="arrow" ){
      png(file=FileName, width=ControlList$Width, height=ControlList$Height, res=ControlList$Res, units=ControlList$Units)
        par( mar=c(2,2,0,0), mgp=c(1.5,0.5,0), tck=-0.02)
        plot( 1, type="n", xlim=c(-1,1)*max(Eigen$values), ylim=c(-1,1)*max(Eigen$values))
        arrows( x0=rep(0,2), y0=rep(0,2), x1=Eigen$vectors[1,]*Eigen$values, y1=Eigen$vectors[2,]*Eigen$values)
      dev.off()
    }

    # Ellipses
    if( type=="ellipse" ){
      rss = function(V) sqrt(sum(V[1]^2+V[2]^2))
      Major_1 = Minor_1 = Major_2 = Minor_2 = NA
      if( !("logkappa1"%in%names(Map)) || !is.na(Map$logkappa1) ){
        Major_1 = Eigen$vectors[,1]*Eigen$values[1] * Report$Range_raw1
        Minor_1 = Eigen$vectors[,2]*Eigen$values[2] * Report$Range_raw1
      }
      if( !("logkappa2"%in%names(Map)) || !is.na(Map$logkappa2) ){
        Major_2 = Eigen$vectors[,1]*Eigen$values[1] * Report$Range_raw2
        Minor_2 = Eigen$vectors[,2]*Eigen$values[2] * Report$Range_raw2
      }
      png(file=FileName, width=ControlList$Width, height=ControlList$Height, res=ControlList$Res, units=ControlList$Units)
        par( mar=c(3,3,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
        Range = 1.1 * c(-1,1) * max(abs( cbind(Major_1,Minor_1, Major_2,Minor_2) ),na.rm=TRUE)
        plot( 1, type="n", xlim=Range, ylim=c(Range[1],Range[2]*1.2), xlab="", ylab="")
        if( !any(is.na(Major_1)) ){
          shape::plotellipse( rx=rss(Major_1), ry=rss(Minor_1), angle=-1*(atan(Major_1[1]/Major_1[2])/(2*pi)*360-90), lcol="green", lty="solid")
        }
        if( !any(is.na(Major_2)) ){
          shape::plotellipse( rx=rss(Major_2), ry=rss(Minor_2), angle=-1*(atan(Major_2[1]/Major_2[2])/(2*pi)*360-90), lcol="black", lty="solid")
        }
        title( "Distance at 10% correlation" )
        mtext(side=1, outer=FALSE, line=2, text="Eastings (km.)")
        mtext(side=2, outer=FALSE, line=2, text="Northings (km.)")
        legend( "top", legend=c("1st linear predictor","2nd linear predictor"), fill=c("green","black"), bty="n")
        #abline( h=0, v=0, lty="dotted")
      dev.off()
    }
  }
}
