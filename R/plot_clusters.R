
#' Plot spatial clusters
#'
#' \code{plot_clusters} plots the clusters on a map
#'
#' Code to plot clusters to look for biogeographic shifts, as well as plotting
#'      averages across all sites-times assigned to a each cluster, for use when
#'      interpreting each cluster.
#'
#' @inheritParams fastcluster::hclust.vector
#' @inheritParams stats::cutree
#'
#' @param var_name Name of object from \code{fit$Report} that is used.  Only
#'        implemented options are \code{"D_gct","Omega1_gc","Omega2_gc","Epsilon1_gct","Epsilon2_gct"}.
#' @param transform_var function to apply to \code{fit$Report[[var_name]]} prior to clustering.
#'        I recommend using log-transform for density, and otherwise using \code{transform_var=identity}.
#' @param method Distance metric. Default \code{method="ward"} is very fast,
#'        but can instead use \code{method="bcdist"} for small problems which calculates
#'        Bray-Curtist dissimilarity using \code{\link[ecodist]{bcdist}} and then applies
#'        Ward clustering
#' @param replace_Inf_with_NA Boolean whether to replace \code{Inf} or \code{-Inf} values
#'        with \code{NA} prior to clustering, as useful sometimes when \code{var_name="D_gct"},
#'        \code{transform_var=log} and replacing nonencounters with zero.
#' @param map_list output from \code{\link{make_map_info output}}
#'
#' @references For details regarding spatial clustering see \url{https://doi.org/10.1093/icesjms/fsac007}
#' @export
plot_clusters <-
function( fit,
          var_name = "D_gct",
          transform_var = log,
          k = 4,
          method = "ward",
          year_labels = fit$year_labels,
          category_names = fit$category_names,
          map_list = NULL,
          working_dir = getwd(),
          file_name = paste0("Class-",var_name),
          file_name2 = paste0("Class-",var_name,"-averages"),
          replace_Inf_with_NA = TRUE,
          size_threshold = 100000,
          col = viridisLite::viridis,
          yaxis_log = TRUE,
          ... ){

  # Informative error
  if( !(var_name %in% c("D_gct","Omega1_gc","Omega2_gc","Epsilon1_gct","Epsilon2_gct")) ){
    stop("Check `var_name`")
  }

  # Change labels
  fit$Report = amend_output( fit,
                             year_labels = year_labels,
                             category_names = category_names )

  # Make map_list if necessary
  if( missing(map_list) ){
    map_list = make_map_info( Region = fit$settings$Region,
                              spatial_list = fit$spatial_list,
                              Extrapolation_List = fit$extrapolation_list )
  }

  # Extract object
  Y_gct = transform_var( strip_units(fit$Report[[var_name]]) )
  if( length(dim(Y_gct))==2 ){
    Y_gct = Y_gct %o% array(1,dimnames=list("Time"=1))
  }

  # Change Inf e.g., from log(0) to NA
  if( replace_Inf_with_NA==TRUE ){
    Y_gct = ifelse( abs(Y_gct)==Inf, NA, Y_gct )
  }

  # Change shape
  Y_z = reshape2:::melt.array( data=Y_gct, varnames=names(dimnames(Y_gct)) )
  Y_zc = reshape2::acast(Y_z, formula = Time + Site ~ Category )

  # Remove NAs prior to clustering
  which_NA = which( apply(Y_zc,MARGIN=1,FUN=function(vec){any(is.na(vec))}) )
  which_notNA = setdiff( 1:nrow(Y_zc), which_NA )
  Yprime_zc = Y_zc[ which_notNA,,drop=FALSE ]

  # Warnings
  if( nrow(Yprime_zc) > size_threshold ){
    warning("Skipping `plot_clusters` ... it will likely not work due to large size")
    return( list("Y_zc"=Y_zc) )
  }
  if( nrow(Yprime_zc) > (size_threshold/10) ) warning("`plot_clusters` will go slowly due to large size")

  # Apply clustering
  if( method == "bcdist" ){
    # Option 1 -- breaks with large sample size
    Dist_zz = ecodist::bcdist(Yprime_zc)  # dist or ecodist::bcdist
    # Dist_zz = dist(Y_zc)
    Hclust = hclust(Dist_zz, method="ward.D2" )
    # Option 2 -- memory issues
    #cluster::agnes( Y_zc )
  }else{
    # Option 3 -- faster and lower memory
    Hclust = fastcluster::hclust.vector( Yprime_zc, method=method )
  }

  # Cut and convert shape
  Classprime_z = cutree( Hclust, k = k )

  # Add back to NAs
  Class_z = rep(NA, nrow(Y_zc))
  Class_z[which_notNA] = Classprime_z

  # back-transform shape
  if( prod(dim(Y_gct)[c(1,3)]) != length(Class_z) ) stop("Check `plot_clusters`")
  Class_gt = array(Class_z, dim=dim(Y_gct)[c(1,3)], dimnames=dimnames(Y_gct)[c(1,3)])

  # Make plot
  plot_variable(
    Y_gt = Class_gt,
    map_list = map_list,
    file_name = file_name,
    working_dir = working_dir,
    #format = format,
    panel_labels = colnames(Class_gt),
    col = col,
    ...
  )

  # Plot cluster memberships
  Ybar_kc = apply( Y_zc, MARGIN=2, FUN=function(y_z,class_z){tapply(y_z,INDEX=class_z,FUN=mean)}, class_z=Class_z )
  if(yaxis_log==TRUE){f=exp}else{f=identity}
  png( file=file.path(working_dir,paste0(file_name2,".png")), width=6, height=6, units="in", res=200 )
    par( mar=c(7,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
    matplot( y = f(t(Ybar_kc)),
             #x = factor(colnames(Ybar_kc)),
             #las = 2,
             xaxt = "n",
             col = col(k),
             type = "l",
             lwd = 2,
             xlab = "",
             ylab = "Cluster average",
             lty = "solid",
             log = ifelse(yaxis_log==TRUE,"y","") )
    axis( side = 1,
          at = seq_len(ncol(Ybar_kc)),
          labels = colnames(Ybar_kc),
          las = 2 )
  dev.off()

  # Return stuff
  Return = list("Y_zc"=Y_zc, "Class_gt"=Class_gt, "Ybar_kc"=Ybar_kc)
  return( invisible(Return) )
}
