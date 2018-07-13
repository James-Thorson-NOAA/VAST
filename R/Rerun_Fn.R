
#' Explore counter-factual scenario
#'
#' \code{Rerun_Fn} re-builds a model while fixing a subset of parameters at zero
#'
#' @param parhat0 parameter values to use by default (presumably maximum-likelihood estimates)
#' @param turnoff_pars character-vector of parameters to turn off
#' @param loc_x location for each knot used to calculate center-of-gravity for counter-factual runs with colnames c('E_km','N_km','Lon','Lat')
#' @param cov_to_turnoff which covariates to turn off, as indicated by their order in \code{TmbData$X_xtp} (which only matters if "gamma1_ctp" or "gamma2_ctp" is in \code{turnoff_pars})
#' @param calculate_COG Boolean whether to calculate COG for each run
#' @param figname name for figure to plot density in counter-factual scenario
#' @inheritParams Build_TMB_Fn
#' @param MapDetails_List output from \code{FishStatsUtils::MapDetails_Fn}
#' @param year_set set of parameters to include
#' @param c_set set of categories to include
#' @param ... additional arguments passed to \code{VAST::Build_TMB_Fn}

#' @return Tagged list
#' \describe{
#'   \item{Report}{Report output for counter-factual run}
#'   \item{NewBuild_List}{Output from \code{VAST::Build_TMB_Fn} using counter-factual parameters}
#' }

#' @export
Rerun_Fn = function( parhat0, turnoff_pars, loc_x, cov_to_turnoff=1:dim(parhat[["gamma2_ctp"]])[3], calculate_COG=TRUE, figname=NULL,
  Map="generate", MapDetails_List=NULL, year_set=1:ncol(parhat0[["beta1_ct"]]), c_set=1:nrow(parhat0[["beta1_ct"]]), ... ){

  # Local function -- calculate center of gravity
  Calc_COG = function( z_x, B_xt ){
    # Set-up
    if( is.vector(z_x) ) z_x = matrix( z_x, ncol=1 )
    P_xt = B_xt / ( rep(1,nrow(B_xt)) %o% apply(B_xt, MARGIN=2, FUN=sum) )

    # Loop through dimensions
    COG_t = NULL
    for( i in 1:ncol(z_x)){
      COG_t = cbind( COG_t, apply( P_xt, MARGIN=2, FUN=weighted.mean, x=z_x[,i]))
    }
    colnames(COG_t) =colnames(z_x)

    return( COG_t )
  }
  # Local function -- plot density maps
  plot_density = function( Report, MapDetails_List, type, year_set, c_set ){
    if("D_xct" %in% names(Report)) D_xcy = Report$D_xct
    if("D_xcy" %in% names(Report)) D_xcy = Report$D_xcy
    if( type=="year" ){
      Dim = c("Nrow"=ceiling(sqrt(length(year_set)))); Dim = c(Dim,"Ncol"=ceiling(length(year_set)/Dim['Nrow']))
      for( i in 1:(length(c_set)-1) ){
        Mat = log(D_xcy[,i,])
        Zlim = range( log(D_xcy[,i,]) )
        FishStatsUtils:::PlotMap_Fn( MappingDetails=MapDetails_List[["MappingDetails"]], Mat=Mat, zlim=Zlim, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName="", Year_Set=year_set, Rescale=FALSE, Rotate=MapDetails_List[["Rotate"]], Format="", Res=NA, mfrow=Dim, mar=c(0,0,2,0), oma=c(3.5,3.5,4,0), Cex=MapDetails_List[["Cex"]], textmargin=textmargin, add=FALSE, pch=20, zone=MapDetails_List[["Zone"]] )
        mtext( side=3, outer=TRUE, text=paste0("Size bounds:",c_set[i]," to ",c_set[i+1]), cex=1.5, line=1 )
      }
    }
    if( type=="category" ){
      Dim = c("Nrow"=ceiling(sqrt(length(c_set)-1))); Dim = c(Dim,"Ncol"=ceiling((length(c_set)-1)/Dim['Nrow']))
      for( i in 1:length(year_set) ){
        Mat = log(D_xcy[,,i])
        Zlim = range( log(D_xcy[,,i]) )
        FishStatsUtils:::PlotMap_Fn( MappingDetails=MapDetails_List[["MappingDetails"]], Mat=Mat, zlim=Zlim, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName="", Year_Set=c_set[-1], Rescale=FALSE, Rotate=MapDetails_List[["Rotate"]], Format="", Res=NA, mfrow=Dim, mar=c(0,0,2,0), oma=c(3.5,3.5,4,0), Cex=MapDetails_List[["Cex"]], textmargin=textmargin, add=FALSE, pch=20, zone=MapDetails_List[["Zone"]] )
        mtext( side=3, outer=TRUE, text=paste0("Year:",year_set[i]), cex=1.5, line=1 )
      }
    }    #
  }

  # Check for errors
  if( !all(turnoff_pars %in% c("Epsiloninput1_sft","Epsiloninput2_sft","beta1_ct","beta2_ct","gamma1_ctp","gamma2_ctp")) ) stop("Some parameters in turnoff_pars are not as expected")
  if( !all(turnoff_pars %in% names(parhat0)) ) stop("Some parameters in turnoff_pars are not in parhat0")

  # Eliminate spatio-temporal variation and differences in age-composition
  parhat_mod = parhat0
  if( "Epsiloninput1_sft" %in% turnoff_pars ) parhat_mod[["Epsiloninput1_sft"]][] = 0
  if( "Epsiloninput2_sft" %in% turnoff_pars ) parhat_mod[["Epsiloninput2_sft"]][] = 0
  if( "beta1_ct" %in% turnoff_pars ) parhat_mod[["beta1_ct"]][] = outer( rowMeans(ParHat[["beta1_ct"]]), rep(1,ncol(ParHat[["beta1_ct"]])) )
  if( "beta2_ct" %in% turnoff_pars ) parhat_mod[["beta2_ct"]][] = outer( rowMeans(ParHat[["beta2_ct"]]), rep(1,ncol(ParHat[["beta2_ct"]])) )
  which_cov_to_turnoff = intersect(cov_to_turnoff, 1:dim(parhat_mod[["gamma1_ctp"]])[3])
  if( "gamma1_ctp" %in% turnoff_pars ) parhat_mod[["gamma1_ctp"]][,,which_cov_to_turnoff] = 0
  if( "gamma2_ctp" %in% turnoff_pars ) parhat_mod[["gamma2_ctp"]][,,which_cov_to_turnoff] = 0

  # Rebuild
  NewBuild_List = Build_TMB_Fn( "Parameters"=parhat_mod, Map=Map, ... )
  Obj = NewBuild_List[["Obj"]]
  Report = Obj$report( )
  Return = list("Report"=Report, "NewBuild_List"=NewBuild_List)

  # calculate COG
  if( calculate_COG==TRUE ){
    # Figure out column names to use
    ColNames = c('E_km','N_km')
    if( all(c('Lon','Lat') %in% colnames(loc_x)) ) ColNames = c(ColNames,c('Lon','Lat'))
    # Total
    B_xt = apply(Report$Index_xcyl[,,,1,drop=FALSE], MARGIN=c(1,3), FUN=sum )
    Return[["COG_t"]] = Calc_COG( z_x=loc_x[,ColNames], B_xt=B_xt )
    # By category
    Return[["COG_ct"]] = NULL
    for( cI in 1:dim(Report$Index_xcyl)[2] ){
      Return[["COG_ct"]] = abind::abind( Return[["COG_ct"]], Calc_COG(z_x=loc_x[,ColNames], B_xt=Report$Index_xcyl[,cI,,1]), along=3)
    }
    Return[["COG_ct"]] = aperm(Return[["COG_ct"]], c(3,1,2))
  }

  # Plot
  if( !is.null(figname) & !is.null(MapDetails_List) ){
    message( "Starting plot for counter-factual run" )
    ThorsonUtilities::save_fig( file=figname, width=8, height=12, type="pdf", onefile=TRUE )
      plot_density( Report=Report, MapDetails_List=MapDetails_List, type="year", year_set=Year_Set, c_set=c_set )
    dev.off()
  }

  # Return
  return( Return )
}
