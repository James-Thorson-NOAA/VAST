


#' @title
#' Plot shifts in range edges
#'
#' @description
#' \code{plot_range_edge} plots range edges
#'
#' @inheritParams plot_biomass_index
#' @inheritParams sample_variable
#' @param working_dir Directory for plots
#' @param quantiles vector specifying quantiles to use for calculating range edges
#' @param calculate_relative_to_average Boolean, whether to calculate edge in UTM coordinates (default),
#'        or instead calculate relative to median across all years. The latter reduces standard errors,
#'        and is appropriate when checking significance for comparison across years for a single species.
#'        The former (default) is appropriate for checking significance for comparison across species.
#'
#' @references For details regarding multivariate index standardization and expansion see \url{https://doi.org/10.22541/au.160331933.33155622/v1}
#' @export
plot_range_edge <-
function( Sdreport,
          Obj,
          year_labels = NULL,
          years_to_plot = NULL,
          strata_names = NULL,
          category_names = NULL,
          working_dir = getwd(),
          quantiles = c(0.05,0.95),
          n_samples = 100,
          interval_width = 1,
          width = NULL,
          height = NULL,
          calculate_relative_to_average = FALSE,
          seed = 123456,
          ...){

  # Unpack
  Report = Obj$report()
  TmbData = Obj$env$data

  # Informative errors
  if(is.null(Sdreport)) stop("Sdreport is NULL; please provide Sdreport")
  if( any(quantiles<0) | any(quantiles>1) ) stop("Please provide `quantiles` between zero and one")
  if( all(TmbData$Z_gm==0) ) stop("Please re-run with 'Options['Calculate_Range']=TRUE' to calculate range edges")
  if( n_samples<10 ) stop("`n_samples` must be at least 10 for any chance of meaningful results for `plot_range_edge`")
  if( n_samples<100 ) warning("Package author recommends `n_samples`>=100 for `plot_range_edge`")

  # Which parameters
  if( "ln_Index_tl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # SpatialDeltaGLMM
    stop("Not implemente")
  }
  if( "ln_Index_ctl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # VAST Version < 2.0.0 or >= 3.6.0
    D_name = "D_gct"
  }
  if( "ln_Index_cyl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # VAST Version >= 2.0.0 and < 3.6.0
    TmbData[["n_t"]] = nrow(TmbData[["t_yz"]])
    D_name = "D_gcy"
  }

  # Default inputs
  if( is.null(year_labels)) year_labels = 1:TmbData$n_t
  if( is.null(years_to_plot) ) years_to_plot = 1:TmbData$n_t
  if( is.null(strata_names) ) strata_names = 1:TmbData$n_l
  if( is.null(category_names) ) category_names = 1:TmbData$n_c
  if( is.null(colnames(TmbData$Z_gm)) ){
    m_labels = paste0("axis",1:ncol(TmbData$Z_gm))
  }else{
    m_labels = colnames(TmbData$Z_gm)
  }

  ##### Local function
  D_gctr = sample_variable( Sdreport=Sdreport, Obj=Obj, variable_name=D_name, n_samples=n_samples, seed=seed )
  if(any(D_gctr==Inf)) stop("`sample_variable` in `plot_range_edge` is producing D=Inf; please use `n_samples=0` to avoid error")
  # If some are NA, check to see what variable is going haywire
  if( FALSE ){
    D_gctr = sample_variable( Sdreport=Sdreport, Obj=Obj, variable_name="L_epsilon2_cf", n_samples=n_samples, seed=seed )
  }

  # Calculate quantiles from observed and sampled densities D_gcy
  E_zctm = array(NA, dim=c(length(quantiles),dim(Report[[D_name]])[2:3],ncol(TmbData$Z_gm)) )
  E_zctmr = array(NA, dim=c(length(quantiles),dim(Report[[D_name]])[2:3],ncol(TmbData$Z_gm),n_samples) )
  Mean_cmr = array(NA, dim=c(dim(Report[[D_name]])[2],ncol(TmbData$Z_gm),n_samples) )
  prop_zctm = array(NA, dim=c(dim(Report[[D_name]])[1:3],ncol(TmbData$Z_gm)) )
  prop_zctmr = array(NA, dim=c(dim(Report[[D_name]])[1:3],ncol(TmbData$Z_gm),n_samples) )
  for( rI in 0:n_samples ){
  for( mI in 1:ncol(TmbData$Z_gm) ){
    order_g = order(TmbData$Z_gm[,mI], decreasing=FALSE)
    if(rI==0) prop_zctm[,,,mI] = apply( Report[[D_name]], MARGIN=2:3, FUN=function(vec){cumsum(vec[order_g])/sum(vec)} )
    if(rI>=0) prop_zctmr[,,,mI,rI] = apply( D_gctr[,,,rI,drop=FALSE], MARGIN=2:3, FUN=function(vec){cumsum(vec[order_g])/sum(vec)} )

    # Calculate edge
    for( cI in 1:dim(E_zctm)[2] ){
      if(rI>=1){
        if( calculate_relative_to_average==TRUE ){
          Mean_cmr[cI,mI,rI] = weighted.mean( as.vector(TmbData$Z_gm[,mI]%o%rep(1,dim(Report[[D_name]])[3])), w=as.vector(D_gctr[,cI,,rI]) )
        }else{
          Mean_cmr[cI,mI,rI] = 0
        }
      }
      for( zI in 1:dim(E_zctm)[1] ){
      for( tI in 1:dim(E_zctm)[3] ){
        if(rI==0){
          index_tmp = which.min( (prop_zctm[,cI,tI,mI]-quantiles[zI])^2 )
          E_zctm[zI,cI,tI,mI] = TmbData$Z_gm[order_g[index_tmp],mI]
        }
        if(rI>=1){
          index_tmp = which.min( (prop_zctmr[,cI,tI,mI,rI]-quantiles[zI])^2 )
          E_zctmr[zI,cI,tI,mI,rI] = TmbData$Z_gm[order_g[index_tmp],mI] - Mean_cmr[cI,mI,rI]
        }
      }}
    }
  }}
  SE_zctm = apply( E_zctmr, MARGIN=1:4, FUN=sd )
  Edge_zctm = abind::abind( "Estimate"=E_zctm, "Std. Error"=SE_zctm, along=5 )
  dimnames(Edge_zctm)[[1]] = paste0("quantile_",quantiles)

  # Plot cumulative distribution
  # matplot( x=Z_zm[,mI], y=prop_zcym[,1,,mI], type="l" )

  # Plot edges
  for( mI in 1:dim(Edge_zctm)[4] ){
    Index_zct = array(Edge_zctm[,,,mI,'Estimate'],dim(Edge_zctm)[1:3])
    sd_Index_zct = array(Edge_zctm[,,,mI,'Std. Error'],dim(Edge_zctm)[1:3])
    plot_index( Index_ctl = aperm(Index_zct,c(2,3,1)),
                sd_Index_ctl = aperm(sd_Index_zct,c(2,3,1)),
                year_labels = year_labels,
                years_to_plot = years_to_plot,
                strata_names = quantiles,
                category_names = category_names,
                DirName = working_dir,
                PlotName = paste0("RangeEdge_",m_labels[mI],".png"),
                Yrange = c(NA,NA),
                interval_width = interval_width,
                width = width,
                height = height,
                xlab = "Year",
                ylab = paste0("Quantiles (",m_labels[mI],")") )
  }

  # Return list of stuff
  Return = list( "year_labels"=year_labels, "Edge_zctm"=Edge_zctm )
  return( invisible(Return) )
}

