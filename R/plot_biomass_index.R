
#' @title
#' Plot index of abundance
#'
#' @description
#' \code{plot_biomass_index} plots an index proportional to population abundance
#'
#' @inheritParams plot_maps
#' @param DirName Directory for saving plot and table
#' @param PlotName Name for plot
#' @param interval_width width for confidence intervals
#' @param strata_names names for spatial strata
#' @param category_names names for categories (if using package \code{`VAST`})
#' @param use_biascorr Boolean, whether to use bias-corrected estimates if available
#' @param plot_legend Add legend for labelling colors
#' @param plot_log Boolean, whether to plot y-axis in log-scale
#' @param width plot width in inches
#' @param height plot height in inches
#' @param TmbData Formatted data inputs, from \code{\link[VAST]{make_data}}
#' @param ... Other inputs to `par()`
#'
#' @return Return Tagged list of output
#' \describe{
#'   \item{Table}{table of index estimates by stratum and year, e.g., for including in an assessment model}
#' }
#'
#' @references For details regarding spatio-temporal index standardization see \url{https://doi.org/10.1093/icesjms/fsu243}
#' @references For details regarding fishing mortality and biomass reference points in models generating those, see \url{https://doi.org/10.1111/faf.12398}
#' @export
plot_biomass_index <-
function( fit,
          DirName = getwd(),
          PlotName = "Index",
          interval_width = 1,
          years_to_plot = NULL,
          category_names = NULL,
          year_labels = NULL,
          strata_names = NULL,
          use_biascorr = TRUE,
          plot_legend = TRUE,
          total_area_km2 = NULL,
          plot_log = FALSE,
          width = NULL,
          height = NULL,
          create_covariance_table = FALSE,
          Yrange = c(ifelse(plot_log==TRUE,NA,0),NA),
          TmbData = fit$data_list,
          Sdreport = fit$parameter_estimates$SD,
          extrapolation_list = fit$extrapolation_list,
          ... ){

  # Which parameters
  if( "ln_Index_ctl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # VAST Version < 2.0.0
    index_name = "Index_ctl"
    log_index_name = "ln_Index_ctl"
  }else if( "ln_Index_cyl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # VAST Version >= 2.0.0
    index_name = "Index_cyl"
    log_index_name = "ln_Index_cyl"
    TmbData[["n_t"]] = nrow(TmbData[["t_yz"]])
  }else{
    stop("`plot_biomass_index` is not compatible with your version")
  }

  # Add t_i if missing (e.g., from VAST V2.8.0 through V9.3.0)
  if( !("t_i" %in% names(TmbData)) ){
    TmbData$t_i = TmbData$t_iz[,1]
  }

  # Check and implement units and labels
  fit$Report = amend_output( fit = fit,
                             year_labels = year_labels,
                             category_names = category_names,
                             strata_names = strata_names,
                             extrapolation_list = extrapolation_list )

  # Informative errors
  if(is.null(Sdreport)) stop("Sdreport is NULL; please provide Sdreport")

  # Add in t_yz if missing (e.g., from earlier version of VAST, or SpatialDeltaGLMM)
  if( !("t_yz" %in% names(TmbData)) ){
    TmbData$t_yz = matrix(1:TmbData$n_t - 1, ncol=1)
  }

  # Fill in missing
  mfrow = c( ceiling(sqrt(TmbData$n_c)), ceiling(TmbData$n_c/ceiling(sqrt(TmbData$n_c))) )
  if( is.null(width)) width = mfrow[2] * 3
  if( is.null(height)) height = mfrow[1] * 3
  if( is.null(years_to_plot) ) years_to_plot = 1:TmbData$n_t

  # Logical check
  if( "unbiased"%in%names(Sdreport) ){
    if( all(is.na(Sdreport$unbiased$value)) ){
      stop("You appear to be using bias-correction, but all values are NA. Please report problem to package author.")
    }
  }

  # Defaults
  if( "treat_nonencounter_as_zero" %in% names(TmbData$Options_list$Options) ){
    treat_missing_as_zero = TmbData$Options_list$Options["treat_nonencounter_as_zero"]
  }else{
    treat_missing_as_zero = FALSE
  }

  # Objects
  # Could be moved to amend_output in future versions
  SD = TMB::summary.sdreport(Sdreport)
  if( !"report" %in% names(as.list(args(TMB:::as.list.sdreport))) ){
    warning( "package `TMB` should be updated to easily access standard errors")
  }
  par_SE = TMB:::as.list.sdreport( Sdreport, what="Std. Error", report=TRUE )
  par_hat = TMB:::as.list.sdreport( Sdreport, what="Estimate", report=TRUE )
  if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
    par_biascorrect = TMB:::as.list.sdreport( Sdreport, what="Est. (bias.correct)", report=TRUE )
    for( int in seq_len(length(par_hat)) ){
      par_hat[[int]] = ifelse( is.na(par_biascorrect[[int]]), par_hat[[int]], par_biascorrect[[int]] )
    }
  }

  # Fix at zeros any years-category combinations with no data
  # Could be moved to amend_output in future versions ... but requires moving par_hat and par_SE to amend_output
  if( treat_missing_as_zero==TRUE ){
    # Determine year-category pairs with no data
    Num_ctl = abind::adrop(TmbData$Options_list$metadata_ctz[,,'num_notna',drop=FALSE], drop=3) %o% rep(1,TmbData$n_l)
    # Replace values with 0 (estimate) and NA (standard error)
    par_hat[[index_name]] = ifelse(Num_ctl==0, 0, par_hat[[index_name]])
    par_SE[[index_name]] = ifelse(Num_ctl==0, 0, par_SE[[index_name]])
    par_hat[[log_index_name]] = ifelse(Num_ctl==0, 0, par_hat[[log_index_name]])
    par_SE[[log_index_name]] = ifelse(Num_ctl==0, 0, par_SE[[log_index_name]])
  }

  # Assign units after fixing values to zero
  for( int in seq_len(length(par_hat)) ){
    if( names(par_hat)[int] %in% names(fit$Report) ){
      dimnames(par_SE[[int]]) = dimnames(par_hat[[int]]) = dimnames(fit$Report[[names(par_hat)[int]]])
    }
    if("units" %in% class(fit$Report[[names(par_hat)[int]]])) units(par_hat[[int]]) = units(fit$Report[[names(par_hat)[int]]])
  }
  if( any(is.na(par_hat)) | any(is.na(par_SE)) ){
    stop( "Problem: Standard errors contain NAs")
  }

  # Plot biomass and Bratio
  Plot_suffix = ""
  if( "Bratio_ctl" %in% names(par_hat) ) Plot_suffix = c( Plot_suffix, "-Bratio" )
  for( plotI in 1:length(Plot_suffix) ){
    if( Plot_suffix[plotI]=="" ){
      Array_ctl = par_hat[[index_name]]
      log_Array_ctl = par_SE[[log_index_name]]
    }
    if( Plot_suffix[plotI]=="-Bratio" ){
      Array_ctl = par_hat[["Bratio_ctl"]]
      log_Array_ctl = par_SE[["ln_Bratio_ctl"]]
    }
    plot_index( Index_ctl = Array_ctl,
                sd_Index_ctl = log_Array_ctl,
                years_to_plot = years_to_plot,
                category_names = dimnames(Array_ctl)[[1]],
                year_labels = dimnames(Array_ctl)[[2]],
                strata_names = dimnames(Array_ctl)[[3]],
                DirName = DirName,
                PlotName = paste0(PlotName,Plot_suffix[plotI],".png"),
                interval_width = interval_width,
                width = width,
                height = height,
                xlab = "Year",
                ylab = make_unit_label( u = units(Array_ctl), lab = "Index", parse = FALSE ),
                scale = "log",
                plot_args = list("log" = ifelse(plot_log==TRUE,"y","")),
                Yrange = Yrange )
  }

  # Plot
  if( "Fratio_ct" %in% names(par_hat) ){
    Array_ct = par_hat[["Fratio_ct"]]
      Array_ct = ifelse( strip_units(Array_ct)==0, NA, Array_ct )
    plot_index( Index_ctl = Array_ct,
                sd_Index_ctl = par_SE[["Fratio_ct"]],
                years_to_plot = years_to_plot,
                category_names = dimnames(Array_ctl)[[1]],
                year_labels = dimnames(Array_ctl)[[2]],
                strata_names = dimnames(Array_ctl)[[3]],
                DirName = DirName,
                PlotName = paste0(PlotName,"-Fratio.png"),
                scale = "uniform",
                interval_width = interval_width,
                width = width,
                height = height,
                xlab = "Year",
                ylab = "Fishing ratio" )
  }

  # Plot stock status
  if( all(c("Fratio_ct","Bratio_ctl") %in% names(par_hat)) ){
    Par = list( mar=c(2,2,1,0), mgp=c(2,0.5,0), tck=-0.02, yaxs="i", oma=c(1,2,0,0), mfrow=mfrow, ... )
    Col = colorRampPalette(colors=c("blue","purple","red"))
    png( file=file.path(DirName,paste0(PlotName,"-Status.png")), width=width, height=height, res=200, units="in")
      par( Par )
      Array1_ct = abind::abind( "Estimate"=matrix(par_hat[["Bratio_ctl"]][,,1],nrow=TmbData$n_c,dimnames=dimnames(par_hat[["Bratio_ctl"]])[1:2]), "Std. Error"=matrix(par_SE[["Bratio_ctl"]][,,1],nrow=TmbData$n_c,dimnames=dimnames(par_hat[["Bratio_ctl"]])[1:2]), along=3 )
      Array1_ct = ifelse( Array1_ct==0, NA, Array1_ct )
      Array2_ct = abind::abind( "Estimate"=par_hat[["Fratio_ct"]], "Std. Error"=par_SE[["Fratio_ct"]], along=3 )
      Array2_ct = ifelse( Array2_ct==0, NA, Array2_ct )
      for( cI in 1:TmbData$n_c ){
        # Calculate y-axis limits
        Xlim = c(0, max(1, Array1_ct[cI,years_to_plot,'Estimate']%o%c(1,1) + Array1_ct[cI,years_to_plot,'Std. Error']%o%c(-interval_width,interval_width),na.rm=TRUE) )
        Ylim = c(0, max(2, Array2_ct[cI,years_to_plot,'Estimate']%o%c(1,1) + Array2_ct[cI,years_to_plot,'Std. Error']%o%c(-interval_width,interval_width),na.rm=TRUE) )
        # Plot stuff
        plot(1, type="n", xlim=Xlim, ylim=Ylim, xlab="", ylab="", main=ifelse(TmbData$n_c>1,category_names[cI],"") )
        points( x=Array1_ct[cI,years_to_plot,'Estimate'], y=Array2_ct[cI,years_to_plot,'Estimate'], col=Col(dim(Array1_ct)[2])[years_to_plot] )
        for( tI in years_to_plot ){
          lines( x=rep(Array1_ct[cI,tI,'Estimate'],2), y=Array2_ct[cI,tI,'Estimate']+Array2_ct[cI,tI,'Std. Error']*c(-interval_width,interval_width), col=Col(dim(Array1_ct)[2])[tI] )
          lines( x=Array1_ct[cI,tI,'Estimate']+Array1_ct[cI,tI,'Std. Error']*c(-interval_width,interval_width), y=rep(Array2_ct[cI,tI,'Estimate'],2), col=Col(dim(Array1_ct)[2])[tI] )
        }
        abline( v=0.4, lty="dotted" )
        abline( h=1, lty="dotted" )
      }        
      legend( "topright", bty="n", fill=c(Col(dim(Array1_ct)[2])[years_to_plot[1]],Col(dim(Array1_ct)[2])[rev(years_to_plot)[1]]), legend=c(dimnames(Array1_ct)[[2]][years_to_plot[1]],dimnames(Array1_ct)[[2]][rev(years_to_plot)[1]]) )
      mtext( side=1:2, text=c("Biomass relative to unfished","Fishing relative to F_40%"), outer=TRUE, line=c(0,0) )
    dev.off()
  }

  # Write to file
  Table = cbind( expand.grid(dimnames(par_hat[[index_name]])),
    "Units" = make_unit_label( u=units(par_hat[[index_name]]), lab="", parse=FALSE ),
    "Estimate" = as.vector(par_hat[[index_name]]),
    "Std. Error for Estimate" = as.vector(par_SE[[index_name]]),
    "Std. Error for ln(Estimate)" = as.vector(par_SE[[log_index_name]]) )
  write.csv( Table, file=file.path(DirName,"Index.csv"), row.names=FALSE)

  # Return stuff
    # Necessary to provide "log_Index_ctl" and "Index_ctl" for use in calculate_proportion, which has been fixed for zeros here
  Return = list( "Table"=Table,
    "log_Index_ctl" = abind::abind("Estimate"=par_hat[[log_index_name]], "Std. Error"=par_SE[[log_index_name]], along=4),
    "Index_ctl" = abind::abind("Estimate"=par_hat[[index_name]], "Std. Error"=par_SE[[index_name]], along=4) )

  # Extract and save covariance
  if( "cov"%in%names(Sdreport) & create_covariance_table==TRUE ){
    DF = expand.grid(dimnames(par_hat[[index_name]]))
    Which = which( names(Sdreport$value)==index_name )
    Cov = Sdreport$cov[Which,Which]
    Corr = cov2cor(Cov) - diag(nrow(Cov))
    rowcolDF = cbind( "RowNum"=row(Corr)[lower.tri(Corr,diag=TRUE)], "ColNum"=col(Corr)[lower.tri(Corr,diag=TRUE)] )
    Table = cbind( DF[rowcolDF[,'ColNum'],], DF[rowcolDF[,'RowNum'],] )
    colnames(Table) = paste0(colnames(Table), rep(c(1,2),each=3))
    Table = cbind( Table, "Correlation"=cov2cor(Cov)[lower.tri(Corr,diag=TRUE)], "Covariance"=Cov[lower.tri(Corr,diag=TRUE)] )
    Table = cbind( Table, "Index1"=Index_ctl[as.matrix(cbind(DF[rowcolDF[,'ColNum'],],1))], "Index2"=Index_ctl[as.matrix(cbind(DF[rowcolDF[,'RowNum'],],1))] )
    WhichZero = which( (Table[,'Index1']*Table[,'Index2']) == 0 )
    Table[WhichZero,c('Correlation','Covariance')] = 0
    Return = c( Return, "Table_of_estimated_covariance"=Table )
  }
  #if( !is.null(Bratio_ctl)) Return = c( Return, list("Bratio_ctl"=Bratio_ctl) )
  #if( !is.null(log_Bratio_ctl)) Return = c( Return, list("log_Bratio_ctl"=log_Bratio_ctl) )
  #if( !is.null(Fratio_ct)) Return = c( Return, list("Fratio_ct"=Fratio_ct) )

  return( invisible(Return) )
}
