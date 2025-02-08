



#' Copy of FishStatsUtils::calculate_knot_areas
#'
#' Included for continuity when using old scripts
#'
#' Please use \code{?FishStatsUtils::calculate_knot_areas} to see list of arguments
#' and usage
#' @param ... Arguments to be passed to \code{calculate_knot_areas}
#' @export
Calc_Polygon_Areas_and_Polygons_Fn = function( ... ){
  .Deprecated( new="FishStatsUtils::calculate_knot_areas" )
  calculate_knot_areas( ... )
}

#' Copy of FishStatsUtils::make_mesh
#'
#' Included for continuity when using old scripts
#'
#' Please use \code{?FishStatsUtils::make_mesh} to see list of arguments
#' and usage
#' @param ... Arguments to be passed to \code{make_mesh}
#' @export
Calc_Anisotropic_Mesh = function( ... ){
  .Deprecated( new="FishStatsUtils::make_mesh" )
  make_mesh( ... )
}

#' Plot variance of GMRF knots
#'
#' \code{map_hypervariance} Plot variance of GMRF knots
#'
#' @param report Report from TmbObj (e.g., TmbList[["Obj"]])
#' @param Spatial_List Output from FishStatsUtils::Spatial_Information_Fn()
#' @param method Choose whether to plot anisotropic or isotropic covariance matrix
map_hypervariance <- function(report, Spatial_List, method){
  solveSubset <- function(Q) {
    require(Matrix)
    require(TMB)
    L <- Cholesky(Q, super=TRUE, perm=TRUE)
    invQ <- .Call("tmb_invQ", L, PACKAGE = "TMB")
    iperm <- invPerm(L@perm + 1L)
    invQ[iperm, iperm, drop=TRUE]
  }

  if (method == "anisotropic") {
    mesh <- Spatial_List$MeshList$anisotropic_mesh
  } else if (method == "isotropic") {
    mesh <- Spatial_List$MeshList$isotropic_mesh
  }


  # Q1 is the covariance matrix for the GMRF of the encounter model
  # Q2 is the covariance matrix for the GMRF of the catch-rate model
  # Mesh is the knots associated with the GMRF

  diag_sigma <- matrix(data = NA, nrow = length(diag(report$Q1)), ncol = 2)
  for( i in 1:2 ){
    Sigma = solveSubset( report[[ c("Q1","Q2")[i] ]] )
    diag_sigma[,i] <- log(diag(Sigma))
  }

  boundary_id <- mesh$segm$bnd$idx
  df2plot0 <-
    data.frame(E_km = mesh$loc[,1],
               N_km = mesh$loc[,2],
               sigma_q1 = diag_sigma[,1],
               sigma_q2 = diag_sigma[,2],
               label = rep(c("Encounter GMRF", "Catch-rate GMRF"),
                           each = length(RelVar[,1])))
  df2plot <-
    tidyr::gather(df2plot0, variable, value, -E_km, -N_km, -label)

  p <-
    ggplot2::ggplot(df2plot, ggplot2::aes(x = E_km, y = N_km)) +
    ggplot2::geom_point(ggplot2::aes(color = value)) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::facet_wrap(~label) +
    ggplot2::ggtitle("Log-variance of GMRF knots")

  print(p)
}


#' Explore spatio-temporal covariance
#'
#' \code{summarize_covariance} plots and returns spatio-temporal covariance among categories
#'
#' @inheritParams plot_overdispersion
#' @param category_order integer-vector to re-order categories while plotting
#' @param category_names character-vector listing name for each category
#' @param plotdir directory (absolute path) for plots
#' @param figname character for first-part of figure names
#' @param plotTF integer-vector (length 4) specifying which covariances to plot (recommended: use \code{plotTF=FieldConfig})
#' @param plot_cor Boolean, whether to plot correlation or covariance
#' @param mgp passed to \code{par}
#' @param tck passed to \code{par}
#' @param oma passed to \code{par}
#' @param ... passed to \code{ThorsonUtilities::save_fig}
summarize_covariance <-
function( Report,
          Data,
          ParHat,
          SD=NULL,
          category_order = 1:Data$n_c,
          category_names = 1:Data$n_c,
          plotdir = paste0(getwd(),"/"),
          figname = "Cov",
          plotTF = NULL,
          plot_cor = TRUE,
          mgp = c(2,0.5,0),
          tck = -0.02,
          oma = c(0,5,2,0),
          ...){

  #
  .Deprecated( new="FishStatsUtils::plot_similarities" )

  # Adds intercept defaults to FieldConfig if missing
  if( is.vector(Data[["FieldConfig"]]) && length(Data[["FieldConfig"]])==4 ){
    Data[["FieldConfig"]] = rbind( matrix(Data[["FieldConfig"]],ncol=2,dimnames=list(c("Omega","Epsilon"),c("Component_1","Component_2"))), "Beta"=c("Beta1"=-2,"Beta2"=-2) )
  }else{
    if( !is.matrix(Data[["FieldConfig"]]) || !all(dim(Data[["FieldConfig"]])==c(3,2)) ){
      stop("`FieldConfig` has the wrong dimensions in `Summarize_Covariance`")
    }
  }

  # Add default for plotTF, or coerce Data$FieldConfig to a vector
  if( is.null(plotTF) ){
    plotTF = as.vector( Data[["FieldConfig"]]>0 )
  }else{
    plotTF = as.vector(plotTF)
  }

  # Object to return
  Return = list()

  # Extract
  for(i in which(Data[["FieldConfig"]]>=0) ){
    Par_name = c("omega1", "epsilon1", "beta1", "omega2", "epsilon2", "beta2")[i]
    L_name = paste0("L_",Par_name,"_z")

    # Extract estimates and standard errors
    if( !is.null(SD) ){
      # Object to build
      sd_summary = summary(SD)
      Slot_name = paste0("lowercov_uppercor_",Par_name)
      if( Slot_name %in% rownames(sd_summary) ){
        # Extract covariances
        Cor = Cov = Mat = ThorsonUtilities::Extract_SE( SD=SD, parname=Slot_name, columns=1:2, Dim=c(Data$n_c,Data$n_c) )
        dimnames(Cor) = dimnames(Cov) = list( category_names, category_names, c("Estimate","Std.Error") )
        # Cor
        Cor[,,1][lower.tri(Cor[,,1])] = t(Mat[,,1])[lower.tri(Mat[,,1])]
        diag(Cor[,,1]) = 1
        Cor[,,2][lower.tri(Cor[,,2])] = t(Mat[,,2])[lower.tri(Mat[,,2])]
        diag(Cor[,,2]) = NA
        # Cov
        Cov[,,1][upper.tri(Cov[,,1])] = t(Mat[,,1])[upper.tri(Mat[,,1])]
        Cov[,,2][upper.tri(Cov[,,2])] = t(Mat[,,2])[upper.tri(Mat[,,2])]
      }else{
        Cov = Cor = NULL
      }
    }else{
      Cov = Cor = NULL
    }

    # Extract estimates
    if( is.null(Cov) | is.null(Cor) ){
      Cov = Cor = array( NA, dim=c(Data$n_c,Data$n_c,2), dimnames=list(category_names,category_names,c("Estimate","Std.Error") ) )
      Cov[,,'Estimate'] = calc_cov( L_z=ParHat[[L_name]], n_f=as.vector(Data[["FieldConfig"]])[i], n_c=Data$n_c )
      Cor[,,'Estimate'] = cov2cor( Cov[,,'Estimate'] )
    }                       #

    # Add to return
    List = list( Cor, Cov )
    names(List) = paste0(c("Cor_","Cov_"), Par_name)
    Return = c( Return, List )
  }

  # Plot covariances
  if( !is.null(figname) ){
    # Work out dimensions
    Dim = c(3,2)
    if( sum(ifelse(plotTF>0,1,0))==1 ) Dim = c(1,1)
    if( all(ifelse(plotTF>0,1,0)==c(1,1,0,0,0,0)) | all(ifelse(plotTF>0,1,0)==c(0,0,1,1,0,0)) ) Dim=c(1,2)
    if( all(ifelse(plotTF>0,1,0)==c(1,0,1,0,0,0)) | all(ifelse(plotTF>0,1,0)==c(0,1,0,1,0,0)) ) Dim=c(2,1)

    # Conversion function
    if(plot_cor==TRUE){
      convert = function(Cov) ifelse(is.na(cov2cor(Cov)),0,cov2cor(Cov))
    }else{
      convert = function(Cov) ifelse(is.na(Cov),0,Cov)
    }

    # Plot analytic
    ThorsonUtilities::save_fig( file=paste0(plotdir,figname,"--Analytic.png"), width=Dim[2]*4+1, height=Dim[1]*4, ... )
      par(mfrow=Dim, mar=c(0,1,1,0), mgp=mgp, tck=tck, oma=oma)
      for(i in 1:6 ){      #
        if( i %in% which(plotTF>0) ){
          Cov_cc = calc_cov( L_z=ParHat[c('L_omega1_z','L_epsilon1_z','L_beta1_z','L_omega2_z','L_epsilon2_z','L_beta2_z')][[i]], n_f=as.vector(Data[["FieldConfig"]])[i], n_c=Data$n_c )
          plot_cov( Cov=convert(Cov_cc)[category_order,category_order], names=list(category_names[category_order],NA)[[ifelse(i==1|i==3|Dim[2]==1,1,2)]], names2=list(1:nrow(Cov_cc),NA)[[ifelse(i==1|i==2,1,2)]], digits=1, font=2 )
          #if(i==1 | Dim[1]==1) mtext(side=3, text="Spatial", line=1.5, font=2)
          #if(i==2 | Dim[1]==1) mtext(side=3, text="Spatio-temporal", line=1.5, font=2)
          #if(i==2 | (Dim[2]==1&i==1)) mtext(side=4, text=ifelse(length(Data$ObsModel)==1||Data$ObsModel[2]==0,"Encounter probability","Component #1"), line=0.5, font=2)
          #if(i==4 | (Dim[2]==1&i==3)) mtext(side=4, text=ifelse(length(Data$ObsModel)==1||Data$ObsModel[2]==0,"Positive catch rate","Component #2"), line=0.5, font=2)
        }
        #if( length(Return[[paste0( "Cov_", c("omega1", "epsilon1", "omega2", "epsilon2")[i])]])==0 ){
        #  Return[[paste0( "Cov_", c("omega1", "epsilon1", "omega2", "epsilon2")[i])]] = Cov_cc
        #  if( !is.null(Cov_cc)) Return[[paste0( "Cor_", c("omega1", "epsilon1", "omega2", "epsilon2")[i])]] = cov2cor(Cov_cc)
        #}
      }
    dev.off()

#    # Plot sample
#    ThorsonUtilities::save_fig( file=paste0(plotdir,figname,"--Sample.png"), width=Dim[2]*4+1, height=Dim[1]*4, ... )
#      par(mfrow=Dim, mar=c(0,1,1,0), mgp=mgp, tck=tck, oma=oma)
#      for(i in which(plotTF>0) ){
#        if(i==1) Cov_cc = cov(Report$Omega1_sc)
#        if(i==2) Cov_cc = cov(apply(Report$Epsilon1_sct,MARGIN=2,FUN=as.vector))
#        if(i==3) Cov_cc =
#        if(i==4) Cov_cc = cov(Report$Omega2_sc)
#        if(i==5) Cov_cc = cov(apply(Report$Epsilon2_sct,MARGIN=2,FUN=as.vector))
#        if(i==6)
#        plot_cov( Cov=convert(Cov_cc)[category_order,category_order], names=list(category_names[category_order],NA)[[ifelse(i==1|i==3|Dim[2]==1,1,2)]], names2=list(1:nrow(Cov_cc),NA)[[ifelse(i==1|i==2,1,2)]], digits=1, font=2 )
#        if(i==1 | Dim[1]==1) mtext(side=3, text="Spatial", line=1.5, font=2)
#        if(i==2 | Dim[1]==1) mtext(side=3, text="Spatio-temporal", line=1.5, font=2)
#        if(i==2 | (Dim[2]==1&i==1)) mtext(side=4, text=ifelse(length(Data$ObsModel)==1||Data$ObsModel[2]==0,"Encounter probability","Component #1"), line=0.5, font=2)
#        if(i==4 | (Dim[2]==1&i==3)) mtext(side=4, text=ifelse(length(Data$ObsModel)==1||Data$ObsModel[2]==0,"Positive catch rate","Component #2"), line=0.5, font=2)
#      }
#    dev.off()
  }

  # Return
  return( invisible(Return) )
}

#' Plot covariance matrix
#'
#' \code{plot_cov} plots and formats a covariance matrix
#'
#' @param Cov matrix (covariance or correlation) used for plotting
#' @param zlim numeric-vector (length 2) defining bounds for color-scale of covariance
#' @param names Labels for y-axis
#' @param names2 Labels for x-axis (on top of covariance)
#' @param ncolors Number of colors for color-scale
#' @param digits Number of digits for text-labels of covariance
#' @param ... passed to \code{text} for labelling covariances
plot_cov = function( Cov, zlim=NULL, names=1:nrow(Cov), names2=names, ncolors=21, digits=2, ... ){
  Col = colorRampPalette(colors=c("blue","white","red"))
  if(is.null(zlim) ) zlim = c(-1,1)*max(abs(Cov))
  image(z=Cov[1:nrow(Cov),nrow(Cov):1], x=seq(0,1,length=nrow(Cov)), y=seq(0,1,length=nrow(Cov)), col=Col(ncolors), xaxt="n", yaxt="n", zlim=zlim, yaxt="n", ylab="", xlab="" )
  if(length(names)>1) axis(side=2, at=seq(0,1,length=nrow(Cov)), labels=rev(names), las=1)
  if(length(names2)>1) axis(side=3, at=seq(0,1,length=nrow(Cov)), labels=names2, las=ifelse(all(nchar(names2)==1),1,2) )
  for(i in 1:nrow(Cov)){
  for(j in 1:nrow(Cov)){
    Label = formatC(Cov[i,j],digits=digits,format="f")
    text( y=seq(1,0,length=nrow(Cov))[i], x=seq(0,1,length=nrow(Cov))[j], labels=Label, ...)
  }}
  box()
}

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
#' @inheritParams make_model
#' @param MapDetails_List output from \code{FishStatsUtils::MapDetails_Fn}
#' @param year_set set of parameters to include
#' @param c_set set of categories to include
#' @param ... additional arguments passed to \code{Build_TMB_Fn}
#' @return Tagged list
#' \describe{
#'   \item{Report}{Report output for counter-factual run}
#'   \item{NewBuild_List}{Output from \code{Build_TMB_Fn} using counter-factual parameters}
#' }
Rerun_Fn = function( parhat0, turnoff_pars, loc_x, cov_to_turnoff=1:dim(parhat0[["gamma2_ctp"]])[3], calculate_COG=TRUE, figname=NULL,
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
  #if( !is.null(figname) & !is.null(MapDetails_List) ){
  #  message( "Starting plot for counter-factual run" )
  #  ThorsonUtilities::save_fig( file=figname, width=8, height=12, type="pdf", onefile=TRUE )
  #    plot_density( Report=Report, MapDetails_List=MapDetails_List, type="year", year_set=Year_Set, c_set=c_set )
  #  dev.off()
  #}

  # Return
  return( Return )
}

#' Calculate stability metrics
#'
#' \code{Coherence} calculates buffering (`phi`) and coherence (`psi`)
#'
#' @inheritParams plot_overdispersion
#' @inheritParams make_data
#' @param covhat estimated covariance used for calculating coherence
#' @return Tagged list containing measures of synchrony
#' \describe{
#'   \item{phi_xz}{Synchrony index for each site (x) and each period (row of \code{yearbounds_zz})}
#'   \item{phi_z}{weighted-average of \code{phi_xz} for each period, weighted by average community-abundance at each site in that period}
#'   \item{psi}{Measure of proportion of variance explained by leading eigen-vectors}
#'   \item{L_c}{Cholesky decomposition of \code{covhat}}
#' }
calc_coherence = function( Report, Data, covhat=NULL, yearbounds_zz=matrix(c(1,Data$n_t),nrow=1) ){

  ##################
  # Synchrony
  ##################
  n_z = nrow(yearbounds_zz)

  # Category-specific density ("D"), mean and variance
  varD_xcz = array( NA, dim=c(Data$n_x,Data$n_c,n_z) )

  # Community density ("C") summing across categories
  C_xt = array(0, dim=c(Data$n_x,Data$n_t))
  varC_xz = array(NA, dim=c(Data$n_x,n_z))

  # Area-weighted total biomass ("B") for each category
  B_ct = array(0, dim=c(Data$n_c,Data$n_t))
  B_t = rep(0, Data$n_t)

  # Proportion of total biomass ("P") for each station
  P_xz = array(NA, dim=c(Data$n_x,n_z))

  # Synchrony indices
  phi_xz = array(NA, dim=c(Data$n_x,n_z))
  phi_z = rep(0, n_z)

  # Temporary variables
  temp_xt = array(NA, dim=c(Data$n_x,Data$n_t))
  temp_c = rep(NA, Data$n_c)

  # Derived quantities
  for( tI in 1:Data$n_t ){
    for( cI in 1:Data$n_c ){
      for( xI in 1:Data$n_x ){
        C_xt[xI,tI] = C_xt[xI,tI] + Report$D_xct[xI,cI,tI]
        B_ct[cI,tI] = B_ct[cI,tI] + Report$D_xct[xI,cI,tI] * Data$a_xl[xI,1]
      }
      B_t[tI] = B_t[tI] + B_ct[cI,tI]
    }
  }

  # Loop through periods (only using operations available in TMB, i.e., var, mean, row, col, segment)
  for( zI in 1:n_z){
  for( xI in 1:Data$n_x){
    # Variance for each category
    for( cI in 1:Data$n_c ){
      for( tI in 1:Data$n_t ) temp_xt[xI,tI] = Report$D_xct[xI,cI,tI]
      varD_xcz[xI,cI,zI] = var(temp_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]])
    }
    # Variance for combined across categories
    varC_xz[xI,zI] = var(C_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]])
    # Proportion in each category
    P_xz[xI,zI] = Data$a_xl[xI,1] * mean( C_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]] / B_t[yearbounds_zz[zI,1]:yearbounds_zz[zI,2]] )
    # Synchrony index
    for( cI in 1:Data$n_c ) temp_c[cI] = varD_xcz[xI,cI,zI]
    phi_xz[xI,zI] = varC_xz[xI,zI] / sum(sqrt(temp_c))^2
    phi_z[zI] = phi_z[zI] + (phi_xz[xI,zI] * P_xz[xI,zI])
  }}

  ##################
  # Coherence
  ##################

  # Coherence index
  if( !is.null(covhat) ){
    L_c = eigen(covhat)$values
    psi = 2 * (mean(cumsum(L_c)/sum(L_c))-0.5)
  }else{
    L_c = psi = NULL
  }

  # Return stuff
  Return = list("phi_xz"=phi_xz, "phi_z"=phi_z, "psi"=psi, "P_xz"=P_xz, "L_c"=L_c)
  return( Return )
}


#' Calculate synchrony among species or locations
#'
#' \code{calc_synchrony} calculates synchrony
#'
#' @inheritParams plot_overdispersion
#' @inheritParams make_data
#' @return Tagged list containing measures of synchrony
#' \describe{
#'   \item{phi_xz}{Synchrony index for each site (x) and each period (row of \code{yearbounds_zz})}
#'   \item{phi_z}{weighted-average of \code{phi_xz} for each period, weighted by average community-abundance at each site in that period}
#' }
calc_synchrony = function( Report, Data, yearbounds_zz=matrix(c(1,Data$n_t),nrow=1) ){

  # Index lengths
  n_z = nrow(yearbounds_zz)
  n_t = Data$n_t
  n_c = Data$n_c
  n_x = Data$n_x
  n_y = nrow(Data$t_yz)

  # Extract elements
  pow = function(a,b) a^b
  D_xcy = Report$D_xcy
  a_xl = Data$a_xl

  # Density ("D") or area-expanded total biomass ("B") for each category (use B when summing across sites)
  D_xy = matrix( 0, nrow=n_x, ncol=n_y );
  B_cy = matrix( 0, nrow=n_c, ncol=n_y );
  B_y = rep( 0, n_y );
  # Mean
  meanD_xcz = array( 0, dim=c(n_x,n_c,n_z) );
  meanD_xz = matrix( 0, nrow=n_x, ncol=n_z );
  meanB_cz = matrix( 0, nrow=n_c, ncol=n_z );
  meanB_z = rep( 0, n_z );
  # Sample variance in category-specific density ("D") and biomass ("B")
  varD_xcz = array( 0, dim=c(n_x,n_c,n_z) );
  varD_xz = matrix( 0, nrow=n_x, ncol=n_z );
  varB_cz = matrix( 0, nrow=n_c, ncol=n_z );
  varB_z = rep( 0, n_z );
  varB_xbar_z = rep( 0, n_z );
  varB_cbar_z = rep( 0, n_z );
  maxsdD_xz = matrix( 0, nrow=n_x, ncol=n_z );
  maxsdB_cz = matrix( 0, nrow=n_c, ncol=n_z );
  maxsdB_z = rep( 0, n_z );
  # Proportion of total biomass ("P") for each location or each category
  propB_xz = matrix( 0, nrow=n_x, ncol=n_z );
  propB_cz = matrix( 0, nrow=n_c, ncol=n_z );
  # Synchrony indices
  phi_xz = matrix( 0, nrow=n_x, ncol=n_z );
  phi_cz = matrix( 0, nrow=n_c, ncol=n_z );
  phi_xbar_z = rep( 0, n_z );
  phi_cbar_z = rep( 0, n_z );
  phi_z = rep( 0, n_z );

  # Calculate total biomass for different categories
  for( yI in 1:n_y ){
    for( cI in 1:n_c ){
      for( xI in 1:n_x ){
        D_xy[xI,yI] = D_xy[xI,yI] + D_xcy[xI,cI,yI];
        B_cy[cI,yI] = B_cy[cI,yI] + D_xcy[xI,cI,yI] * a_xl[xI,1];
        B_y[yI] = B_y[yI] + D_xcy[xI,cI,yI] * a_xl[xI,1];
      }
    }
  }
  # Loop through periods (only using operations available in TMB, i.e., var, mean, row, col, segment)
  for( zI in 1:n_z ){
    for( xI in 1:n_x ){
      # Variance for biomass in each category, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
      for( cI in 1:n_c ){
        temp_mean = 0;
        for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) meanD_xcz[xI,cI,zI] = meanD_xcz[xI,cI,zI] + D_xcy[xI,cI,yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
        for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ){
          varD_xcz[xI,cI,zI] = varD_xcz[xI,cI,zI] + pow(D_xcy[xI,cI,yI]-meanD_xcz[xI,cI,zI],2) / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]);
        }
      }
      # Variance for combined biomass across categories, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) meanD_xz[xI,zI] = meanD_xz[xI,zI] + D_xy[xI,yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
        varD_xz[xI,zI] = varD_xz[xI,zI] + pow(D_xy[xI,yI]-meanD_xz[xI,zI],2) / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]);
      }
    }
    for( cI in 1:n_c ){
      # Variance for combined biomass across sites, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) meanB_cz[cI,zI] = meanB_cz[cI,zI] + B_cy[cI,yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
        varB_cz[cI,zI] = varB_cz[cI,zI] + pow(B_cy[cI,yI]-meanB_cz[cI,zI],2) / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]);
      }
    }
    # Variance for combined biomass across sites and categories, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
    for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) meanB_z[zI] = meanB_z[zI] + B_y[yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
    for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
      varB_z[zI] = varB_z[zI] + pow(B_y[yI]-meanB_z[zI],2) / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]);
    }
    # Proportion in each site
    for( xI in 1:n_x ){
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
        propB_xz[xI,zI] = propB_xz[xI,zI] + a_xl[xI,1] * D_xy[xI,yI] / B_y[yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
      }
    }
    # Proportion in each category
    for( cI in 1:n_c ){
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
        propB_cz[cI,zI] = propB_cz[cI,zI] + B_cy[cI,yI] / B_y[yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
      }
    }
    # Species-buffering index (calculate in Density so that areas with zero area are OK)
    for( xI in 1:n_x ){
      for( cI in 1:n_c ){
        maxsdD_xz[xI,zI] = maxsdD_xz[xI,zI] + pow(varD_xcz[xI,cI,zI], 0.5);
      }
      phi_xz[xI,zI] = varD_xz[xI,zI] / pow( maxsdD_xz[xI,zI], 2);
      varB_xbar_z[zI] = varB_xbar_z[zI] + pow(a_xl[xI,1],2) * varD_xz[xI,zI] * propB_xz[xI,zI];
      phi_xbar_z[zI] = phi_xbar_z[zI] + phi_xz[xI,zI] * propB_xz[xI,zI];
    }
    # Spatial-buffering index
    for( cI in 1:n_c ){
      for( xI in 1:n_x ){
        maxsdB_cz[cI,zI] = maxsdB_cz[cI,zI] + a_xl[xI,1] * pow(varD_xcz[xI,cI,zI], 0.5);
      }
      phi_cz[cI,zI] = varB_cz[cI,zI] / pow( maxsdB_cz[cI,zI], 2);
      varB_cbar_z[zI] = varB_cbar_z[zI] + varB_cz[cI,zI] * propB_cz[cI,zI];
      phi_cbar_z[zI] = phi_cbar_z[zI] + phi_cz[cI,zI] * propB_cz[cI,zI];
    }
    # Spatial and species-buffering index
    for( cI in 1:n_c ){
      for( xI in 1:n_x ){
        maxsdB_z[zI] = maxsdB_z[zI] + a_xl[xI,1] * pow(varD_xcz[xI,cI,zI], 0.5);
      }
    }
    phi_z[zI] = varB_z[zI] / pow( maxsdB_z[zI], 2);
  }

  # Return stuff
  Return = list( "phi_z"=phi_z, "phi_xz"=phi_xz, "phi_cz"=phi_cz, "phi_xbar_z"=phi_xz, "phi_cbar_z"=phi_cz, "B_y"=B_y, "B_cy"=B_cy, "D_xy"=D_xy, "D_xcy"=D_xcy)
  return( Return )
}


#' Calculate stability metrics
#'
#' \code{Coherence} calculates buffering (`phi`) and coherence (`psi`)
#'
#' @inheritParams plot_overdispersion
#' @inheritParams make_data
#' @param covhat estimated covariance used for calculating coherence
#' @return Tagged list containing measures of synchrony
#' \describe{
#'   \item{phi_xz}{Synchrony index for each site (x) and each period (row of \code{yearbounds_zz})}
#'   \item{phi_z}{weighted-average of \code{phi_xz} for each period, weighted by average community-abundance at each site in that period}
#'   \item{psi}{Measure of proportion of variance explained by leading eigen-vectors}
#'   \item{L_c}{Cholesky decomposition of \code{covhat}}
#' }
Coherence = function( Report, Data, covhat=NULL, yearbounds_zz=matrix(c(1,Data$n_t),nrow=1) ){

  ##################
  # Synchrony
  ##################
  n_z = nrow(yearbounds_zz)

  # Category-specific density ("D"), mean and variance
  varD_xcz = array( NA, dim=c(Data$n_x,Data$n_c,n_z) )

  # Community density ("C") summing across categories
  C_xt = array(0, dim=c(Data$n_x,Data$n_t))
  varC_xz = array(NA, dim=c(Data$n_x,n_z))

  # Area-weighted total biomass ("B") for each category
  B_ct = array(0, dim=c(Data$n_c,Data$n_t))
  B_t = rep(0, Data$n_t)

  # Proportion of total biomass ("P") for each station
  P_xz = array(NA, dim=c(Data$n_x,n_z))

  # Synchrony indices
  phi_xz = array(NA, dim=c(Data$n_x,n_z))
  phi_z = rep(0, n_z)

  # Temporary variables
  temp_xt = array(NA, dim=c(Data$n_x,Data$n_t))
  temp_c = rep(NA, Data$n_c)

  # Derived quantities
  for( tI in 1:Data$n_t ){
    for( cI in 1:Data$n_c ){
      for( xI in 1:Data$n_x ){
        C_xt[xI,tI] = C_xt[xI,tI] + Report$D_xct[xI,cI,tI]
        B_ct[cI,tI] = B_ct[cI,tI] + Report$D_xct[xI,cI,tI] * Data$a_xl[xI,1]
      }
      B_t[tI] = B_t[tI] + B_ct[cI,tI]
    }
  }

  # Loop through periods (only using operations available in TMB, i.e., var, mean, row, col, segment)
  for( zI in 1:n_z){
  for( xI in 1:Data$n_x){
    # Variance for each category
    for( cI in 1:Data$n_c ){
      for( tI in 1:Data$n_t ) temp_xt[xI,tI] = Report$D_xct[xI,cI,tI]
      varD_xcz[xI,cI,zI] = var(temp_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]])
    }
    # Variance for combined across categories
    varC_xz[xI,zI] = var(C_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]])
    # Proportion in each category
    P_xz[xI,zI] = Data$a_xl[xI,1] * mean( C_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]] / B_t[yearbounds_zz[zI,1]:yearbounds_zz[zI,2]] )
    # Synchrony index
    for( cI in 1:Data$n_c ) temp_c[cI] = varD_xcz[xI,cI,zI]
    phi_xz[xI,zI] = varC_xz[xI,zI] / sum(sqrt(temp_c))^2
    phi_z[zI] = phi_z[zI] + (phi_xz[xI,zI] * P_xz[xI,zI])
  }}

  ##################
  # Coherence
  ##################

  # Coherence index
  if( !is.null(covhat) ){
    L_c = eigen(covhat)$values
    psi = 2 * (mean(cumsum(L_c)/sum(L_c))-0.5)
  }else{
    L_c = psi = NULL
  }

  # Return stuff
  Return = list("phi_xz"=phi_xz, "phi_z"=phi_z, "psi"=psi, "P_xz"=P_xz, "L_c"=L_c)
  return( Return )
}


#' K-fold crossvalidation
#'
#' \code{Crossvalidate_Fn} runs a k-fold crossvalidation analysis on a fitted model object
#'
#' @param record_dir, a directory with writing privilege where the runs are stored
#' @param parhat, a tagged list of parameters from the fitted model, e.g., as generated by \code{parhat <- obj$env$parList()} after covergence
#' @param original_data, a tagged list of data for the VAST model
#' @param group_i, a vector of positive integers, indicating the k-fold group for each observation (default=NULL, which generates a new group_i with even probability)
#' @param kfold, the number of crossvalidation batches used (default=10)
#' @param skip_finished boolean specifying whether to rerun (skip_finished==FALSE) or skip (skip_finished==TRUE) previously completed runs (Default=FALSE)
#' @param skip_finished boolean specifying whether to rerun (skip_finished==FALSE) or skip (skip_finished==TRUE) previously completed runs (Default=FALSE)
#' @param newtonsteps number of extra newton steps to take after optimization (alternative to \code{loopnum})
#' @param ... Additional arguments to pass to \code{VAST::Build_TMB_Fn}
#' @return Results a matrix with total predictive negative log-likelihood for each crossvalidation partition, and number of crossvalidation samples for that partition
Crossvalidate_Fn = function(record_dir, parhat, original_data, group_i=NULL, kfold=10, newtonsteps=1, skip_finished=FALSE, ... ){
  # Lump observations into groups
  if( is.null(group_i) || length(group_i)!=original_data$n_i ){
    message( "Generating group_i" )
    Group_i = sample( x=1:kfold, size=original_data$n_i, replace=TRUE )
  }else{
    message( "Using input group_i" )
    Group_i = group_i
  }
  save( Group_i, file=paste0(record_dir,"Group_i.RData"))

  # Results
  Results = array(NA, dim=c(kfold,3), dimnames=list(paste0("k=",1:kfold),c("predictive_negloglike","number_of_crossvalidation_samples","predictive_nll_per_sample")) )

  # Loop through
  for(i in 1:kfold){
    # Directory
    CrossvalidationDir = paste0(record_dir,"k=",i,"/")
    dir.create( CrossvalidationDir )

    # Run
    if( "Save.RData"%in%list.files(CrossvalidationDir) & skip_finished==TRUE ){
      load(file=paste0(CrossvalidationDir,"Save.RData"))
    }else{
      # Modify data
      Data = original_data
      Data$PredTF_i = ifelse(Group_i==i,1,0)

      # Build new one
      TmbList = VAST::Build_TMB_Fn("TmbData"=Data, "Parameters"=parhat, ...)#, "Random"=NULL)
      #TmbList = VAST::Build_TMB_Fn("TmbData"=Data, "Parameters"=parhat, "RunDir"=TmbDir, "Version"=Version, "loc_x"=loc_x, "RhoConfig"=RhoConfig, "TmbDir"=TmbDir, "Use_REML"=Use_REML) #, "Map"=Save$Map, "Random"=NULL)

      # Extract objects
      Obj = TmbList[["Obj"]]
      TmbList[["Upper"]][grep("logkappa",names(TmbList[["Upper"]]))] = Inf

      # Run model
      Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], newtonsteps=newtonsteps )  # , rel.tol=1e-20

      # Reports
      Report = Obj$report( Obj$env$last.par.best )
      ParHat = Obj$env$parList(Opt$par)

      # Save stuff
      Save = list("Opt"=Opt, "Report"=Report, "ParHat"=ParHat, "TmbData"=Data, "Map"=TmbList$Map)
      save(Save, file=paste0(CrossvalidationDir,"Save.RData"))
      capture.output( Opt, file=paste0(CrossvalidationDir,"Opt.txt"))
    }

    # Load results
    Results[i,c("predictive_negloglike","number_of_crossvalidation_samples")] = c( Save$Report$pred_jnll, sum(Save$TmbData$PredTF_i))
    Results[i,"predictive_nll_per_sample"] = Results[i,"predictive_negloglike"] / Results[i,"number_of_crossvalidation_samples"]
  }

  # Return output
  return( Results)
}


plot_eigen = function( Cov, which2plot=1:min(3,ncol(Cov)), names=1:nrow(Cov), digits=2, las=2, add=FALSE, ... ){
  Eigen = eigen(Cov)
  rownames( Eigen$vectors ) = names
  if(add==FALSE) par( mfrow=par()$mfrow, oma=par()$oma, mar=par()$mar, tck=par()$tck )
  for( cI in which2plot ){
    plot( 1, type="n", xlim=c(0.5,ncol(Cov)+0.5), ylim=c(-1,1.2), xlab="", ylab="", xaxt="n", ... )
    for( rI in 1:nrow(Cov)){
      graphics::lines(y=c(0,Eigen$vectors[rI,cI])*sign(Eigen$vectors[which.max(abs(Eigen$vectors[,cI])),cI]), x=rep(rI,2), lwd=5)
    }
    graphics::abline( h=0, lty="dotted" )
    graphics::legend( "top", bty="n", legend=paste0("Eigenvalue = ",formatC(Eigen$values[cI],format="f",digits=2)) )
    axis(1, at=1:ncol(Cov), labels=names, las=las  )
  }
  return( invisible(Eigen) )
}



Plot_Ellipse_Fn = function( filename, Report, yearset=NULL, years2include=NULL, MappingDetails, MapSizeRatio=c("Height(in)"=4,"Width(in)"=4), Xlim, Ylim, ZinUTM=TRUE, zone=NA, ncol_legend=2, Format="png", ... ){
  # Only run if necessary outputs are available
  if( all(c("mean_Z_tl","cov_Z_tl") %in% names(Report)) ){
    # Decide on years
    if( is.null(years2include) ) years2include = 1:nrow(Report$mean_Z_tl)

    # Plot ellipse -- Lat/Lon OR Lat/Depth
    if(Format=="png") png( file=filename, width=MapSizeRatio['Width(in)'], height=MapSizeRatio['Height(in)'], res=200, units="in")
      # Settings
      ColBorder = colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red","darkred"))(nrow(Report$mean_Z_tl)) # c("blue","red")# c("darkblue","blue","lightblue","lightgreen","yellow","orange","red","darkred")
      if( ZinUTM==FALSE ){
        maps::map(MappingDetails[[1]], MappingDetails[[2]], ylim=mean(Ylim)+1*c(-0.5,0.5)*diff(Ylim), xlim=mean(Xlim)+1*c(-0.5,0.5)*diff(Xlim), fill=TRUE, ...) # , orientation=c(mean(y.lim),mean(x.lim),15)
      }else{
        FishStatsUtils:::Plot_States_in_UTM_Fn( MappingDetails=MappingDetails, fillcol=NA, xlim=Xlim, ylim=Ylim, zone=zone, ... )
      }
      # Plot ellipsoids
      for(t in years2include){
        Shape = mixtools::ellipse(mu=Report$mean_Z_tl[t,1:2], sigma=Report$cov_Z_tl[t,1:2,1:2], alpha=1-(1-pnorm(1))*2, newplot=FALSE, draw=FALSE) #, type="l")
        polygon( Shape, border=ColBorder[t], lwd=3 )    # , col=Col[t]
      }
      if( !is.null(yearset)) legend( "top", fill=ColBorder[years2include], bty="n", legend=yearset[years2include], ncol=ncol_legend)
      #for(t in 2:length(Col)) arrows( x0=mean_Z_tl[t-1,1], y0=mean_Z_tl[t-1,2], x1=mean_Z_tl[t,1], y1=mean_Z_tl[t,2], length=0.1 )
    if(Format=="png") dev.off()
  }else{
    message("mean_Z_tl or cov_Z_tl not found in Report")
  }
}

plotSites <-
function( mapObj, x, y, col, cex=2, ... ){
  # Load package
  require(maptools)

  # Plot map
  if( is.list(mapObj) & length(mapObj)>0 ){
    plot(mapObj, type="l", xlim=mapObj[["bbox"]]['x',], ylim=mapObj[["bbox"]]['y',], ...)
  }else{
    plot(1, type="n", xlim=attr(mapObj,"bbox")['x',], ylim=attr(mapObj,"bbox")['y',], ...)
    if( length(mapObj)>1 ) suppressWarnings( plot(mapObj, add=TRUE) )
  }

  #plot site locations
  points(x=x, y=y, pch=20, cex=cex, col=col)
}

Timeseries_Fn <-
function(Report, FileName, Year_Set, ControlList=list("Width"=4*3, "Height"=2*3, "Res"=200, "Units"='in')){
  # Extract stuff from report
  D_it = Report$D_xt[NN_Extrap$nn.idx,]
  R1_it = Report$R1_xt[NN_Extrap$nn.idx,]
  R2_it = Report$R2_xt[NN_Extrap$nn.idx,]
  D_xt = Report$D_xt
  # Plot time series
  png(file=FileName, width=ControlList$Width, height=ControlList$Height, res=ControlList$Res, units=ControlList$Units)
    par(mfrow=c(2,4), oma=c(0,0,0,0), mar=c(3,3,2,0), mgp=c(1.5,0.5,0), tck=-0.02)
    # Entropy
    EntropyFn = function(Vec){ sum(Vec * log(Vec+1e-250)/log(length(Vec)) ) }
    Entropy_t = apply( Report$D_xt, MARGIN=2, FUN=EntropyFn )
    plot( x=Year_Set, y=Entropy_t, type="l", main="Entropy")
    # Variance
    Var_t = apply( Report$D_xt, MARGIN=2, FUN=var )
    plot( x=Year_Set, y=Var_t, type="l", main="Variance", ylim=c(0,max(Var_t)) )
    # CV
    CV_t = apply( Report$D_xt, MARGIN=2, FUN=function(Vec){ sd(Vec)/mean(Vec) } )
    plot( x=Year_Set, y=CV_t, type="l", main="CV", ylim=c(0,max(CV_t)) )
    # Occupancy probability
    Occup_t = colMeans( R1_it )
    plot( x=Year_Set, y=Occup_t, type="l", main="Occup_t", ylim=c(0,1) )
    # Density
    CondDens_t = colMeans( R2_it )
    plot( x=Year_Set, y=CondDens_t, type="l", main="CondDens_t", ylim=c(0,max(CondDens_t)) )
    # Density
    Index_t = colMeans( D_it )
    plot( x=Year_Set, y=Index_t, type="l", main="Index_t", ylim=c(0,max(Index_t)) )
    # correlation between occupancy and abundance
    Cor_t = sapply(1:length(Year_Set), FUN=function(Num){cor(R1_it[,Num],R2_it[,Num],method="spearman")})
    plot( x=Year_Set, y=Cor_t, type="l", main="Cor_t", ylim=c(-1,1) ); abline(h=0, lty="dotted")
  dev.off()
}


#' @title
#' Plot vessel effects
#'
#' @description
#' \code{Vessel_Fn} plots estimated vessel effects for model components
#'
#' @param TmbData Formatted data inputs, from `Data_Fn(...)`
#' @param Sdreport TMB output from `TMB::sdreport(Obj)`
#' @param FileName_VYplot Full filename (including directory) for center-of-gravity plot
#'
#' @return Return Tagged list of output
#'
Vessel_Fn <-
function( TmbData, Sdreport, FileName_VYplot=NULL ){
  Summary = TMB:::summary.sdreport(Sdreport)
  Return = NULL

  if( !any(c("nu1_vt","nu2_vt") %in% rownames(Summary)) ){
    message( "Not plotting vessel effects because none are present" )
  }else{
    message( "\nPlotting vessel effects..." )
    nu_vt = array(NA, dim=c(TmbData$n_v, TmbData$n_t, 2, 2))
    for(vI in 1:TmbData$n_v){
    for(tI in 1:TmbData$n_t){
      Num = (vI-1)*TmbData$n_t + tI
      nu_vt[vI,tI,1,] = Summary[which(rownames(Summary)=="nu1_vt")[Num],]
      nu_vt[vI,tI,2,] = Summary[which(rownames(Summary)=="nu2_vt")[Num],]
    }}

    # Make plot
    if( !is.null(FileName_VYplot) ) jpeg(FileName_VYplot, width=1.5*TmbData$n_t,height=5,res=200,units="in")
      par(mfrow=c(2,TmbData$n_t), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02, oma=c(0,3,0,0))
      for(eI in 1:2){
      for(tI in 1:TmbData$n_t){
        plot(x=1:TmbData$n_v, y=1:TmbData$n_v, type="n", ylim=range( c(nu_vt[,,eI,1]+nu_vt[,,eI,2],nu_vt[,,eI,1]-nu_vt[,,eI,2]) ), xlab="Vessel", ylab="Effect", main=TmbData$Year_Set[tI])
        if(tI==1) mtext( side=2, outer=FALSE, line=2, text=c("Presence/absence","Positive catch rate")[eI])
        for(vI in 1:TmbData$n_v){
          points( x=vI, y=nu_vt[vI,tI,eI,1])
          lines( x=rep(vI,2), y=c(nu_vt[vI,tI,eI,1]-nu_vt[vI,tI,eI,2],nu_vt[vI,tI,eI,1]+nu_vt[vI,tI,eI,2]))
        }
      }}
    if( !is.null(FileName_VYplot) ) dev.off()

    # Add stuff
    Return = c(Return, list("nu_vt"=nu_vt) )
  }

  # Return stuff
  return( Return )
}


#' Calculate weights for bilinear interpolation
#'
#' \code{bilinear_interp} calculates three coefficients, used to calculate height via blinear interpolation from height and location of three neighbors
#'
#' @param xyz1, 2D location and height of first neighbor
#' @param xyz2, 2D location and height of second neighbor
#' @param xyz3, 2D location and height of third neighbor
#' @param xypred, 2D location of location to interpolate
#' @return Tagged list of useful output
#' \describe{
#'   \item{zpred}{height at location \code{xypred}}
#'   \item{phi}{Coefficients for interpolation}
#' }
bilinear_interp = function( xyz1, xyz2, xyz3, xypred ){
  # Make constaint matrix
  #  1st row:  sum to one
  #  2nd row:  sum to x locations
  #  3rd row:  sum to y location
  B = rbind( rep(1,3), cbind(xyz1[1:2],xyz2[1:2],xyz3[1:2]) )

  # Calculate weighting vector that satisfies these constraints
  phi = solve(B) %*% c(1,xypred)

  # Calculate interpolated value
  zpred = c(xyz1[3],xyz2[3],xyz3[3]) %*% phi
  Return = list("phi"=phi, "zpred"=zpred)
  return( Return )
}
# bilinear_interp( xyz1=c(0,0,0), xyz2=c(1,0,1), xyz3=c(0,1,2), xypred=c(0.25,0.5))


Plot_States_in_UTM_Fn = function( MappingDetails, Rotate=0, fillcol=NA, zone=NA, ... ){
  Map = maps::map(MappingDetails[[1]], MappingDetails[[2]], plot=FALSE, fill=TRUE) # , orientation=c(mean(y.lim),mean(x.lim),15)
  Tmp1 = na.omit( cbind('PID'=cumsum(is.na(Map$x)), 'POS'=1:length(Map$x), 'X'=Map$x, 'Y'=Map$y ))
  # Convert_LL_to_UTM_Fn
  attr(Tmp1,"projection") = "LL"
  attr(Tmp1,"zone") = zone
  tmpUTM = suppressMessages(PBSmapping::convUL(Tmp1))                                                         #$
  coordinates(tmpUTM) = c("X","Y")
  tmp <- maptools::elide( tmpUTM, rotate=Rotate)
  # Plot map
  plot(1, pch="", ... )
  lev = levels(as.factor(tmp@data$PID))
  for(levI in 1:(length(lev))) {
    indx = which(tmpUTM$PID == lev[levI])
    polygon(tmp@coords[indx,'x'], tmp@coords[indx,'y'], col=fillcol)
  }
}


Plot_range_quantiles = function( Data_Extrap, Report, TmbData, a_xl, NN_Extrap, Year_Set=NULL, Prob_vec=c(0.10,0.5,0.90), FileName_Quantiles=paste0(getwd(),"/Range_boundary.png") ){
  # Default inputs
  if( is.null(Year_Set)) Year_Set = 1:TmbData$n_t

  # Plot range boundaries
  calc_quantile = function(x, w, prob=0.5, doplot=FALSE){
    DF = data.frame(x,w)[order(x),]
    DF = cbind(DF, 'cumsum'=cumsum(DF$w)/sum(w))
    Pred = rep(NA,length(prob))
    for(i in 1:length(prob)){
      Between = which( DF$cumsum>prob[i] )[1] - 1:0
      Pred[i] = DF[Between[1],'x'] + (DF[Between[2],'x']-DF[Between[1],'x'])*(prob[i]-DF[Between[1],'cumsum'])/(DF[Between[2],'cumsum']-DF[Between[1],'cumsum'])
    }
    if(doplot==TRUE){
      plot( x=DF$x, y=DF$cumsum, type="l" )
      abline(v=Pred)
    }
    return( Pred )
  }

  # Extrapolation locations
  Data_Extrap_Range = cbind( Data_Extrap[,c('Lat','Lon','N_km','E_km')], 'Include'=ifelse(strip_units(Data_Extrap[,'Area_km2'])>0, TRUE, FALSE) )
  # Add and rescale density
  for(t in 1:length(Year_Set)){
    Data_Extrap_Range = cbind( Data_Extrap_Range, Report$Index_xtl[NN_Extrap$nn.idx,t,1] * (strip_units(Data_Extrap[,'Area_km2']) / a_xl[NN_Extrap$nn.idx,1]) )
    colnames( Data_Extrap_Range )[ncol(Data_Extrap_Range)] = paste0("Year_",Year_Set[t])
    Data_Extrap_Range[,paste0("Year_",Year_Set[t])] = ifelse( strip_units(Data_Extrap[,'Area_km2'])==0 & a_xl[NN_Extrap$nn.idx,1]==0, 0, Data_Extrap_Range[,paste0("Year_",Year_Set[t])])
  }

  # Plot
  png( file=FileName_Quantiles, width=6.5, height=6.5, res=200, units="in")
    par( mfcol=c(2,2), mar=c(0,2,2,0), mgp=c(1.75,0.25,0), tck=-0.02, oma=c(4,0,0,0))
    for(z in 1:4){
      Range_Quantile = array( NA, dim=c(length(Year_Set),3), dimnames=list(NULL,c("min","mid","max")) )
      for(t in 1:nrow(Range_Quantile)) Range_Quantile[t,] = calc_quantile( x=Data_Extrap_Range[,c('Lat','Lon','N_km','E_km')[z]], w=Data_Extrap_Range[,paste0('Year_',Year_Set[t])], prob=c(0.05,0.5,0.95), doplot=FALSE )
      matplot( y=Range_Quantile, x=Year_Set, type="l", col="black", lty="solid", lwd=2, xlab="", ylab="", main=c('Latitude','Longitude','Nothings','Eastings')[z], xaxt="n" )
      if(z %in% c(2,4)) axis(1)
    }
    mtext( side=1, text="Year", outer=TRUE, line=2)
  dev.off()

  # Return
  Return = list( "Data_Extrap_Range"=Data_Extrap_Range )
  return( invisible(Return) )
}


#' Convert from Lat-Long to Eastings-Northings using WGS
#'
#' \code{Convert_LL_to_EastNorth_Fn} converts from Latitude-Longitude to World Geodetic System Eastings-Northings for a given location
#'
#' @param Lat vector of latitudes
#' @param Lon vector of longitudes
#' @param crs EPSG reference for coordinate reference system (CRS) defining Eastings-Northings after transformation
#' @return A data frame with the following columns
#' \describe{
#'   \item{E_km}{The eastings for each value of Lon (in kilometers)}
#'   \item{N_km}{The northings for each value of Lat (in kilometers)}
#' }
Convert_LL_to_EastNorth_Fn <-
function( Lon, Lat, crs=NA ){
  # SEE:  https://github.com/nwfsc-assess/geostatistical_delta-GLMM/issues/25#issuecomment-345825230

  # Attach package
  #require(rgdal)
  #on.exit( detach("package:rgdal") )
  require(sp)
  on.exit( detach("package:sp") )

  # Transform
  dstart<-data.frame(lon=Lon, lat=Lat) # that's the object
  coordinates(dstart) <- c("lon", "lat")
  proj4string(dstart) <- CRS("+init=epsg:4326") # that's the lat long projection
  CRS.new <- CRS(crs) # that's the eastings and northings projection
  dstart.t <- spTransform(dstart, CRS.new) # here's where you transform

  # Clean up
  dstart.t = cbind( "E_km"=dstart.t@coords[,"lon"]/1000, "N_km"=dstart.t@coords[,'lat']/1000 )
  attr(dstart.t,"zone") = crs

  # Return results
  return( dstart.t )
}


#' Convert from Lat-Long to UTM
#'
#' \code{Convert_LL_to_UTM_Fn} converts from Latitude-Longitude to Universal Transverse Mercator projections for a given location
#'
#' @param Lat vector of latitudes
#' @param Lon vector of longitudes
#' @param zone UTM zone (integer between 1 and 60) or alphanumeric CRS code used by package rgdal to convert latitude-longitude coordinates to projection in kilometers; \code{zone=NA} uses UTM and automatically detects the appropriate zone
#' @param flip_around_dateline boolean specifying whether to flip Lat-Lon locations around the dateline, and then retransform back (only useful if Lat-Lon straddle the dateline)
#' @return A data frame with the following columns
#' \describe{
#'   \item{X}{The UTM eastings for each value of Lon}
#'   \item{Y}{The UTM northings measured from the equator for each Lat}
#' }
Convert_LL_to_UTM_Fn <-
function( Lon, Lat, zone=NA, flip_around_dateline=FALSE ){

  # Convert
  # if zone=NA or NULL, then it automatically detects appropriate zone
  Tmp = cbind('PID'=1,'POS'=1:length(Lon),'X'=Lon,'Y'=Lat)
  if( flip_around_dateline==TRUE ) Tmp[,'X'] = ifelse( Tmp[,'X']>0, Tmp[,'X']-180, Tmp[,'X']+180)
  attr(Tmp,"projection") = "LL"
  attr(Tmp,"zone") = zone
  tmpUTM = PBSmapping::convUL(Tmp)                                                         #$
  if( !is.na(zone)) message("convUL: For the UTM conversion, used zone ",zone," as specified")

  # Return results
  return( tmpUTM )
}

#' @title
#' Plot maps with areal results
#'
#' @description
#' \code{PlotMap_Fn} is a hidden function to plot a map and fill in regions with colors to represent intensity in an areal-interpretion of model results
#'
#' @inheritParams plot_maps
#' @param plot_legend_fig Boolean, whether to plot a separate figure for the heatmap legend or not
#' @param land_color color for filling in land (use \code{land_color=rgb(0,0,0,alpha=0)} for transparent land)
#' @param ... arguments passed to \code{par}
#'
#' @details
#' This function was necessary to build because \code{mapproj::mapproject} as used in \code{maps::map} has difficulties with both rotations (for most projections) and
#' truncating the cocuntry boundaries within the plotting region (which \code{mapproj::mapproject} appears to do prior to projection,
#' so that the post-projection is often missing boundaries that are within the plotting rectangle).  I use rectangular projections by default, but Lamberts or Albers conformal
#' projections would also be useful for many cases.
PlotMap_Fn <-
function(MappingDetails, Mat, PlotDF, MapSizeRatio=c('Width(in)'=4,'Height(in)'=4), Xlim, Ylim, FileName=paste0(getwd(),"/"), Year_Set,
         Rescale=FALSE, Rotate=0, Format="png", Res=200, zone=NA, Cex=0.01, textmargin="", add=FALSE, pch=15,
         outermargintext=c("Eastings","Northings"), zlim=NULL, Col=NULL,
         Legend=list("use"=FALSE, "x"=c(10,30), "y"=c(10,30)), mfrow=c(1,1), plot_legend_fig=TRUE, land_color="grey", ignore.na=FALSE,
         map_style="rescale", ...){

  # Warning
  warning( "`PlotMap_Fn` is soft-deprecated, please use `plot_variable` instead for many improvements")

  # Check for problems
  if( length(Year_Set) != ncol(Mat) ){
    warning( "Year_Set and `ncol(Mat)` don't match: Changing Year_Set'")
    Year_Set = 1:ncol(Mat)
  }

  # Transform to grid or other coordinates
  Mat = Mat[PlotDF[,'x2i'],,drop=FALSE]
  Which = which( PlotDF[,'Include']>0 )
  if( Rescale!=FALSE ) Mat = Mat / outer(rep(Rescale,nrow(Mat)), colMeans(Mat[Which,]))

  # Plotting functions
  f = function(Num, zlim=NULL){
    if( is.null(zlim)) Return = ((Num)-min(Num,na.rm=TRUE))/max(diff(range(Num,na.rm=TRUE)),0.01)
    if( !is.null(zlim)) Return = ((Num)-zlim[1])/max(diff(zlim),0.01)
    return( Return )
  }
  if( is.null(Col)) Col = colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))
  if( is.function(Col)) Col = Col(50)

  # Plot
  Par = list( mfrow=mfrow, ...)
  if(Format=="png"){
    png(file=paste0(FileName, ".png"),
        width=Par$mfrow[2]*MapSizeRatio['Width(in)'],
        height=Par$mfrow[1]*MapSizeRatio['Height(in)'], res=Res, units='in')
  }
  if(Format=="jpg"){
    jpeg(file=paste0(FileName, ".jpg"),
         width=Par$mfrow[2]*MapSizeRatio['Width(in)'],
         height=Par$mfrow[1]*MapSizeRatio['Height(in)'], res=Res, units='in')
  }
  if(Format%in%c("tif","tiff")){
    tiff(file=paste0(FileName, ".tif"),
         width=Par$mfrow[2]*MapSizeRatio['Width(in)'],
         height=Par$mfrow[1]*MapSizeRatio['Height(in)'], res=Res, units='in')
  }
    if(add==FALSE) par( Par )          # consider changing to Par=list() input, which overloads defaults a la optim() "control" input
    for(tI in 1:length(Year_Set)){
      if( is.null(MappingDetails) ){
        plot(1, type="n", ylim=Ylim, xlim=Xlim, main="", xlab="", ylab="", mar=par()$mar )#, main=Year_Set[t])
        points(x=PlotDF[Which,'Lon'], y=PlotDF[Which,'Lat'], col=Col[ceiling(f(Mat[Which,],zlim=zlim)[,t]*(length(Col)-1))+1], cex=0.01)
      }else{
        # If not rotating:  Use simple plot
        if( Rotate==0 ){
          Col_Bin = ceiling( f(Mat[Which,,drop=FALSE],zlim=zlim)[,tI]*(length(Col)-1) ) + 1
          if( map_style=="rescale" ){
            # Make plot size using plot(.)
            plot( 1, type="n", ylim=mean(Ylim)+c(-0.5,0.5)*diff(Ylim), xlim=mean(Xlim)+c(-0.5,0.5)*diff(Xlim), xaxt="n", yaxt="n", xlab="", ylab="" )
            points(x=PlotDF[Which,'Lon'], y=PlotDF[Which,'Lat'], col=Col[Col_Bin], cex=Cex, pch=pch)
            Map = maps::map(MappingDetails[[1]], MappingDetails[[2]], plot=FALSE)
            map( Map, add=TRUE ) #, col=land_color, fill=TRUE )  # ->  Using land-fill color produces weird plotting artefacts
          }else{
            # Make plot size using map(.)
            map(MappingDetails[[1]], MappingDetails[[2]], ylim=mean(Ylim)+c(-0.5,0.5)*diff(Ylim), xlim=mean(Xlim)+c(-0.5,0.5)*diff(Xlim), plot=TRUE, fill=FALSE, mar=c(0.1,0.1,par("mar")[3],0), myborder=0 )
            points(x=PlotDF[Which,'Lon'], y=PlotDF[Which,'Lat'], col=Col[Col_Bin], cex=Cex, pch=pch)
            map(MappingDetails[[1]], MappingDetails[[2]], plot=TRUE, col=land_color, fill=TRUE, add=TRUE )
          }
        }
        # If rotating:  Record all polygons; Rotate them and all points;  Plot rotated polygons;  Plot points
        if( Rotate!=0 ){
          # Extract map features
          boundary_around_limits = 3
          Map = maps::map(MappingDetails[[1]], MappingDetails[[2]], plot=FALSE, ylim=mean(Ylim)+boundary_around_limits*c(-0.5,0.5)*diff(Ylim), xlim=mean(Xlim)+boundary_around_limits*c(-0.5,0.5)*diff(Xlim), fill=TRUE) # , orientation=c(mean(y.lim),mean(x.lim),15)
          Tmp1 = na.omit( cbind('PID'=cumsum(is.na(Map$x)), 'POS'=1:length(Map$x), 'X'=Map$x, 'Y'=Map$y, matrix(0,ncol=length(Year_Set),nrow=length(Map$x),dimnames=list(NULL,Year_Set))) )
          TmpLL = rbind( Tmp1, cbind('PID'=max(Tmp1[,1])+1, 'POS'=1:length(Which)+max(Tmp1[,2]), 'X'=PlotDF[Which,'Lon'], 'Y'=PlotDF[Which,'Lat'], Mat[Which,]) )
          tmpUTM = TmpLL
          # Convert map to Eastings-Northings
          if( is.numeric(zone) ){
            tmpUTM[,c('X','Y')] = as.matrix(Convert_LL_to_UTM_Fn( Lon=TmpLL[,'X'], Lat=TmpLL[,'Y'], zone=zone, flip_around_dateline=ifelse(MappingDetails[[1]]%in%c("world2","world2Hires"),FALSE,FALSE) )[,c('X','Y')])
          }else{
            tmpUTM[,c('X','Y')] = as.matrix(Convert_LL_to_EastNorth_Fn( Lon=TmpLL[,'X'], Lat=TmpLL[,'Y'], crs=zone )[,c('E_km','N_km')])
          }
          # Rotate map features and maps simultaneously
          tmpUTM = data.frame(tmpUTM)
          sp::coordinates(tmpUTM) = c("X","Y")
          tmpUTM_rotated <- maptools::elide( tmpUTM, rotate=Rotate)
          plot( 1, type="n", xlim=range(tmpUTM_rotated@coords[-c(1:nrow(Tmp1)),'x']), ylim=range(tmpUTM_rotated@coords[-c(1:nrow(Tmp1)),'y']), xaxt="n", yaxt="n" )
          Col_Bin = ceiling( f(tmpUTM_rotated@data[-c(1:nrow(Tmp1)),-c(1:2),drop=FALSE],zlim=zlim)[,tI]*(length(Col)-1) ) + 1
          if( ignore.na==FALSE && any(Col_Bin<1 | Col_Bin>length(Col)) ) stop("zlim doesn't span the range of the variable")
          points(x=tmpUTM_rotated@coords[-c(1:nrow(Tmp1)),'x'], y=tmpUTM_rotated@coords[-c(1:nrow(Tmp1)),'y'], col=Col[Col_Bin], cex=Cex, pch=pch)
          # Plot map features
          lev = levels(as.factor(tmpUTM_rotated@data$PID))
          for(levI in 1:(length(lev)-1)) {
            indx = which(tmpUTM$PID == lev[levI])
            if( var(sign(TmpLL[indx,'Y']))==0 ){
              polygon(x=tmpUTM_rotated@coords[indx,'x'], y=tmpUTM_rotated@coords[indx,'y'], col=land_color)
            }else{
              warning( "Skipping map polygons that straddle equator, because PBSmapping::convUL doesn't work for these cases" )
            }
          }
        }
      }
      title( Year_Set[tI], line=0.1, cex.main=ifelse(is.null(Par$cex.main), 1.8, Par$cex.main), cex=ifelse(is.null(Par$cex.main), 1.8, Par$cex.main) )
      box()
    }
    # Include legend
    if( Legend$use==TRUE ){
      FishStatsUtils:::smallPlot( FishStatsUtils:::Heatmap_Legend(colvec=Col, heatrange=list(range(Mat[Which,],na.rm=TRUE),zlim)[[ifelse(is.null(zlim),1,2)]], dopar=FALSE), x=Legend$x, y=Legend$y, mar=c(0,0,0,0), mgp=c(2,0.5,0), tck=-0.2, font=2 )  #
    }
    # Margin text
    if(add==FALSE) mtext(side=1, outer=TRUE, outermargintext[1], cex=1.75, line=par()$oma[1]/2)
    if(add==FALSE) mtext(side=2, outer=TRUE, outermargintext[2], cex=1.75, line=par()$oma[2]/2)
  if(Format %in% c("png","jpg","tif","tiff")) dev.off()
  # Legend
  if( plot_legend_fig==TRUE ){
    if(Format=="png"){
      png(file=paste0(FileName, "_Legend.png",sep=""),
          width=1, height=2*MapSizeRatio['Height(in)'], res=Res, units='in')
    }
    if(Format=="jpg"){
      jpeg(file=paste0(FileName, "_Legend.jpg",sep=""),
           width=1, height=2*MapSizeRatio['Height(in)'], res=Res, units='in')
    }
    if(Format%in%c("tif","tiff")){
      tiff(file=paste0(FileName, "_Legend.tif",sep=""),
           width=1, height=2*MapSizeRatio['Height(in)'], res=Res, units='in')
    }
    if(Format %in% c("png","jpg","tif","tiff")){
      FishStatsUtils:::Heatmap_Legend( colvec=Col, heatrange=list(range(Mat,na.rm=TRUE),zlim)[[ifelse(is.null(zlim),1,2)]], textmargin=textmargin )
      dev.off()
    }
  }
  return( invisible(list("Par"=Par)) )
}

Heatmap_Legend <-
function( colvec, heatrange, textmargin=NULL, labeltransform="uniform", dopar=TRUE, side=4 ){
  if(dopar==TRUE) par( xaxs="i", yaxs="i", mar=c(1,0,1,2+ifelse(is.null(textmargin),0,1.5)), mgp=c(1.5,0.25,0), tck=-0.02 )
  N = length(colvec)
  Y = seq(heatrange[1], heatrange[2], length=N+1)
  plot( 1, type="n", xlim=c(0,1), ylim=heatrange, xlab="", ylab="", main="", xaxt="n", yaxt="n", cex.main=1.5, xaxs="i", yaxs="i")
  for( i in 1:N) polygon( x=c(0,1,1,0), y=Y[c(i,i,i+1,i+1)], col=colvec[i], border=NA)
  if( labeltransform=="uniform" ) Labels = pretty(heatrange)
  if( labeltransform=="inv_log10" ) Labels = 10^pretty(heatrange)
  axis(side=4, at=pretty(heatrange), labels=Labels )
  if(!is.null(textmargin)) mtext(side=side, text=textmargin, line=2, cex=1.5, las=0)
  #mtext(side=1, outer=TRUE, line=1, "Legend")
}


#' Inset small plot within figure
#'
#' Inset plot with margins, background and border (based on: https://github.com/cran/berryFunctions/blob/master/R/smallPlot.R)
#'
#' @return parameters of small plot, invisible.
smallPlot <- function( expr, x=c(5,70), y=c(50,100), x1,y1,x2,y2, mar=c(12, 14, 3, 3), mgp=c(1.8, 0.8, 0),
  bg=par("bg"), border=par("fg"), las=1, resetfocus=TRUE, ...){

  # Input check:                               #  y1 | P1       |
  if(missing(x1)) x1 <- min(x, na.rm=TRUE)     #     |          |
  if(missing(x2)) x2 <- max(x, na.rm=TRUE)     #  y2 |       P2 |
  if(missing(y1)) y1 <- max(y, na.rm=TRUE)     #     ------------
  if(missing(y2)) y2 <- min(y, na.rm=TRUE)     #       x1    x2

  # catch outside plot:
  if(x1<0)  {x1 <- 0;   warning("x (",x1,") set to 0.")}
  if(y2<0)  {y2 <- 0;   warning("y (",y2,") set to 0.")}
  if(x2>100){x2 <- 100; warning("x (",x2,") set to 100.")}
  if(y1>100){y1 <- 100; warning("y (",y1,") set to 100.")}

  # control for 0:1 input:
  if(diff(range(x, na.rm=TRUE)) < 1  |  diff(range(y, na.rm=TRUE)) < 1  ){
    stop("x or y was probably given as coodinates between 0 and 1. They must be between 0 and 100.")
  }

  # old parameters to be restored at exit:
  op <- par(no.readonly=TRUE)

  # inset plot: background, border
  par(plt=c(x1, x2, y2, y1)/100, new=TRUE, mgp=mgp) # plt / fig
  plot.new() # code line from ade4::add.scatter
  u <- par("usr")
  rect(u[1], u[3], u[2], u[4], col=bg, border=border)

  # inset plot: margins
  par(plt=c(x1+mar[2], x2-mar[4], y2+mar[1], y1-mar[3])/100, new=TRUE, las=las, ...)

  # Actual plot:
  expr

  # par of small plot:
  sp <- par(no.readonly=TRUE)

  # par reset
  if(resetfocus){
    if( par("mfrow")[1]==1 & par("mfrow")[2]==1  ){
      par(op) # ruins multiple figure plots, so:
    }else{
      par(plt=op$plt, new=op$new, mgp=op$mgp, las=op$las)
    }
  }
  return(invisible(sp))
}


#' Format habitat covariate matrix
#'
#' \code{format_covariates} formats an array used by \code{Data_Fn} for habitat (a.k.a. density) covariates
#'
#' @param Lat_e, Latitude for covariate sample e
#' @param Lon_e, Longitude for covariate sample e
#' @param t_e, Time (e.g., year) for covariate sample e
#' @param Cov_ep, matrix of covariates
#' @param Spatial_List, Output from \code{FishStatsUtils::Spatial_Information_Fn}, representing spatial location of knots
#' @param FUN, function used to aggregate observations associated with each knot-time combination
#' @param Year_Set, Set of times \code{t_e} used when generating \code{Cov_xtp}
#' @param na.omit, What to do when some knot-time combination has no observation. Options include \code{"error"} which throw an error, or \code{"time-average"} which fills in the average for other years with observations for each knot
#' @return Tagged list of useful output
#' \describe{
#'   \item{Cov_xtp}{3-dimensional array for use in \code{Data_Fn}}
#'   \item{var_p}{a matrix summarizing the data-level variance, variance among knots, and lost variance when aggregating from data to knots for each covariate}
#' }
format_covariates = function( Lat_e, Lon_e, t_e, Cov_ep, Extrapolation_List, Spatial_List, FUN=mean, Year_Set=min(t_e):max(t_e), na.omit="error" ){

  # Knots in UTM: Spatial_List$loc_x
  # Info for projection:  Extrapolation_List[c('zone','flip_around_dateline')]

  # Backwards compatibility
  if( is.null(Spatial_List$fine_scale) ) Spatial_List$fine_scale = FALSE

  # Determine grids to use
  if( Spatial_List$fine_scale==FALSE ){
    loc_g = Spatial_List$loc_x
  }else{
    if(is.null(Spatial_List)) stop("Problem with `Spatial_List`")
    loc_g = Spatial_List$loc_g
  }

  #
  if( is.vector(Cov_ep)) Cov_ep = matrix(Cov_ep,ncol=1)

  # Step 1:  Project Lat_e and Lon_e to same coordinate system as knots
  if( is.numeric(Extrapolation_List$zone) ){
    loc_e = Convert_LL_to_UTM_Fn( Lon=Lon_e, Lat=Lat_e, zone=Extrapolation_List$zone, flip_around_dateline=Extrapolation_List$flip_around_dateline )                                                         #$
    loc_e = cbind( 'E_km'=loc_e[,'X'], 'N_km'=loc_e[,'Y'])
  }else{
    loc_e = Convert_LL_to_EastNorth_Fn( Lon=Lon_e, Lat=Lat_e, crs=Extrapolation_List$zone )
  }

  # Associate each knot with all covariate measurements that are closest to that knot
  if( na.omit %in% c("error","time-average") ){
    # Step 2: Determine nearest knot for each LatLon_e
    NN = RANN::nn2( data=loc_g[,c('E_km','N_km')], query=loc_e[,c('E_km','N_km')], k=1 )$nn.idx[,1]

    # Step 3: Determine average covariate for each knot and year
    Cov_xtp = NULL
    for(pI in 1:ncol(Cov_ep)){
      Cov_xt = tapply( Cov_ep[,pI], INDEX=list(factor(NN,levels=1:nrow(loc_g)), factor(t_e,levels=Year_Set)), FUN=FUN )
      Cov_xtp = abind::abind( Cov_xtp, Cov_xt, along=3 )
    }

    # Step 4: QAQC
    if( any(is.na(Cov_xtp)) ){
      if( na.omit=="error" ){
        stop("No covariate value for some combination of knot and year")
      }
      if( na.omit=="time-average"){
        Cov_xtp = ifelse( is.na(Cov_xtp), aperm(apply(Cov_xtp,MARGIN=c(1,3),FUN=mean, na.rm=TRUE)%o%rep(1,length(Year_Set)),c(1,3,2)), Cov_xtp )
      }
    }

    # Step 5: Calculate variance lost when aggregating to knots
    KnotVar_p = apply( Cov_xtp, MARGIN=3, FUN=function(vec){var(as.vector(vec))} )
    DataVar_p = apply( Cov_ep, MARGIN=2, FUN=var )
    LostVar_p = 1 - KnotVar_p / DataVar_p
    var_p = matrix( c(KnotVar_p,DataVar_p,LostVar_p), nrow=3, byrow=TRUE)
    dimnames(var_p) = list( c("Variance_among_knots","Variance_among_observations","Proportion_of_variance_lost"), colnames(Cov_ep) )
    for(pI in 1:ncol(Cov_ep)){
      message( "Projecting to knots lost ", formatC(100*LostVar_p[pI],format="f",digits=1), "% of variance for variable ", ifelse(is.null(colnames(Cov_ep)[pI]),pI,colnames(Cov_ep)[pI]))
    }
  }

  # Associate each knot with the X closest measurements in a given year
  if( is.numeric(na.omit) ){
    Cov_xtp = array(NA, dim=c(nrow(loc_g),length(Year_Set),ncol(Cov_ep)), dimnames=list(NULL,Year_Set,colnames(Cov_ep)) )

    for(tI in 1:length(Year_Set)){
      locprime_e = loc_e[which(t_e==Year_Set[tI]),c('E_km','N_km'),drop=FALSE]
      if( nrow(locprime_e) == 0 ) stop("No measurements for year ", Year_Set[tI] )
      NN = RANN::nn2( data=locprime_e[,c('E_km','N_km')], query=loc_g[,c('E_km','N_km')], k=na.omit )$nn.idx
      for(xI in 1:dim(Cov_xtp)[1] ){
      for(pI in 1:dim(Cov_xtp)[3] ){
        Cov_xtp[xI,tI,pI] = FUN( Cov_ep[which(t_e==Year_Set[tI]),pI][NN[xI,]] )
      }}
    }
  }

  # Calculate interpolated value
  Return = list( Cov_xtp=Cov_xtp )
  if( na.omit %in% c("error","time-average") ){
    Return[["var_p"]] = var_p
  }
  return( Return )
}


#' Simulate data
#'
#' \code{Geostat_Sim} simulates data for use when testing \code{VAST}
#'
#' @param Sim_Settings an optional tagged list for specifying input parameters (see function code to determine settings)
#' @inheritParams make_spatial_info
#' @param Data_Geostat A data frame with column headers \code{c('Lon','Lat','Year','Vessel','AreaSwept_km2')} containing sample design to mimic
#' @param standardize_fields Boolean, whether to ensure that random fields have sample mean and standard deviation equal to their inputted values
#' @return Return Tagged list of output
#' \describe{
#'   \item{Data_Geostat}{Simulated data for analysis}
#'   \item{B_tl}{True biomass for each year and stratum}
#' }
Geostat_Sim <-
function(Sim_Settings, Extrapolation_List, Data_Geostat=NULL, MakePlot=FALSE, DateFile=paste0(getwd(),"/"), standardize_fields=FALSE ){
  # Terminology
  # O: Spatial component; E: spatiotemporal component
  # O1: presence-absence; O2: positive catch rate
  # P1 : linear predictor of presence/absence in link-space
  # R1: predictor of presence/absence in natural-space (transformed out of link-space)
  # v_i: vessel for sample i
  # t_i: year for sample i

  # Warning
  warning( "`Geostat_Sim` is deprecated, please use `obj$env$simulate(.)` instead for supported simulation method")

  # Specify default values
  Settings = list("beta1_mean"=0, "beta2_mean"=0, "beta1_slope"=0, "beta2_slope"=0, "beta1_sd"=0, "beta2_sd"=0, "Nyears"=10, "Nsamp_per_year"=600,
    "Depth1_km"=0, "Depth1_km2"=0, "Dist1_sqrtkm"=0, "Depth2_km"=0, "Depth2_km2"=0, "Dist2_sqrtkm"=0,
    "SigmaO1"=0.5, "SigmaO2"=0.5, "SigmaE1"=0.1, "SigmaE2"=0.1, "SigmaV1"=0, "SigmaV2"=0, "SigmaVY1"=0, "SigmaVY2"=0,
    "Range1"=1000, "Range2"=500, "SigmaM"=1, "ObsModel"=c(2,0),
    "Nages"=1, "M"=Inf, "K"=Inf, "Linf"=1, "W_alpha"=1, "W_beta"=3, "Selex_A50_mean"=0, "Selex_A50_sd"=0, "Selex_Aslope"=Inf )

  # Replace defaults with provided values (if any)
  for( i in seq_len(length(Sim_Settings)) ){
    if(names(Sim_Settings)[i] %in% names(Settings)){
      Settings[[match(names(Sim_Settings)[i],names(Settings))]] = Sim_Settings[[i]]
    }
  }

  start_time = Sys.time()

  # Local functions
  RFsim = function(model, x, y, standardize=TRUE){
    if( model_O1@par.general$var>0 ){
      vec = RandomFields::RFsimulate(model=model, x=x, y=y)@data[,1]
      if(standardize==TRUE) vec = (vec - mean(vec)) / sd(vec) * sqrt(model@par.general$var)
    }else{
      vec = rep(0,length(x))
    }
    return(vec)
  }
  rPoisGam = function( n=1, shape, scale, intensity ){
    Num = rpois( n=n, lambda=intensity )
    Biomass = rep(0, length=n)
    for(i in which(Num>0) ) Biomass[i] = sum( rgamma(n=Num[i], shape=shape, scale=scale) )
    return( Biomass )
  }

  # Initialize random-field objects
  model_O1 = RandomFields::RMgauss(var=Settings[['SigmaO1']]^2, scale=Settings[['Range1']])
  model_E1 = RandomFields::RMgauss(var=Settings[['SigmaE1']]^2, scale=Settings[['Range1']])
  model_O2 = RandomFields::RMgauss(var=Settings[['SigmaO2']]^2, scale=Settings[['Range2']])
  model_E2 = RandomFields::RMgauss(var=Settings[['SigmaE2']]^2, scale=Settings[['Range2']])

  # Define sampling locations
  if( is.null(Data_Geostat) ){
    s_i = sample(1:nrow(Extrapolation_List$Data_Extrap), size=Settings[['Nyears']]*Settings[['Nsamp_per_year']]) #, nrow=Settings[['Nsamp_per_year']], ncol=Settings[['Nyears']])
    t_i = rep( 1:Settings[['Nyears']], each=Settings[['Nsamp_per_year']])
    v_i = rep(rep( 1:4, each=Settings[['Nsamp_per_year']]/4), Settings[['Nyears']])
    w_i = rep(0.01, length(s_i))
  }else{
    NN_domain = RANN::nn2( data=Extrapolation_List$Data_Extrap[,c('Lon','Lat')], query=Data_Geostat[,c('Lon','Lat')], k=1)
    s_i = NN_domain$nn.idx[,1]
    t_i = Data_Geostat[,'Year'] - min(Data_Geostat[,'Year']) + 1
    v_i = as.numeric(factor(Data_Geostat[,'Vessel']))
    w_i = Data_Geostat[,'AreaSwept_km2']
  }

  # Duplicate locations if Nages>1
  if( Settings[["Nages"]]>=2 ){
    a_i = rep( 1:Settings[["Nages"]], times=length(s_i) )
    s_i = rep( s_i, each=Settings[["Nages"]] )
    t_i = rep( t_i, each=Settings[["Nages"]] )
    v_i = rep( v_i, each=Settings[["Nages"]] )
    w_i = rep( w_i, each=Settings[["Nages"]] )
  }else{
    a_i = rep( 1, times=length(s_i) )
  }

  # Generate locations, years, vessel for each sample i
  loc_s = Extrapolation_List$Data_Extrap[,c('E_km','N_km','Lat','Lon')]
  loc_i = loc_s[s_i,]

  # Generate intercepts
  exp_beta1_t = Settings[['beta1_mean']] + Settings[['beta1_slope']] * (1:max(t_i)-mean(1:max(t_i))/2)
  beta1_t = rnorm( max(t_i), mean=exp_beta1_t, sd=Settings[['beta1_sd']] )
  exp_beta2_t = Settings[['beta2_mean']] + Settings[['beta2_slope']] * (1:max(t_i)-mean(1:max(t_i))/2)
  beta2_t = rnorm( max(t_i), mean=exp_beta2_t, sd=Settings[['beta2_sd']] )

  # Simulate spatial and spatio-temporal components
  O1_s = RFsim(model=model_O1, x=loc_s[,'E_km'], y=loc_s[,'N_km'], standardize=standardize_fields)
  O2_s = RFsim(model=model_O2, x=loc_s[,'E_km'], y=loc_s[,'N_km'], standardize=standardize_fields)
  message( "Finished variation that is constant over time" )
  E1_st = E2_st = matrix(NA, nrow=nrow(loc_s), ncol=max(t_i) )
  for(t in 1:max(t_i)){
    E1_st[,t] = RFsim(model=model_E1, x=loc_s[,'E_km'], y=loc_s[,'N_km'], standardize=standardize_fields)
    E2_st[,t] = RFsim(model=model_E2, x=loc_s[,'E_km'], y=loc_s[,'N_km'], standardize=standardize_fields)
    message( "Finished variation for year ", t )
  }

  # Extract covariates
  X_sj = array( 0, dim=c(nrow(Extrapolation_List$Data_Extrap),3), dimnames=list(NULL,c('Depth_km','Depth_km2','Rock_dist_')) )
  for( colI in 1:ncol(X_sj)){
    if(colnames(X_sj)[colI] %in% colnames(Extrapolation_List$Data_Extrap)) X_sj[,colI] = as.matrix(Extrapolation_List$Data_Extrap[,c('Depth_km','Depth_km2','Rock_dist_')[colI]])
    X_sj = ifelse( is.na(X_sj), outer(rep(1,nrow(X_sj)),colMeans(X_sj,na.rm=TRUE)), X_sj)
  }

  # Simulate vessel effects
  Vessel1_vy = array( rnorm( n=max(v_i)*max(t_i), mean=0, sd=Settings[['SigmaVY1']]), dim=c(max(v_i),max(t_i)) )
  Vessel2_vy = array( rnorm( n=max(v_i)*max(t_i), mean=0, sd=Settings[['SigmaVY2']]), dim=c(max(v_i),max(t_i)) )
  Vessel1_v = array( rnorm( n=max(v_i), mean=0, sd=Settings[['SigmaV1']]), dim=c(max(v_i)) )
  Vessel2_v = array( rnorm( n=max(v_i), mean=0, sd=Settings[['SigmaV2']]), dim=c(max(v_i)) )

  # Calculate covariate effect
  eta1_s = as.matrix(X_sj) %*% unlist(Settings[c('Depth1_km','Depth1_km2','Dist1_sqrtkm')])
  eta2_s = as.matrix(X_sj) %*% unlist(Settings[c('Depth2_km','Depth2_km2','Dist2_sqrtkm')])

  # If Nages==1, then simulate as biomass-dynamic model
  if( Settings[["Nages"]]==1 ){
    # Calculate expected values
    P1_st = P2_st = matrix(NA, nrow=nrow(loc_s), ncol=max(t_i) )
    for(t in 1:max(t_i)){
      P1_st[,t] = beta1_t[t] + O1_s + E1_st[,t] + eta1_s
      P2_st[,t] = beta2_t[t] + O2_s + E2_st[,t] + eta2_s
    }

    # Calculate linear predictors
    P1_i = P1_st[cbind(s_i,t_i)] + Vessel1_vy[cbind(v_i,t_i)] + Vessel1_v[v_i]
    P2_i = P2_st[cbind(s_i,t_i)] + Vessel2_vy[cbind(v_i,t_i)] + Vessel2_v[v_i]

    # Calculate true spatio-temporal biomass and annual abundance
    if( Settings[['ObsModel']][2]==0 ){
      B_ast = array( plogis(P1_st) * exp(P2_st) * outer(Extrapolation_List$Area_km2_x,rep(1,max(t_i))), dim=c(1,dim(P1_st)) )
    }
    if( Settings[['ObsModel']][2] %in% c(1,2) ){
      B_ast = array( exp(P1_st) * exp(P2_st) * outer(Extrapolation_List$Area_km2_x,rep(1,max(t_i))), dim=c(1,dim(P2_st)) )
    }
    B_st = B_ast[1,,]

    # Objects specific to Nages==1
    Return = list( "P1_st"=P1_st, "P2_st"=P2_st )

    # Calculate true spatio-temporal biomass and annual abundance
    if( Settings[['ObsModel']][2]==0 ){
      B_st = plogis(P1_st) * exp(P2_st) * outer(Extrapolation_List$Area_km2_x,rep(1,max(t_i)))
    }
    if( Settings[['ObsModel']][2] %in% c(1,2) ){
      B_st = exp(P1_st) * exp(P2_st) * outer(Extrapolation_List$Area_km2_x,rep(1,max(t_i)))
    }
  }

  # If Nages==1, then simulate as biomass-dynamic model
  if( Settings[["Nages"]]>=2 ){
    if( !(Settings[['ObsModel']][2] %in% c(1,2)) ) stop("If using age-structured simulation, please use 'Settings[['ObsModel']][2]=1'")

    # Deviations in initial age-structure
    E1_sa = matrix(NA, nrow=nrow(loc_s), ncol=Settings[["Nages"]] )
    for(a in 2:Settings[["Nages"]]){
      E1_sa[,a] = RFsim(model=model_E1, x=loc_s[,'E_km'], y=loc_s[,'N_km'], standardize=standardize_fields)
      message( "Finished variation for age ", a, " in year 1" )
    }

    # Biomass at age
    L_a = Settings[["Linf"]] * (1 - exp(-Settings[["K"]] * 1:Settings[["Nages"]]) )
    W_a = Settings[["W_alpha"]] * L_a^Settings[["W_beta"]]

    # Selectivity at age
    Selex_A50_v = rnorm( max(v_i), mean=Settings[["Selex_A50_mean"]], sd=Settings[["Selex_A50_sd"]] )
    Selex_av = NULL
    for( v in 1:max(v_i)) Selex_av = cbind( Selex_av, plogis(1:Settings[["Nages"]], location=Selex_A50_v[v], scale=Settings[["Selex_Aslope"]]) )

    # Numbers-density at age (N_ast) and Individual weight at age (W_ast) for age 1
    N_ast = W_ast = array(NA, dim=c(Settings[["Nages"]],nrow(loc_s),max(t_i)) )
    for(t in 1:max(t_i)){
      # Individual weight at age
      W_ast[,,t] = W_a %o% exp( beta2_t[t] + O2_s + E2_st[,t] + eta2_s )
      # Numbers-density at age
      N_ast[1,,t] = exp( beta1_t[t] + O1_s + E1_st[,t] + eta1_s )
      for( a in 2:Settings[["Nages"]] ){
        # Initial age-structure
        if( t==1 ){
          N_ast[a,,t] = exp( beta1_t[t] + O1_s + E1_sa[,a] + eta1_s - Settings[["M"]]*(a-1)  )
        }
        if( t>=2 ){
          N_ast[a,,t] = exp(-Settings[["M"]]) * N_ast[a-1,,t-1]
        }
      }
    }

    # Calculate linear predictors
    P1_i = log(Selex_av[cbind(a_i,v_i)] * N_ast[cbind(a_i,s_i,t_i)]) + Vessel1_vy[cbind(v_i,t_i)] + Vessel1_v[v_i]
    P2_i = log(W_ast[cbind(a_i,s_i,t_i)]) + Vessel2_vy[cbind(v_i,t_i)] + Vessel2_v[v_i]

    # Calculate true AVAILABLE spatio-temporal biomass and annual abundance
    Selex_a = plogis(1:Settings[["Nages"]], location=Settings[["Selex_A50_mean"]], scale=Settings[["Selex_Aslope"]])
    B_ast = N_ast * W_ast * (Selex_a %o% Extrapolation_List$Area_km2_x %o% rep(1,max(t_i)))
    B_st = apply( B_ast, MARGIN=2:3, FUN=sum )
    B_at = apply( B_ast, MARGIN=c(1,3), FUN=sum )

    # Objects specific to Nages>=2
    Return = list( "L_a"=L_a, "W_a"=W_a, "Selex_av"=Selex_av, "E1_sa"=E1_sa, "N_ast"=N_ast, "W_ast"=W_ast, "B_ast"=B_ast, "B_at"=B_at )
  }

  # Calculate linked-predictors
  if( Settings[['ObsModel']][2]==0 ){
    R1_i = plogis( P1_i )
    R2_i = w_i * exp( P2_i )
  }
  if( Settings[['ObsModel']][2]==1 ){
    R1_i = 1 - exp( -1*w_i*exp(P1_i) )
    R2_i = w_i*exp(P1_i) / R1_i * exp(P2_i)
  }
  if( Settings[['ObsModel']][2]==2 ){
    R1_i = w_i * exp(P1_i)
    R2_i = exp(P2_i)
  }

  # Simulate catch and assemble data frame
  if( Settings[['ObsModel']][1]==2 ){
    b1_i = rbinom( n=length(R1_i), size=1, prob=R1_i )
    b2_i = rlnorm( n=length(R2_i), meanlog=log(R2_i)-Settings[['SigmaM']]^2/2, sdlog=Settings[['SigmaM']])
    b_i = b1_i * b2_i
  }
  if( Settings[['ObsModel']][1]==8 ){
    if( Settings[['ObsModel']][2]!=2 ) stop("Must use 'ObsModel=c(8,2)'")
    b_i = rep(NA,length(R1_i))
    for(i in 1:length(b_i)) b_i[i] = rPoisGam( n=1, shape=Settings[['SigmaM']], scale=R2_i[i], intensity=R1_i[i] )
  }
  if( Settings[['ObsModel']][1]==10 ){
    if( Settings[['ObsModel']][2]!=2 ) stop("Must use 'ObsModel=c(10,2)'")
    b_i = rep(NA,length(R1_i))                            # R1_i(i)*R2_i(i), R1_i(i), invlogit(SigmaM(c_i(i),0))+1.0
    for(i in 1:length(b_i)) b_i[i] = tweedie::rtweedie(n=1, mu=R1_i[i]*R2_i[i], phi=R1_i[i], power=plogis(Settings[['SigmaM']])+1.0)
  }
  Data_Geostat = cbind( "Catch_KG"=b_i, "Year"=t_i, "Vessel"=v_i, "Age"=a_i, "AreaSwept_km2"=w_i, "Lat"=loc_i[,'Lat'], "Lon"=loc_i[,'Lon'] )

  # Calculate true center-of-gravity (COG)
  COG_tm = array(NA, dim=c(max(t_i),4), dimnames=list(NULL,c("Lat","Lon","E_km","N_km")))
  for( tI in 1:nrow(COG_tm)){
  for( mI in 1:ncol(COG_tm)){
    COG_tm[tI,mI] = weighted.mean( x=Extrapolation_List$Data_Extrap[,colnames(COG_tm)[mI]], w=B_st[,tI] )
  }}

  # Calculate stratified biomass
  B_tl = array(NA, dim=c(max(t_i),ncol(Extrapolation_List$a_el)), dimnames=list(NULL,colnames(Extrapolation_List$a_el)))
  for( l in 1:ncol(B_tl) ){
    B_tl[,l] = colSums( B_st * outer(Extrapolation_List$a_el[,l]/Extrapolation_List$Area_km2_x,rep(1,max(t_i))), na.rm=TRUE )
  }

  # Timing
  time_for_simulation = Sys.time()-start_time
  message( "Total time: ", time_for_simulation )

  # Objects shared by both
  Return = c( Return, list("Data_Geostat"=data.frame(Data_Geostat), "B_tl"=B_tl, "B_st"=B_st, "COG_tm"=COG_tm, "beta1_t"=beta1_t,
    "beta2_t"=beta2_t, "O1_s"=O1_s, "O2_s"=O2_s, "E1_st"=E1_st, "Vessel1_vy"=Vessel1_vy, "Vessel2_vy"=Vessel2_vy, "eta1_s"=eta1_s, "eta2_s"=eta2_s,
    "Vessel1_v"=Vessel1_v, "Vessel2_v"=Vessel2_v, "E2_st"=E2_st, "s_i"=s_i, "t_i"=t_i, "a_i"=a_i, "w_i"=w_i,
    "P1_i"=P1_i, "P2_i"=P2_i, "R1_i"=R1_i, "R2_i"=R2_i,
    "time_for_simulation"=time_for_simulation, "Sim_Settings"=Settings) )
  if( exists("b1_i")) Return = c(Return, list("b1_i"=b1_i, "b2_i"=b2_i))
  return( Return )
}

#' Diagnostic QQ function
#'
#' @param TmbData TMB Model input data list
#' @param Report TMB Model output data list
#' @param save_dir Directory to save plots, if NULL is specified then do not save plot
#' @param FileName_PP If NULL is specified then do not save this type of plot
#' @param FileName_Phist If NULL is specified then do not save this type of plot
#' @param FileName_QQ If NULL is specified then do not save this type of plot
#' @param FileName_Qhist If NULL is specified then do not save this type of plot
#' @return A list containing results for each specified categories
plot_quantile_diagnostic <- function(TmbData,
                  Report,
                  DateFile=paste0(getwd(),"/"),
                  save_dir=paste0(DateFile,"/QQ_Fn/"),
                  FileName_PP="Posterior_Predictive",
                  FileName_Phist="Posterior_Predictive-Histogram",
                  FileName_QQ="Q-Q_plot",
                  FileName_Qhist="Q-Q_hist"){

    # Retrieve data based on model type
    if("n_e" %in% names(TmbData)){
        # VAST Versions 3.0.0, 4.0.0: names n_e, e_i and ObsModel_ez are added to support for error distns
        n_e <- TmbData$n_e
        e_i <- TmbData$e_i
        ObsModel_ez <- TmbData$ObsModel_ez
        sigmaM <- Report$SigmaM
    }else if("n_c" %in% names(TmbData)){
        # VAST Version < 3.0.0: names n_c, c_i and ObsModel are used to group by categories
        n_e <- TmbData$n_c
        e_i <- TmbData$c_i
        ObsModel_ez <- rep(1, n_e) %o% TmbData$ObsModel
        sigmaM <- Report$SigmaM
        if(is.vector(sigmaM))
            # VAST Versions 1.0.0 and 1.1.0
            sigmaM <- rep(1, n_e) %o% Report$SigmaM
    }else{
        # SpatialDeltaGLMM
        n_e <- 1
        e_i <- rep(0,TmbData$n_i)
        ObsModel_ez <- matrix(TmbData$ObsModel, nrow = 1)
        sigmaM <- matrix(Report$SigmaM, nrow = 1)
    }
    # else if("n_p" %in% names(TmbData)){
    #     # MIST: not supported, as MIST uses different definitions for ObsModel
    #     return(list(message="The function does not support MIST yet."))
    # }

    # Check data
    if(nlevels(as.factor(e_i))!=n_e) stop("Error in e_i: nlevels does not agree with n_e")
    if(nrow(as.matrix(ObsModel_ez))!=n_e) stop("Error in ObsModel_ez: nrow does not agree with n_e")
    if(nrow(as.matrix(sigmaM))!=n_e) stop("Error in sigmaM: nrow does not agree with n_e")

    # Check save_dir
    dir.create(save_dir, recursive=TRUE, showWarnings = FALSE)
    if(!dir.exists(save_dir)) stop(paste0("Wrong directory, cannot save plots: ", save_dir))

    # Return list
    Return <- vector("list", length = n_e)

    # utility function
    pow = function(a,b) a^b

    # Loop through each group (plot functions remain unchanged from previous version)
    for(i_e in 1:n_e){

      # Generate plot names
      if(!is.null(save_dir)){
        if(!is.null(FileName_PP)) save_PP=paste0(save_dir,"/",FileName_PP,"-",i_e,".jpg")
        if(!is.null(FileName_Phist)) save_Phist=paste0(save_dir,"/",FileName_Phist, "-",i_e,".jpg")
        if(!is.null(FileName_QQ)) save_QQ=paste0(save_dir,"/",FileName_QQ,"-",i_e,".jpg")
        if(!is.null(FileName_Qhist)) save_Qhist=paste0(save_dir,"/",FileName_Qhist,"-",i_e,".jpg")
      }

      # Find where b_i > 0 within category i_e
      Which = which(strip_units(TmbData$b_i) > 0 & e_i == (i_e-1))
      Q = rep(NA, length(Which)) # vector to track quantiles for each observation
      y = array(NA, dim=c(length(Which),1000)) # matrix to store samples
      pred_y = var_y = rep(NA, length(Which) ) # vector to track quantiles for each observation

      # Calculate pred_y
      # I can't use R2_i anymore because interpretation changed around March 9, 2017 (due to area-swept change in Poisson-process and Tweedie functions)
      # However, I CAN use P2_i, which has a stable definition over time (as a linear predictor)
      if( length(ObsModel_ez[i_e,])>=2 && ObsModel_ez[i_e,2]==2 ){
        Return[[i_e]] = list("type"=ObsModel_ez[i_e,], message="QQ not set up for Tweedie distribution")
        next
      }
      if( !(ObsModel_ez[i_e,1] %in% c(1,2)) ){
        Return[[i_e]] = list("type"=ObsModel_ez[i_e,], message="QQ not working except for when using a Gamma or Lognormal distribution")
        next
      }
      if( length(ObsModel_ez[i_e,])==1 || ObsModel_ez[i_e,2]%in%c(0,3) ){
        for(ObsI in 1:length(Which)){
          pred_y[ObsI] = strip_units(TmbData$a_i[Which[ObsI]]) * exp(Report$P2_i[Which[ObsI]])
        }
      }
      if( length(ObsModel_ez[i_e,])>=2 && ObsModel_ez[i_e,2]%in%c(1,4) ){
        for(ObsI in 1:length(Which)){
          if(sigmaM[e_i[Which[ObsI]]+1,3]!=1) stop("`QQ_Fn` will not work with Poisson-link delta model across all VAST versions given values for turned-off parameters")
          R1_i = 1 - exp( -1 * sigmaM[e_i[Which[ObsI]]+1,3] * strip_units(TmbData$a_i[Which[ObsI]]) * exp(Report$P1_i[Which[ObsI]]) )
          pred_y[ObsI] = strip_units(TmbData$a_i[Which[ObsI]]) * exp(Report$P1_i[Which[ObsI]]) / R1_i * exp(Report$P2_i[Which[ObsI]]);
        }
      }

      # Simulate quantiles for different distributions: Loop through observations
      for(ObsI in 1:length(Which)){
        if(ObsModel_ez[i_e,1]==1){
          y[ObsI,] = rlnorm(n=ncol(y), meanlog=log(pred_y[ObsI])-pow(sigmaM[i_e,1],2)/2, sdlog=sigmaM[i_e,1])   # Plotting in log-space
          Q[ObsI] = plnorm(q=strip_units(TmbData$b_i[Which[ObsI]]), meanlog=log(pred_y[ObsI])-pow(sigmaM[i_e,1],2)/2, sdlog=sigmaM[i_e,1])
        }
        if(ObsModel_ez[i_e,1]==2){
          b = pow(sigmaM[i_e, 1],2) * pred_y[ObsI];
          y[ObsI,] = rgamma(n=ncol(y), shape=1/pow(sigmaM[i_e,1],2), scale=b)
          Q[ObsI] = pgamma(q=strip_units(TmbData$b_i[Which[ObsI]]), shape=1/pow(sigmaM[i_e,1],2), scale=b)
        }
      }

      # Make plot while calculating posterior predictives
      #if(!is.null(FileName_PP) & !is.null(save_dir)) jpeg(save_PP, width=10, height=3, res=200, units="in")
      #par(mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
      #plot(TmbData$b_i[Which], ylab="", xlab="", log="y", main="", col="blue")

      # Add results to plot: Loop through observations
      for(ObsI in 1:length(Which)){
        var_y[ObsI] = var( y[ObsI,] )
        Quantiles = quantile(y[ObsI,],prob=c(0.025,0.25,0.75,0.975))
        lines(x=c(ObsI,ObsI), y=Quantiles[2:3], lwd=2)
        lines(x=c(ObsI,ObsI), y=Quantiles[c(1,4)], lwd=1,lty="dotted")
        if( strip_units(TmbData$b_i[Which[ObsI]])>max(Quantiles) | strip_units(TmbData$b_i[Which[ObsI]])<min(Quantiles)){
          points(x=ObsI,y=TmbData$b_i[Which[ObsI]],pch=4,col="red",cex=2)
        }
      }
      if(!is.null(FileName_PP) & !is.null(save_dir)) dev.off()

      # Q-Q plot
      if(!is.null(FileName_Phist) & !is.null(save_dir)) jpeg(save_Phist, width=4, height=4, res=200, units="in")
      par(mfrow=c(1,1), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
      Qtemp = na.omit(Q)
      Order = order(Qtemp)
      plot(x=seq(0,1,length=length(Order)), y=Qtemp[Order], main="Q-Q plot", xlab="Uniform", ylab="Empirical", type="l", lwd=3)
      abline(a=0,b=1)
      if(!is.null(FileName_Phist) & !is.null(save_dir)) dev.off()

      # Aggregate predictive distribution
      #if(!is.null(FileName_QQ) & !is.null(save_dir)) jpeg(save_QQ, width=4, height=4, res=200, units="in")
      #par(mfrow=c(1,1), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
      #hist( log(y), main="Aggregate predictive dist.", xlab="log(Obs)", ylab="Density")
      #if(!is.null(FileName_QQ) & !is.null(save_dir)) dev.off()

      # Quantile histogram
      #if(!is.null(FileName_Qhist) & !is.null(save_dir)) jpeg(save_Qhist, width=4, height=4, res=200, units="in")
      #par(mfrow=c(1,1), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
      #hist(na.omit(Q), main="Quantile_histogram", xlab="Quantile", ylab="Number")
      #if(!is.null(FileName_Qhist) & !is.null(save_dir)) dev.off()

      # Return stuff
      Return[[i_e]] = list("type"=ObsModel_ez[i_e,], "Q"=Q, "var_y"=var_y, "pred_y"=pred_y, "y"=y )
    }

    if(length(Return)==1) Return <- Return[[1]] # single species model
    return( Return )
}

#' @title
#' Plot Pearson residuals on map
#'
#' @description
#' \code{plot_residuals} shows average Pearson residual for every knot for encounter probability and positive catch rate components
#'
#' @inheritParams plot_maps
#'
#' @param Lat_i Latitude for every observation \code{i}
#' @param Lon_i Longitude for every observation \code{i}
#' @param Q Output from \code{QQ_Fn}
#' @param savdir directory to use when saving results
#' @param TmbData a tagged list of data inputs generated by \code{make_data}
#' @param zrange lower and upper bounds for plotted value of Pearson residuals, used to change plotting scale to eliminate influence of outliers on plotting scale
#' @param ... arguments passed to \code{FisshStatsUtils::plot_variable}
#'
#' @return A tagged list of Pearson residuals
#' \describe{
#'   \item{Q1_xt}{Matrix of average residuals for encounter/non-encounter component by site \code{x} and year \code{t}}
#'   \item{Q2_xt}{Matrix of average residuals for positive-catch-rate component by site \code{x} and year \code{t}}
#' }
plot_residuals = function( Lat_i, Lon_i, TmbData, Report, Q, projargs='+proj=longlat',
         working_dir=paste0(getwd(),"/"), spatial_list, extrapolation_list,
         Year_Set=NULL, Years2Include=NULL, zrange, ... ){

  ##################
  # Presence-absence
  # http://data.princeton.edu/wws509/notes/c3s8.html
  ##################

  # Add t_iz if missing (e.g., from earlier version of VAST, or SpatialDeltaGLMM)
  if( !("t_iz" %in% names(TmbData)) ){
    TmbData$t_iz = matrix( TmbData$t_i, ncol=1 )
  }

  # Add in t_yz if missing (e.g., from earlier version of VAST, or SpatialDeltaGLMM)
  if( !("t_yz" %in% names(TmbData)) ){
    TmbData$t_yz = matrix(1:TmbData$n_t - 1, ncol=1)
  }

  # Extract binomial stuff for encounter-nonencounter data
  exp_rate_xy = obs_rate_xy = total_num_xy = exp_num_xy = obs_num_xy = matrix(NA, nrow=spatial_list$n_x, ncol=nrow(TmbData$t_yz) )
  for( yI in 1:nrow(TmbData$t_yz) ){
    which_i_in_y = ( TmbData$t_iz == outer(rep(1,TmbData$n_i),TmbData$t_yz[yI,]) )
    which_i_in_y = which( apply(which_i_in_y,MARGIN=1,FUN=all) )
    if( length(which_i_in_y)>0 ){
      exp_rate_xy[,yI] = tapply( Report$R1_i[which_i_in_y], INDEX=factor(spatial_list$knot_i[which_i_in_y],levels=1:spatial_list$n_x), FUN=mean )
      obs_rate_xy[,yI] = tapply( strip_units(TmbData$b_i[which_i_in_y])>0, INDEX=factor(spatial_list$knot_i[which_i_in_y],levels=1:spatial_list$n_x), FUN=mean )
      total_num_xy[,yI] = tapply( strip_units(TmbData$b_i[which_i_in_y]), INDEX=factor(spatial_list$knot_i[which_i_in_y],levels=1:spatial_list$n_x), FUN=length )
    }else{
      total_num_xy[,yI] = 0
    }
    exp_num_xy = exp_rate_xy * total_num_xy
    obs_num_xy = obs_rate_xy * total_num_xy
  }

  # Method #1 -- Binomial cumulative function
  #Q1_xt = pbinom( obs_num_xt, size=total_num_xt, prob=exp_rate_xt )

  # Method #2 -- Pearson residuals
  Q1_xy = (obs_num_xy - exp_num_xy) / sqrt(  exp_num_xy*(total_num_xy-exp_num_xy)/total_num_xy )

  # Method #3 -- Deviance residuals
  #Q1_xt = 2*(obs_num_xt*log(obs_num_xt/exp_num_xt) + (total_num_xt-obs_num_xt)*log((total_num_xt-obs_num_xt)/(total_num_xt-exp_num_xt)))
  #Q1_xt = ifelse( obs_num_xt==0, 2*((total_num_xt-obs_num_xt)*log((total_num_xt-obs_num_xt)/(total_num_xt-exp_num_xt))), Q1_xt)
  #Q1_xt = ifelse( (total_num_xt-obs_num_xt)==0, 2*(obs_num_xt*log(obs_num_xt/exp_num_xt)), Q1_xt)

  ##################
  # Positive catch rates
  ##################

  # Extract quantile for positive catch rates
  #Q_i = Q[["Q"]]
  which_pos = which( strip_units(TmbData$b_i)>0 )
  bvar_ipos = bpred_ipos = NULL
  # Univariate Q interface
  if( all(c("var_y","pred_y") %in% names(Q)) ){
    bvar_ipos = Q[["var_y"]]   # Change name to avoid naming-convention of y with reporting-interval
    bpred_ipos = Q[["pred_y"]]  # Change name to avoid naming-convention of y with reporting-interval
  }
  # Multivariate Q interface
  if( all(c("var_y","pred_y") %in% names(Q[[1]])) ){
    bvar_ipos = bpred_ipos = rep(NA, length=length(which_pos))
    for(i_e in 1:length(Q)){
      which_pos_and_e = which(TmbData$e_i[which_pos]==(i_e-1))
      bvar_ipos[which_pos_and_e] = Q[[i_e]][["var_y"]]
      bpred_ipos[which_pos_and_e] = Q[[i_e]][["pred_y"]]
    }
  }
  if( is.null(bvar_ipos) & is.null(bpred_ipos) ){
    stop("Something is wrong with `Q` input")
  }

  ### Method #1 -- chi-squared transformation of cumulative function
  # Convert to Chi-squared distribution
  #Chisq_x = tapply( Q_i, INDEX=spatial_list$knot_i[which(TmbData$b_i>0)], FUN=function(vec){sum(-2*log(vec))} )
  # Calculate d.f. for Chi-squared distribution
  #DF_x = 2*tapply( Q_i, INDEX=spatial_list$knot_i[which(TmbData$b_i>0)], FUN=length )
  # Convert back to uniform distribution
  #P_x = pchisq( q=Chisq_x, df=DF_x )

  ### Method #2 -- average quantile
  #Q2_x = tapply( Q_i, INDEX=spatial_list$knot_i[which(TmbData$b_i>0)], FUN=mean )
  #Q2_xt = tapply( Q_i, INDEX=list(spatial_list$knot_i[which(TmbData$b_i>0)],TmbData$t_i[which(TmbData$b_i>0)]), FUN=mean )

  ### Method #3 -- Pearson residuals
  sum_obs_xy = sum_exp_xy = var_exp_xy = matrix(NA, nrow=spatial_list$n_x, ncol=nrow(TmbData$t_yz) )
  for( yI in 1:nrow(TmbData$t_yz) ){
    which_i_in_y = ( TmbData$t_iz == outer(rep(1,TmbData$n_i),TmbData$t_yz[yI,]) )
    which_i_in_y = which( apply(which_i_in_y,MARGIN=1,FUN=all) )
    which_i_in_y_and_pos = intersect( which_i_in_y, which_pos )
    which_ipos_in_y = ( TmbData$t_iz[which_pos,] == outer(rep(1,length(which_pos)),TmbData$t_yz[yI,]) )
    which_ipos_in_y = which( apply(which_ipos_in_y,MARGIN=1,FUN=all) )
    if( length(which_i_in_y_and_pos)>0 ){
      sum_obs_xy[,yI] = tapply( strip_units(TmbData$b_i[which_i_in_y_and_pos]), INDEX=factor(spatial_list$knot_i[which_i_in_y_and_pos],levels=1:spatial_list$n_x), FUN=sum )
      sum_exp_xy[,yI] = tapply( bpred_ipos[which_ipos_in_y], INDEX=factor(spatial_list$knot_i[which_i_in_y_and_pos],levels=1:spatial_list$n_x), FUN=sum )
      var_exp_xy[,yI] = tapply( bvar_ipos[which_ipos_in_y], INDEX=factor(spatial_list$knot_i[which_i_in_y_and_pos],levels=1:spatial_list$n_x), FUN=sum )
    }
  }
  Q2_xy = (sum_obs_xy - sum_exp_xy) / sqrt(var_exp_xy)

  #################
  # Plots
  #################

  if( !is.null(working_dir) ){
    for( zI in 1:2 ){
      Q_xy = list( Q1_xy, Q2_xy )[[zI]]
      if( !missing(zrange) ){
        Q_xy = ifelse( Q_xy<zrange[1], zrange[1], Q_xy )
        Q_xy = ifelse( Q_xy>zrange[2], zrange[2], Q_xy )
        zlim = zrange
      }else{
        zlim = c(-1,1) * ceiling(max(abs(Q_xy),na.rm=TRUE))
      }
      #Q_xt = ifelse( abs(Q_xt)>3, 3*sign(Q_xt), Q_xt )

      # Plot settings
      Col = colorRampPalette(colors=c("blue","white","red"))
      textmargin = "Pearson residual"
      plot_code = c("pearson_residuals_1", "pearson_residuals_2")[zI]

      # Spatial information
      # Aggregate residual-values to knots regardless of value for fine_scale
      x2i = spatial_list$NN_Extrap$nn.idx[,1]
      Include = strip_units(extrapolation_list[["Area_km2_x"]])>0 & strip_units(extrapolation_list[["a_el"]][,1])>0
      DF = cbind( extrapolation_list$Data_Extrap[,c('Lon','Lat')], "x2i"=x2i, "Include"=Include )

      # Fill in labels
      if( is.null(Year_Set) ) Year_Set = 1:ncol(Q_xy)
      if( is.null(Years2Include) ) Years2Include = 1:ncol(Q_xy)

      # Make plots
      plot_args = plot_variable( Y_gt=ifelse(is.na(Q_xy),mean(zlim),Q_xy), map_list=list("PlotDF"=DF), projargs=projargs, working_dir=working_dir,
        panel_labels=Year_Set[Years2Include], file_name=plot_code, zlim=zlim, col=Col, ... )
    }
  }

  #################
  # Returns
  #################

  Return = list( "Q1_xy"=Q1_xy, "Q2_xy"=Q2_xy  )  # , "obs_num_xy"=obs_num_xy, "exp_num_xy"=exp_num_xy, "total_num_xy"=total_num_xy
  return( invisible(Return) )
}


#' @title
#' Check predicted encounter probability against observed encounter frequency
#'
#' @description
#' \code{plot_encounter_diagnostic} is a diagnostic function for checking validity of the encounter-probability component of a spatio-temporal model
#'
#' @inheritParams plot_maps
#' @inheritParams plot_biomass_index
#' @param ... arguments passed to \code{par}
#'
#' @return Return Tagged list of output
#' \describe{
#'   \item{Diag_i}{Diagnostic output for each sample \code{i}}
#'   \item{Diag_z}{Diagnostic output for each bin \code{z}}
#' }
plot_encounter_diagnostic = function( Report, Data_Geostat, cutpoints_z=seq(0,1,length=21), interval_width=1.96, DirName=paste0(getwd(),"/"),
  PlotName="Diag--Encounter_prob.png", ... ){

  # Get bin for each datum
  z_i = cut( Report$R1_i, breaks=cutpoints_z, include.lowest=TRUE )
  midpoints_z = rowMeans( cbind(cutpoints_z[-1],cutpoints_z[-length(cutpoints_z)]) )

  # Get encounter frequency for each bin
  freq_z = tapply( ifelse(Data_Geostat[,'Catch_KG']>0,1,0), INDEX=z_i, FUN=mean )

  # Get expectation given model
  num_z = tapply( Report$R1_i, INDEX=z_i, FUN=length )
  mean_z = tapply( Report$R1_i, INDEX=z_i, FUN=mean )
  var_z = tapply( Report$R1_i, INDEX=z_i, FUN=function(vec){sum(vec*(1-vec))} )
  sd_mean_z = sqrt(var_z / num_z^2)

  # Plot
  Par = list( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02, ... )
  png( file=paste0(DirName,"/",PlotName), width=5, height=5, res=200, units="in")
    par( Par )
    plot( x=midpoints_z, y=freq_z, pch=20, cex=1.2, xlim=c(0,1), ylim=c(0,1), xlab="Predicted encounter probability", ylab="Observed encounter frequency" )
    plot_lines( x=midpoints_z[which(!is.na(mean_z))], y=mean_z[which(!is.na(mean_z))], ybounds=(mean_z%o%c(1,1)+sd_mean_z%o%c(-interval_width,interval_width))[which(!is.na(mean_z)),,drop=FALSE], lwd=2, bounds_type="shading", col_bounds=rgb(1,0,0,0.2), col="red" )
    abline(a=0, b=1, lty="dotted", lwd=2 )
    legend( "topleft", legend=c("Observed","Predicted"), fill=c("black","red"), bty="n")
  dev.off()

  # Return stuff
  Return = NULL
  Return[["Diag_i"]] = cbind("b_i"=Data_Geostat[,'Catch_KG'], "z_i"=z_i)
  Return[["Diag_z"]] = cbind("midpoints_z"=midpoints_z, "freq_z"=freq_z, "num_z"=num_z, "mean_z"=mean_z, "var_z"=var_z, "sd_mean_z"=sd_mean_z)
  return( invisible(Return) )
}

