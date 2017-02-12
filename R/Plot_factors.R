
#' Plot factor-decomposition of covariance
#'
#' \code{Plot_factors} plots factor loadings, average spatial factors, and spatio-temporal factors
#'
#' @inheritParams Plot_Overdispersion
#' @inheritParams Summarize_Covariance
#' @param Year_Set plotting-names for time dimension
#' @param mapdetails_list output from \code{SpatialDeltaGLMM::MapDetails_Fn}
#' @param Dim_year Plotting dimension (row,column) for plot of years (default: square with sufficient size for number of years)
#' @param Dim_species Plotting dimension (row,column) for plot of categories (default: square with sufficient size for number of categories)
#' @param plotdir directory for saving plots
#' @param land_color color for filling in land (use \code{land_color=rgb(0,0,0,alpha=0)} for transparent land)

#' @export
Plot_factors = function( Report, ParHat, Data, SD, Year_Set=NULL, category_names=NULL,
  mapdetails_list=NULL, Dim_year=NULL, Dim_species=NULL, plotdir=paste0(getwd(),"/"), land_color="grey" ){

  # Fill in missing inputs
  if( "D_xct" %in% names(Report) ){
    if( is.null(Year_Set) ) Year_Set = 1:dim(Report$D_xct)[3]
    if( is.null(category_names) ) category_names = 1:dim(Report$D_xct)[2]
  }
  if( "D_xcy" %in% names(Report) ){
    if( is.null(Year_Set) ) Year_Set = 1:dim(Report$D_xcy)[3]
    if( is.null(category_names) ) category_names = 1:dim(Report$D_xcy)[2]
    Report[["D_xct"]] = Report[["D_xcy"]]
  }

  # Dimensions for plotting
  Dim = function( num ) c(ceiling(sqrt(num)), ceiling(num/ceiling(sqrt(num))) )
  Dim_year = Dim(length(Year_Set))
  Dim_species = Dim(length(category_names))

  # Extract covariance
  #Cov_List = Summarize_Covariance( Report=Report, ParHat=ParHat, Data=Data, SD=SD, category_names=category_names, figname=NULL )

  # Extract loadings matrices (more numerically stable than extracting covariances, and then re-creating Cholesky)
  L_list = vector("list", length=4)
  names(L_list) = c("Omega1", "Epsilon1", "Omega2", "Epsilon2")

  # Loop through
  for(i in 1:4){
    # Variable names
    Par_name = names(L_list)[i]
    if(Par_name == "Omega1") Var_name = "Omegainput1_sf"
    if(Par_name == "Epsilon1") Var_name = "Epsiloninput1_sft"
    if(Par_name == "Omega2") Var_name = "Omegainput2_sf"
    if(Par_name == "Epsilon2") Var_name = "Epsiloninput2_sft"

    # Continue if component is included
    if( Data[["FieldConfig"]][[Par_name]]>0 ){
      # Get loadings matrix
      L_list[[i]] = calc_cov( L_z=ParHat[[paste0("L_",tolower(Par_name),"_z")]], n_f=Data[["FieldConfig"]][[Par_name]], n_c=Data$n_c, returntype="loadings_matrix" )
      rownames(L_list[[i]]) = category_names

      # Get covariance # SpatialDFA::
      # Cov_jj=Cov_List[[paste0("Cov_",tolower(Par_name))]][,,'Estimate']; Psi=Report[[Var_name]]; RotationMethod="PCA"; testcutoff=1e-4
      #Var_rot = Rotate_Fn( Cov_jj=Cov_List[[paste0("Cov_",tolower(Par_name))]][,,'Estimate'], Psi=Report[[Var_name]], RotationMethod="PCA", testcutoff=1e-4 )
      # L_pj=L_list[[i]]; Psi=Report[[Var_name]]; RotationMethod="PCA"; testcutoff=1e-4
      Var_rot = Rotate_Fn( L_pj=L_list[[i]], Psi=Report[[Var_name]], RotationMethod="PCA", testcutoff=1e-4 )

      # Plot loadings
      Dim_factor = Dim(Data[["FieldConfig"]][[Par_name]])
      png( file=paste0(plotdir,"Factor_loadings--",Par_name,".png"), width=Dim_factor[2]*4, height=Dim_factor[1]*4, units="in", res=200 )
        par( mfrow=Dim_factor, mar=c(0,2,2,0) )
        for( cI in 1:Data[["FieldConfig"]][[Par_name]] ) SpatialDFA::PlotLoadings( L_pj=Var_rot$L_pj_rot, whichfactor=cI )
      dev.off()

      # Plot factors
      if( !is.null(mapdetails_list) ){
        # Plot factors by year
        if( Par_name %in% c("Epsilon1","Epsilon2")){
          # plot_set=3; MappingDetails; Report; Sdreport=NULL; Nknots=Inf; PlotDF; MapSizeRatio=c('Width(in)'=4,'Height(in)'=4); Xlim; Ylim; FileName=paste0(getwd(),"/"); Year_Set=NULL; Years2Include=NULL; Rescale=FALSE; Rotate=0; Format="png"; Res=200; zone=NA; Cex=0.01; add=FALSE; category_names=NULL; textmargin=NULL; pch=NULL; Legend=list("use"=FALSE,"x"=c(10,30),"y"=c(10,30)); mfrow=NULL; plot_legend_fig=TRUE
          # plot_set=c(NA,6,NA,7)[i]; MappingDetails=mapdetails_list[["MappingDetails"]]; Report=list("D_xct"=Report$D_xct,"Epsilon1_sct"=Var_rot$Psi_rot,"Epsilon2_sct"=Var_rot$Psi_rot); PlotDF=mapdetails_list[["PlotDF"]]; MapSizeRatio=mapdetails_list[["MapSizeRatio"]]; Xlim=mapdetails_list[["Xlim"]]; Ylim=mapdetails_list[["Ylim"]]; FileName=plotdir; Year_Set=Year_Set; Rotate=mapdetails_list[["Rotate"]]; category_names=paste0("Factor_",1:length(category_names)); mar=c(0,0,2,0); oma=c(1.5,1.5,0,0); Cex=mapdetails_list[["Cex"]]; cex=1.8; mfrow=Dim_year; cex.main=1.0; Legend=mapdetails_list[["Legend"]]; zone=mapdetails_list[["Zone"]]; plot_legend_fig=FALSE; land_color=land_color
          SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=c(NA,6,NA,7)[i], MappingDetails=mapdetails_list[["MappingDetails"]], Report=list("D_xct"=Var_rot$Psi_rot,"Epsilon1_sct"=Var_rot$Psi_rot,"Epsilon2_sct"=Var_rot$Psi_rot), PlotDF=mapdetails_list[["PlotDF"]], MapSizeRatio=mapdetails_list[["MapSizeRatio"]], Xlim=mapdetails_list[["Xlim"]], Ylim=mapdetails_list[["Ylim"]], FileName=plotdir, Year_Set=Year_Set, Rotate=mapdetails_list[["Rotate"]], category_names=paste0("Factor_",1:length(category_names)), mar=c(0,0,2,0), oma=c(1.5,1.5,0,0), Cex=mapdetails_list[["Cex"]], cex=1.8, mfrow=Dim_year, cex.main=1.0, Legend=mapdetails_list[["Legend"]], zone=mapdetails_list[["Zone"]], plot_legend_fig=FALSE, land_color=land_color)
        }  #

        # Plot average factors across years
        Mat_sf = apply(Var_rot$Psi_rot, MARGIN=1:2, FUN=mean)
        SpatialDeltaGLMM:::PlotMap_Fn( MappingDetails=mapdetails_list[["MappingDetails"]], Mat=Mat_sf, PlotDF=mapdetails_list[["PlotDF"]], MapSizeRatio=mapdetails_list[["MapSizeRatio"]], Xlim=mapdetails_list[["Xlim"]], Ylim=mapdetails_list[["Ylim"]], FileName=paste0(plotdir,"Factor_maps--",Par_name), Year_Set=paste0("Factor_",1:ncol(Mat_sf)), Rotate=mapdetails_list[["Rotate"]], zone=mapdetails_list[["Zone"]], mar=c(0,0,2,0), oma=c(2.5,2.5,0,0), Cex=0.01, mfrow=Dim_factor, pch=20, Legend=mapdetails_list[["Legend"]], plot_legend_fig=FALSE, land_color=land_color)
      }
    }else{
      L_list[[i]] = NULL
    }
  }
}
