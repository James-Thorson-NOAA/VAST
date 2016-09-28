
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
#' @param Dim_species Plotting dimension (row,column) for plot of categories (default: square with sufficient size for number of categories)

#' @export
Plot_factors = function( Report, ParHat, Data, SD, Year_Set=1:dim(Report$D_xct)[3], category_names=1:dim(Report$D_xct)[2],
  mapdetails_list=NULL, Dim_year=NULL, Dim_species=NULL, plotdir=paste0(getwd(),"/") ){

  # Dimensions for plotting
  Dim = function( num ) c(ceiling(sqrt(num)), ceiling(num/ceiling(sqrt(num))) )
  Dim_year = Dim(length(Year_Set))
  Dim_species = Dim(length(category_names))

  # Extract covariance
  Cov_List = Summarize_Covariance( Report=Report, ParHat=ParHat, Data=Data, SD=SD, category_names=category_names, figname=NULL )

  # Loop through
  for(i in 1:4){
    # Variable names
    Par_name = c("Omega1", "Epsilon1", "Omega2", "Epsilon2")[i]
    if(Par_name %in% c("Epsilon1","Epsilon2")) Var_name = paste0(Par_name,"_sct")
    if(Par_name %in% c("Omega1","Omega2")) Var_name = paste0(Par_name,"_st")

    # Continue if component is included
    if( FieldConfig[[Par_name]]>0 ){
      # Get covariance
      # Cov_jj=Cov_List[[paste0("Cov_",tolower(Par_name))]][,,1]; Psi=Report[[Var_name]]; RotationMethod="PCA"; testcutoff=1e-4
      Var_rot = SpatialDFA::Rotate_Fn( Cov_jj=Cov_List[[paste0("Cov_",tolower(Par_name))]][,,1], Psi=Report[[Var_name]], RotationMethod="PCA", testcutoff=1e-4 )

      # Plot factors by year
      if( Par_name %in% c("Epsilon1","Epsilon2")){
        SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=c(NA,6,NA,7)[i], MappingDetails=mapdetails_list[["MappingDetails"]], Report=list("D_xct"=Report$D_xct,"Epsilon1_sct"=Var_rot$Psi_rot,"Epsilon2_sct"=Var_rot$Psi_rot), PlotDF=mapdetails_list[["PlotDF"]], MapSizeRatio=mapdetails_list[["MapSizeRatio"]], Xlim=mapdetails_list[["Xlim"]], Ylim=mapdetails_list[["Ylim"]], FileName=plotdir, Year_Set=Year_Set, Rotate=mapdetails_list[["Rotate"]], category_names=paste0("Factor_",1:length(category_names)), mar=c(0,0,2,0), oma=c(1.5,1.5,0,0), Cex=mapdetails_list[["Cex"]], cex=1.8, mfrow=Dim_year, cex.main=1.0, Legend=mapdetails_list[["Legend"]], zone=mapdetails_list[["Zone"]], plot_legend_fig=FALSE)
      }  #

      # Plot loadings
      png( file=paste0(plotdir,"Factor_loadings--",Par_name,".png"), width=Dim_species[2]*2.5, height=Dim_species[1]*2.5, units="in", res=200 )
        par( mfrow=Dim_species, mar=c(0,2,2,0) )
        for( cI in 1:length(Species_Set)) SpatialDFA::PlotLoadings( Var_rot$L_pj_rot, whichfactor=cI )
      dev.off()

      # Plot loadings
      if(Par_name %in% c("Epsilon1","Epsilon2")) Mat_sc = apply(Var_rot$Psi_rot, MARGIN=1:2, FUN=mean)
      if(Par_name %in% c("Omega1","Omega2")) Mat_sc = Var_rot$Psi_rot
      SpatialDeltaGLMM:::PlotMap_Fn( MappingDetails=mapdetails_list[["MappingDetails"]], Mat=Mat_sc, PlotDF=mapdetails_list[["PlotDF"]], MapSizeRatio=mapdetails_list[["MapSizeRatio"]], Xlim=mapdetails_list[["Xlim"]], Ylim=mapdetails_list[["Ylim"]], FileName=paste0(plotdir,"Factor_maps--",Par_name), Year_Set=paste0("Factor_",1:length(Species_Set)), Rotate=mapdetails_list[["Rotate"]], zone=mapdetails_list[["Zone"]], mar=c(0,0,2,0), oma=c(2.5,2.5,0,0), Cex=0.01, mfrow=Dim_species, pch=20, Legend=list(use=TRUE, x=c(65,75), y=c(35,65)), plot_legend_fig=FALSE)
    }  #
  }
}
