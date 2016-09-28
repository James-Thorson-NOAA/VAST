
#' Plot factor-decomposition of covariance
#'
#' \code{Plot_factors} plots factor loadings, average spatial factors, and spatio-temporal factors
#'
#' @inheritParams Plot_Overdispersion
#' @inheritParams Summarize_Covariance
#' @param Year_Set plotting-names for time dimension
#' @param Dim_year Plotting dimension (row,column) for plot of years (default: square with sufficient size for number of years)
#' @param Dim_species Plotting dimension (row,column) for plot of categories (default: square with sufficient size for number of categories)

#' @export
Plot_factors = function( Report, ParHat, Data, SD, Year_Set=1:dim(Report$D_xct)[3], category_names=1:dim(Report$D_xct)[2],
  Dim_year=NULL, Dim_species=NULL, plotdir=paste0(getwd(),"/") ){

  # Dimensions for plotting
  Dim = function( num ) c(ceiling(sqrt(num)), ceiling(num/ceiling(sqrt(num))) )
  Dim_year = Dim(length(Year_Set))
  Dim_species = Dim(length(category_names))

  # Extract covariance
  Cov_List = Summarize_Covariance( Report=Report, ParHat=ParHat, Data=Data, SD=SD, names_set=category_names )

  # Loop through
  for(i in 1:4){
    # Variable names
    Par_name = c("Omega1", "Epsilon1", "Omega2", "Epsilon2")[i]
    if(Par_name %in% c("Epsilon1","Epsilon2")) Var_name = paste0(Par_name,"_sct")
    if(Par_name %in% c("Omega1","Omega2")) Var_name = paste0(Par_name,"_st")

    # Continue if component is included
    if( FieldConfig[[Par_name]]>0 ){
      # Get covariance
      Var_rot = Rotate_Fn( Cov_jj=Cov_List[[paste0("Cov_",tolower(Par_name))]][,,1], Psi=Report[[Var_name]], testcutoff=1e-4 )

      # Plot factors by year
      if( Par_name %in% c("Epsilon1","Epsilon2")){
        PlotResultsOnMap_Fn(plot_set=c(NA,6,NA,7)[i], MappingDetails=MapDetails_List[["MappingDetails"]], Report=list("D_xct"=Report$D_xct,"Epsilon1_sct"=Var_rot$Psi_rot,"Epsilon2_sct"=Var_rot$Psi_rot), PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateDir, Year_Set=Year_Set, Rotate=MapDetails_List[["Rotate"]], category_names=paste0("Factor_",1:length(category_names)), mar=c(0,0,2,0), oma=c(1.5,1.5,0,0), Cex=MapDetails_List[["Cex"]], cex=1.8, mfrow=Dim_year, cex.main=1.0, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], plot_legend_fig=FALSE)
      }  # SpatialDeltaGLMM::

      # Plot loadings
      png( file=paste0(DateDir,"Factor_loadings--",Par_name,".png"), width=Dim_species[2]*2.5, height=Dim_species[1]*2.5, units="in", res=200 )
        par( mfrow=Dim_species, mar=c(0,2,2,0) )
        for( cI in 1:length(Species_Set)) PlotLoadings( Var_rot$L_pj_rot, whichfactor=cI )
      dev.off()

      # Plot loadings
      if(Par_name %in% c("Epsilon1","Epsilon2")) Mat_sc = apply(Var_rot$Psi_rot, MARGIN=1:2, FUN=mean)
      if(Par_name %in% c("Omega1","Omega2")) Mat_sc = Var_rot$Psi_rot
      PlotMap_Fn( MappingDetails=MapDetails_List[["MappingDetails"]], Mat=Mat_sc, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=paste0(DateDir,"Factor_maps--",Par_name), Year_Set=paste0("Factor_",1:length(Species_Set)), Rotate=MapDetails_List[["Rotate"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(2.5,2.5,0,0), Cex=0.01, mfrow=Dim_species, pch=20, Legend=list(use=TRUE, x=c(65,75), y=c(35,65)), plot_legend_fig=FALSE)
    }  # SpatialDeltaGLMM:::
  }
}
