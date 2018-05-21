

#' Calculate the compositional-expansion
#'
#' \code{calculate_proportion} takes output from a VAST run and calculates the proportion of biomass in different categories
#'
#' @param Index output from \code{SpatialDeltaGLMM::PlotIndex_Fn}
#' @inheritParams SpatialDeltaGLMM::PlotIndex_Fn
#' @param ... list of settings to pass to \code{sdreport}
#'
#' @return Tagged list of output
#' \describe{
#'   \item{Prop_ctl}{Proportion of biomass for each category c, time t, and stratum l}
#'   \item{Neff_tl}{Effective sample size (median across categories) for each time t and stratum l}
#' }
#'
#' @export
calculate_proportion = function( TmbData, Index, Year_Set=NULL, Years2Include=NULL, strata_names=NULL, category_names=NULL, plot_legend=TRUE,
  DirName=paste0(getwd(),"/"), PlotName="Proportion.png", interval_width=1, width=6, height=6, ... ){

  # Warnings and errors
  if( !all(TmbData[['FieldConfig']] %in% c(-2,-1)) ){
    stop("Derivation only included for independent categories")
  }
  Index_ctl = array(Index$Index_ctl[,,,'Estimate'],dim=dim(Index$Index_ctl)[1:3])
  SE_Index_ctl = array(Index$Index_ctl[,,,'Std. Error'],dim=dim(Index$Index_ctl)[1:3])

  # Calculate proportions, and total biomass
  Prop_ctl = Index_ctl / outer(rep(1,TmbData$n_c),apply(Index_ctl,MARGIN=2:3,FUN=sum))
  Index_tl = apply(Index_ctl,MARGIN=2:3,FUN=sum)
  SE_Index_tl = sqrt(apply(SE_Index_ctl^2,MARGIN=2:3,FUN=sum,na.rm=TRUE))

  # Approximate variance for proportions, and effective sample size
  Neff_ctl = var_Prop_ctl = array(NA,dim=dim(Prop_ctl))
  for( cI in 1:dim(var_Prop_ctl)[1]){
  for( tI in 1:dim(var_Prop_ctl)[2]){
  for( lI in 1:dim(var_Prop_ctl)[3]){
    # Original version
    #var_Prop_ctl[cI,tI,lI] = Index_ctl[cI,tI,lI]^2/Index_tl[tI,lI]^2 * (SE_Index_ctl[cI,tI,lI]^2/Index_ctl[cI,tI,lI]^2  + SE_Index_tl[tI,lI]^2/Index_tl[tI,lI]^2 )
    # Slightly extended version
    var_Prop_ctl[cI,tI,lI] = Index_ctl[cI,tI,lI]^2/Index_tl[tI,lI]^2 * (SE_Index_ctl[cI,tI,lI]^2/Index_ctl[cI,tI,lI]^2 - 2*SE_Index_ctl[cI,tI,lI]^2/(Index_ctl[cI,tI,lI]*Index_tl[tI,lI]) + SE_Index_tl[tI,lI]^2/Index_tl[tI,lI]^2 )
    # Covert to effective sample size
    Neff_ctl[cI,tI,lI] = Prop_ctl[cI,tI,lI] * (1-Prop_ctl[cI,tI,lI]) / var_Prop_ctl[cI,tI,lI]
  }}}

  # Median effective sample size across categories
  Neff_tl = apply(Neff_ctl, MARGIN=2:3, FUN=median, na.rm=TRUE)

  # Fill in missing
  if( is.null(Year_Set) ) Year_Set = 1:TmbData$n_t
  if( is.null(Years2Include) ) Years2Include = 1:TmbData$n_t
  if( is.null(strata_names) ) strata_names = 1:TmbData$n_l
  if( is.null(category_names) ) category_names = 1:TmbData$n_c

  # Plot
  Par = list( mar=c(2,2,1,0), mgp=c(2,0.5,0), tck=-0.02, yaxs="i", oma=c(2,2,0,0), mfrow=c(ceiling(sqrt(TmbData$n_t)),ceiling(TmbData$n_t/ceiling(sqrt(TmbData$n_t)))), ... )
  png( file=paste0(DirName,"/",PlotName), width=width, height=height, res=200, units="in")
    par( Par )
    for( tI in 1:TmbData$n_t ){
      # Calculate y-axis limits
      Ylim = c(0, max(Prop_ctl[,tI,]%o%c(1,1) + sqrt(var_Prop_ctl[,tI,])%o%c(-interval_width,interval_width),na.rm=TRUE))
      # Plot stuff
      plot(1, type="n", xlim=range(category_names), ylim=1.05*Ylim, xlab="", ylab="", main=ifelse(TmbData$n_t>1,paste0("Year ",Year_Set[tI]),"") )
      for(l in 1:TmbData$n_l){
        SpatialDeltaGLMM:::Plot_Points_and_Bounds_Fn( y=Prop_ctl[,tI,l], x=1:TmbData$n_c+seq(-0.1,0.1,length=TmbData$n_l)[l], ybounds=Prop_ctl[,tI,]%o%c(1,1) + sqrt(var_Prop_ctl[,tI,])%o%c(-interval_width,interval_width), type="b", col=rainbow(TmbData[['n_l']])[l], col_bounds=rainbow(TmbData[['n_l']])[l], ylim=Ylim)
      }
      if(plot_legend==TRUE) legend( "top", bty="n", fill=rainbow(TmbData[['n_l']]), legend=as.character(strata_names), ncol=2 )
    }
    mtext( side=1:2, text=c("Age","Proportion of biomass"), outer=TRUE, line=c(0,0) )
  dev.off()

  # Return stuff
  Return = list("Prop_ctl"=Prop_ctl, "Neff_tl"=Neff_tl, "var_Prop_ctl"=var_Prop_ctl, "Index_tl"=Index_tl, "Neff_ctl"=Neff_ctl)
  return( Return )
}
