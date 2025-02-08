

#' Calculate the compositional-expansion
#'
#' \code{calculate_proportion} takes output from a VAST run and calculates the proportion of biomass in different categories
#'
#' @inheritParams plot_biomass_index
#' @inheritParams make_data
#' @param sample_size_method Method used to calculate the variance in proportions, which is then converted to an approximately equivalent multinomial sample size that can be used as input-sample-size in a subsequent stock assessment model.  Options are:
#' \describe{
#'   \item{\code{sample_size_method="Taylor_series"}}{a Taylor-series approximation to the ratio of X/(X+Y) where X is the category-specific index and Y is the index for for all other categories}
#'   \item{\code{sample_size_method="sample_based"}}{Taking samples from the joint precision of fixed and random effects, calculating proportions for each sample, and then computing the variance across those samples}
#' }
#' The sample-based approximation is expected to have higher variance and therefore lower approximate sample size.  However, it may also have poor performance in cases when variance estimates are imprecise (such that the multivariate-normal approximation to joint precision is poor), and has not been thoroughly groundtested in real-world cases.
#' @param Index output from \code{FishStatsUtils::plot_biomass_index}
#' @param ... list of arguments to pass to \code{plot_index}
#'
#' @return Tagged list of output
#' \describe{
#'   \item{Prop_ctl}{Proportion of biomass for each category c, time t, and stratum l}
#'   \item{Neff_tl}{Effective sample size (median across categories) for each time t and stratum l}
#' }
#'
#' @references For details regarding multivariate index standardization and expansion see \url{https://cdnsciencepub.com/doi/full/10.1139/cjfas-2018-0015}
#' @export
calculate_proportion <-
function( fit,
          TmbData,
          Index,
          Expansion_cz = NULL,
          year_labels = NULL,
          years_to_plot = NULL,
          strata_names = NULL,
          category_names = NULL,
          sample_size_method = c("Taylor_series","sample_based"),
          plot_legend = ifelse(TmbData$n_l>1,TRUE,FALSE),
          DirName = getwd(),
          PlotName = "Proportion.png",
          PlotName2 = "Average.png",
          interval_width = 1,
          width = 6,
          height = 6,
          xlab = "Category",
          ylab = "Proportion",
          n_samples = 250,
          ... ){

  # Warnings and errors
  if( !all(TmbData[['FieldConfig']] %in% c(-3,-2,-1)) ){
    message("Derivation only included for independent categories")
    return( invisible("Not run") )
  }
  Index_ctl = array(Index$Index_ctl[,,,'Estimate'],dim=dim(Index$Index_ctl)[1:3])
  SE_Index_ctl = array(Index$Index_ctl[,,,'Std. Error'],dim=dim(Index$Index_ctl)[1:3])

  if( !is.null(Expansion_cz) ){
    Index_ctl = as.array(Index_ctl[Expansion_cz[,1]==1,,,drop=FALSE])
    SE_Index_ctl = as.array(SE_Index_ctl[Expansion_cz[,1]==1,,,drop=FALSE])
    category_names = category_names[Expansion_cz[,1]==1]
  }

  # Calculate proportions, and total biomass
  Prop_ctl = Index_ctl / outer(rep(1,dim(Index_ctl)[1]),apply(Index_ctl,MARGIN=2:3,FUN=sum))
  Index_tl = apply(Index_ctl,MARGIN=2:3,FUN=sum)
  SE_Index_tl = sqrt(apply(SE_Index_ctl^2,MARGIN=2:3,FUN=sum,na.rm=TRUE))

  if( tolower(sample_size_method[1]) == "taylor_series" ){
    # Approximate variance for proportions, and effective sample size
    Neff_ctl = var_Prop_ctl = array(NA,dim=dim(Prop_ctl))
    for( cI in 1:dim(var_Prop_ctl)[1]){
    for( tI in 1:dim(var_Prop_ctl)[2]){
    for( lI in 1:dim(var_Prop_ctl)[3]){
      # Original version
      #var_Prop_ctl[cI,tI,lI] = Index_ctl[cI,tI,lI]^2/Index_tl[tI,lI]^2 * (SE_Index_ctl[cI,tI,lI]^2/Index_ctl[cI,tI,lI]^2  + SE_Index_tl[tI,lI]^2/Index_tl[tI,lI]^2 )
      # Slightly extended version
      var_Prop_ctl[cI,tI,lI] = Index_ctl[cI,tI,lI]^2/Index_tl[tI,lI]^2 * (SE_Index_ctl[cI,tI,lI]^2/Index_ctl[cI,tI,lI]^2 - 2*SE_Index_ctl[cI,tI,lI]^2/(Index_ctl[cI,tI,lI]*Index_tl[tI,lI]) + SE_Index_tl[tI,lI]^2/Index_tl[tI,lI]^2 )
      var_Prop_ctl[cI,tI,lI] = ifelse( Index_ctl[cI,tI,lI]==0, 0, var_Prop_ctl[cI,tI,lI] )  # If dividing by zero, replace with 0
      # Covert to effective sample size
      #Neff_ctl[cI,tI,lI] = Prop_ctl[cI,tI,lI] * (1-Prop_ctl[cI,tI,lI]) / var_Prop_ctl[cI,tI,lI]
    }}}
  }else{
    Index_ctlr = sample_variable( Sdreport = fit$parameter_estimates$SD, 
                            Obj = fit$tmb_list$Obj, 
                            variable_name = "Index_ctl", 
                            n_samples = n_samples, 
                            sample_fixed = TRUE )
    Prop_ctlr = Index_ctlr / outer( rep(1,dim(Index_ctlr)[1]), apply(Index_ctlr,MARGIN=2:4,FUN=sum) )
    var_Prop_ctl = apply(Prop_ctlr, MARGIN=1:3, FUN=var)
  }

  # Covert to effective sample size
  Neff_ctl = Prop_ctl * (1-Prop_ctl) / var_Prop_ctl

  # Median effective sample size across categories
  Neff_tl = apply(Neff_ctl, MARGIN=2:3, FUN=median, na.rm=TRUE)

  # Plot
  if( !is.na(PlotName) ){
    plot_index( Index_ctl=Prop_ctl, sd_Index_ctl=sqrt(var_Prop_ctl), year_labels=year_labels, years_to_plot=years_to_plot,
      strata_names=strata_names, category_names=category_names, plot_legend=plot_legend,
      DirName=DirName, PlotName=PlotName, interval_width=interval_width, width=width, height=height,
      xlab=xlab, ylab=ylab, scale="uniform", ... )
  }

  # Calculate weighted mean
  sd_Mean_tl = Mean_tl = apply( Prop_ctl, MARGIN=2:3, FUN=function(vec){sum(vec*(1:length(vec)))} )
  for( tI in 1:nrow(sd_Mean_tl)){
  for( lI in 1:ncol(sd_Mean_tl)){
    sd_Mean_tl[tI,lI] = sqrt(sum( var_Prop_ctl[,tI,lI] * (1:dim(var_Prop_ctl)[1] - Mean_tl[tI,lI])^2 ))
  }}

  # Plot
  if( !is.na(PlotName2) ){
    plot_index( Index_ctl=1%o%Mean_tl, sd_Index_ctl=1%o%sd_Mean_tl, year_labels=year_labels, years_to_plot=years_to_plot,
      strata_names=strata_names, category_names=category_names, plot_legend=plot_legend,
      DirName=DirName, PlotName=PlotName2, interval_width=interval_width, width=width, height=height,
      xlab=xlab, ylab="Category", scale="uniform", Yrange=c(NA,NA), ... )     # , Yrange=c(1,dim(var_Prop_ctl)[1])
  }

  # Return stuff
  Return = list("Prop_ctl"=Prop_ctl, "Neff_tl"=Neff_tl, "var_Prop_ctl"=var_Prop_ctl, "Index_tl"=Index_tl, "Neff_ctl"=Neff_ctl,
    "Mean_tl"=Mean_tl, "sd_Mean_tl"=sd_Mean_tl )
  return( invisible(Return) )
}
