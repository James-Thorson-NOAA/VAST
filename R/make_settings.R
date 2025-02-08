
#' Make list of settings
#'
#' \code{make_settings} makes a list of settings for a given purpose
#'
#' This function assembles a default set of user-decisions for a specified modelling purpose. The default settings are guessed based on generic guidance, and should be carefully reviewed for real-world purposes. If the user supplies values for individual settings e.g. \code{FieldConfig}, then these values override the defaults that are provided by interpreting \code{purpose}
#'
#' @inheritParams make_data
#' @inheritParams make_extrapolation_info
#' @inheritParams make_spatial_info
#' @inheritParams Convert_LL_to_UTM_Fn
#' @inheritParams plot_biomass_index
#' @inheritParams fit_tmb
#' @param purpose character indicating what purpose is intended for the model, and therefore what default settings are perhaps appropriate. Many of these have examples at the VAST wiki https://github.com/James-Thorson-NOAA/VAST/wiki. Only currently implemented for:
#' \describe{
#' \item{\code{purpose="index"}}{Index of abundance calculated summing across space to get annual abundance values for each category}
#' \item{\code{purpose="index2"}}{The same as "index" except uses gamma distribution for positive catches, restricts \code{max_cells} to 2000, and uses bias correction. This is currently recommended over "index".}
#' \item{\code{purpose="condition_and_density"}}{Jointly estimate density (numbers per area) and fish condition (relative weight given length, used to predict density-weighted average condition}
#' \item{\code{purpose="MICE"}}{Model of intermediate complexity to estimate species interactions spatially}
#' \item{\code{purpose="ordination"}}{Multivariate ordination to identify similar species by estimating a reduced set of axes that collectively explain variability}
#' \item{\code{purpose="EOF"}}{An empirical orthogonal function analysis to ordinate on years instead of categories as in "ordination". Deprecated, use "EOF2" instead.}
#' \item{\code{purpose="EOF2"}}{Same as "EOF" but uses improved settings that match updates to the package.}
#' \item{\code{purpose="EOF3"}}{Same as "EOF2" but ensures that spatio-temporal factors are zero-centered, such that estimated Omega represents distribution in the "average" year.}
#' }
#' @param use_anisotropy Boolean indicating whether to estimate two additional parameters representing geometric anisotropy
#' @param vars_to_correct a character-vector listing which parameters to include for bias-correction, as passed to \code{\link{fit_tmb}}
#' @param treat_nonencounter_as_zero Boolean indicating whether to treat any year-category combination as having zero biomass when generating abundance indices and resulting compositional estimates
#' @param n_categories number of categories in a multivariate model (only necessary to specify given some values for \code{purpose})
#'
#' @return Tagged list containing default settings for a given purpose, use \code{names} on output to see or modify list of settings.
#'
#' @family wrapper functions
#' @seealso \code{\link[VAST]{VAST}} for general documentation, \code{\link{make_settings}} for generic settings, \code{\link{fit_model}} for model fitting, and \code{\link{plot_results}} for generic plots
#'
#' @references For discussion of some of these options see \url{https://doi.org/10.1016/j.fishres.2018.10.013}
#' @export
make_settings <-
function( n_x,
          purpose = "index",
          Region,
          fine_scale = TRUE,
          strata.limits = data.frame('STRATA'="All_areas"),
          zone=NA,
          FieldConfig,
          RhoConfig,
          OverdispersionConfig,
          ObsModel,
          bias.correct,
          Options,
          use_anisotropy,
          vars_to_correct,
          Version,
          treat_nonencounter_as_zero,
          n_categories,
          VamConfig,
          max_cells,
          knot_method,
          mesh_package ){

  # Get version
  if(missing(Version)) Version = get_latest_version()
  purpose_found = FALSE

  ###################
  # Index standardization
  ###################

  # Deprecated
  if( tolower(purpose) == "index" ){
    purpose_found = TRUE
    warning( "The package author recommends using purpose=`index2` for updated defaults; purpose=`index` is retained for backwards compatibility but not recommended" )
    if( convert_version_name(Version) >= convert_version_name("VAST_v7_0_0") ){
      if(missing(FieldConfig)) FieldConfig = matrix( "IID", ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) )
    }else{
      if(missing(FieldConfig)) FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
    }
    if(missing(RhoConfig)) RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
    if(missing(VamConfig)) VamConfig = c("Method"=0, "Rank"=0, "Timing"=0)
    if(missing(OverdispersionConfig)) OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
    if(missing(ObsModel)) ObsModel = c(1,1)
    if(missing(bias.correct)) bias.correct = TRUE
    if(missing(treat_nonencounter_as_zero)) treat_nonencounter_as_zero = TRUE
    if(missing(Options)) Options =  c( "Calculate_Range"=TRUE, "Calculate_effective_area"=TRUE, "treat_nonencounter_as_zero"=treat_nonencounter_as_zero )
    if(missing(vars_to_correct)) vars_to_correct = c( "Index_cyl", "Index_ctl" )
    if(missing(knot_method)) knot_method = "samples"
    if(missing(max_cells)) max_cells = Inf
    if(missing(mesh_package)) mesh_package = "INLA"
  }

  # Current
  if( tolower(purpose) == "index2" ){
    purpose_found = TRUE
    if( convert_version_name(Version) >= convert_version_name("VAST_v7_0_0") ){
      if(missing(FieldConfig)) FieldConfig = matrix( "IID", ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) )
    }else{
      if(missing(FieldConfig)) FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
    }
    if(missing(RhoConfig)) RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
    if(missing(VamConfig)) VamConfig = c("Method"=0, "Rank"=0, "Timing"=0)
    if(missing(OverdispersionConfig)) OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
    if(missing(ObsModel)) ObsModel = c(2,1)
    if(missing(bias.correct)) bias.correct = TRUE
    if(missing(treat_nonencounter_as_zero)) treat_nonencounter_as_zero = TRUE
    if(missing(Options)) Options =  c( "Calculate_Range"=TRUE, "Calculate_effective_area"=TRUE, "treat_nonencounter_as_zero"=treat_nonencounter_as_zero )
    if(missing(vars_to_correct)) vars_to_correct = c( "Index_cyl", "Index_ctl" )
    if(missing(knot_method)) knot_method = "grid"
    if(missing(max_cells)) max_cells = max( 2000, n_x*10 )
    if(missing(mesh_package)) mesh_package = "INLA"
  }
  if( tolower(purpose) == "index3" ){
    purpose_found = TRUE
    if( convert_version_name(Version) >= convert_version_name("VAST_v7_0_0") ){
      if(missing(FieldConfig)) FieldConfig = matrix( "IID", ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) )
    }else{
      if(missing(FieldConfig)) FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
    }
    if(missing(RhoConfig)) RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
    if(missing(VamConfig)) VamConfig = c("Method"=0, "Rank"=0, "Timing"=0)
    if(missing(OverdispersionConfig)) OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
    if(missing(ObsModel)) ObsModel = c(2,1)
    if(missing(bias.correct)) bias.correct = TRUE
    if(missing(treat_nonencounter_as_zero)) treat_nonencounter_as_zero = TRUE
    if(missing(Options)) Options =  c( "Calculate_Range"=TRUE, "Calculate_effective_area"=TRUE, "treat_nonencounter_as_zero"=treat_nonencounter_as_zero )
    if(missing(vars_to_correct)) vars_to_correct = c( "Index_cyl", "Index_ctl" )
    if(missing(knot_method)) knot_method = "grid"
    if(missing(max_cells)) max_cells = max( 2000, n_x*10 )
    if(missing(mesh_package)) mesh_package = "fmesher"
  }

  ###################
  # Condition and density
  ###################
  if( tolower(purpose) == "condition_and_density" ){
    purpose_found = TRUE
    if( convert_version_name(Version) >= convert_version_name("VAST_v7_0_0") ){
      if(missing(FieldConfig)) FieldConfig = matrix( c(2,2,"IID",0,0,"IID"), ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) )
    }else{
      if(missing(FieldConfig)) FieldConfig = c("Omega1"=2, "Epsilon1"=2, "Omega2"=0, "Epsilon2"=0)
    }
    if(missing(RhoConfig)) RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
    if(missing(VamConfig)) VamConfig = c("Method"=0, "Rank"=0, "Timing"=0)
    if(missing(OverdispersionConfig)) OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
    if(missing(ObsModel)) ObsModel = c(1,4)
    if(missing(bias.correct)) bias.correct = TRUE
    if(missing(treat_nonencounter_as_zero)) treat_nonencounter_as_zero = FALSE
    if(missing(Options)) Options =  c( "Calculate_Range"=FALSE, "Calculate_effective_area"=FALSE, "Calculate_Cov_SE"=TRUE,
      "treat_nonencounter_as_zero"=treat_nonencounter_as_zero, "Project_factors"=TRUE )
    if(missing(vars_to_correct)) vars_to_correct = c( "Index_cyl", "Index_ctl" )
    if(missing(knot_method)) knot_method = "samples"
    if(missing(max_cells)) max_cells = Inf
    if(missing(mesh_package)) mesh_package = "INLA"
  }

  ###################
  # Spatial model of intermediate complexity for ecosystems (MICE-in-space)
  ###################
  if( tolower(purpose) %in% c("mice","interactions") ){
    purpose_found = TRUE
    if(missing(n_categories)) stop("Must supply `n_categories` when using purpose==`mice`")
    if( convert_version_name(Version) >= convert_version_name("VAST_v7_0_0") ){
      if(missing(FieldConfig)) FieldConfig = matrix( c(n_categories,n_categories,"IID",n_categories,n_categories,"IID"), ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) )
    }else{
      if(missing(FieldConfig)) FieldConfig = c("Omega1"=n_categories, "Epsilon1"=n_categories, "Omega2"=n_categories, "Epsilon2"=n_categories)
    }
    if(missing(RhoConfig)) RhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=6)
    if(missing(VamConfig)) VamConfig = c("Method"=2, "Rank"=n_categories, "Timing"=1)
    if(missing(OverdispersionConfig)) OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
    if(missing(ObsModel)) ObsModel = c(1,1)
    if(missing(bias.correct)) bias.correct = TRUE
    if(missing(Options)) Options =  c( "Calculate_Range"=FALSE, "Calculate_effective_area"=FALSE, "Calculate_Cov_SE"=FALSE,
      "Calculate_Fratio"=TRUE, "Estimate_B0"=TRUE )
    if(missing(vars_to_correct)) vars_to_correct = c( "Index_cyl", "Bratio_cyl", "Index_ctl", "Bratio_ctl" )
    if(missing(knot_method)) knot_method = "samples"
    if(missing(max_cells)) max_cells = Inf
    if(missing(mesh_package)) mesh_package = "INLA"
  }

  ###################
  # Spatial model for ordinating species
  ###################
  if( tolower(purpose) %in% c("ordination") ){
    purpose_found = TRUE
    if(missing(n_categories)) stop("Must supply `n_categories` when using purpose==`ordination`")
    if( convert_version_name(Version) >= convert_version_name("VAST_v7_0_0") ){
      if(missing(FieldConfig)) FieldConfig = matrix( c(n_categories,n_categories,n_categories,n_categories,n_categories,n_categories), ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) )
    }else{
      if(missing(FieldConfig)) FieldConfig = c("Omega1"=n_categories, "Epsilon1"=n_categories, "Omega2"=n_categories, "Epsilon2"=n_categories)
    }
    if(missing(RhoConfig)) RhoConfig = c("Beta1"=4, "Beta2"=4, "Epsilon1"=4, "Epsilon2"=4)
    if(missing(VamConfig)) VamConfig = c("Method"=0, "Rank"=0, "Timing"=0)
    if(missing(OverdispersionConfig)) OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
    if(missing(ObsModel)) ObsModel = c(1,1)
    if(missing(bias.correct)) bias.correct = FALSE
    if(missing(Options)) Options =  c( "Calculate_Range"=FALSE, "Calculate_effective_area"=FALSE, "Calculate_Cov_SE"=TRUE, "Project_factors"=TRUE )
    if(missing(vars_to_correct)) vars_to_correct = c( "Index_cyl", "Index_ctl" )
    if(missing(knot_method)) knot_method = "samples"
    if(missing(max_cells)) max_cells = Inf
    if(missing(mesh_package)) mesh_package = "INLA"
  }

  ###################
  # Spatial EOF analysis
  ###################

  if( tolower(purpose) %in% c("eof") ){
    purpose_found = TRUE
    warning( "The package author recommends using purpose=`eof3` for improved model specification; purpose=`eof` is retained for backwards compatibility but not recommended" )
    if(missing(n_categories)) stop("Must supply `n_categories` when using purpose==`eof`")
    if( convert_version_name(Version) >= convert_version_name("VAST_v7_0_0") ){
      if(missing(FieldConfig)) FieldConfig = matrix( c(0,n_categories,"IID", 0,0,"IID"), ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) )
    }else{
      if(missing(FieldConfig)) FieldConfig = c("Omega1"=0, "Epsilon1"=n_categories, "Omega2"=0, "Epsilon2"=0)
    }
    if(missing(RhoConfig)) RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
    if(missing(VamConfig)) VamConfig = c("Method"=0, "Rank"=0, "Timing"=0)
    if(missing(OverdispersionConfig)) OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
    if(missing(ObsModel)) ObsModel = c(1,1)
    if(missing(bias.correct)) bias.correct = FALSE
    if(missing(Options)) Options =  c( "Calculate_Range"=FALSE, "Calculate_effective_area"=FALSE, "Calculate_Cov_SE"=TRUE, "Project_factors"=TRUE )
    if(missing(vars_to_correct)) vars_to_correct = c( "Index_cyl", "Index_ctl" )
    if(missing(knot_method)) knot_method = "samples"
    if(missing(max_cells)) max_cells = Inf
    if(missing(mesh_package)) mesh_package = "INLA"
  }

  if( tolower(purpose) %in% c("eof2") ){
    purpose_found = TRUE
    warning( "The package author recommends using purpose=`eof3` for improved model specification; purpose=`eof2` is retained for backwards compatibility but not recommended" )
    if( convert_version_name(Version) < convert_version_name("VAST_v9_2_0") ){
      stop("Cannot use purpose=`eof2` with VAST versions < 9.2.0")
    }
    if(missing(FieldConfig)) FieldConfig = matrix( c("IID","Identity","IID",2, 0,0,"IID","Identity"), ncol=2, nrow=4, dimnames=list(c("Omega","Epsilon","Beta","Epsilon_year"),c("Component_1","Component_2")) )
    if(missing(RhoConfig)) RhoConfig = c("Beta1"=0, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
    if(missing(VamConfig)) VamConfig = c("Method"=0, "Rank"=0, "Timing"=0)
    if(missing(OverdispersionConfig)) OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
    if(missing(ObsModel)) ObsModel = c(2,1)
    if(missing(bias.correct)) bias.correct = FALSE
    if(missing(Options)) Options =  c( "Calculate_Range"=FALSE, "Calculate_effective_area"=FALSE, "Calculate_Cov_SE"=TRUE, "Project_factors"=TRUE )
    if(missing(vars_to_correct)) vars_to_correct = c( "Index_cyl", "Index_ctl" )
    if(missing(knot_method)) knot_method = "grid"
    if(missing(max_cells)) max_cells = max( 2000, n_x*10 )
    if(missing(mesh_package)) mesh_package = "INLA"
  }

  if( tolower(purpose) %in% c("eof3") ){
    purpose_found = TRUE
    if( convert_version_name(Version) < convert_version_name("VAST_v13_0_0") ){
      stop("Cannot use purpose=`eof2` with VAST versions < 13.0.0")
    }
    if(missing(FieldConfig)) FieldConfig = matrix( c("IID","Identity","IID",2, 0,0,"IID","Identity"), ncol=2, nrow=4, dimnames=list(c("Omega","Epsilon","Beta","Epsilon_year"),c("Component_1","Component_2")) )
    if(missing(RhoConfig)) RhoConfig = c("Beta1"=0, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
    if(missing(VamConfig)) VamConfig = c("Method"=0, "Rank"=0, "Timing"=0)
    if(missing(OverdispersionConfig)) OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
    if(missing(ObsModel)) ObsModel = c(2,1)
    if(missing(bias.correct)) bias.correct = FALSE
    if(missing(Options)) Options =  c( "Calculate_Range"=FALSE, "Calculate_effective_area"=FALSE, "Calculate_Cov_SE"=TRUE, "Project_factors"=TRUE, "zerosum_penalty"=1 )
    if(missing(vars_to_correct)) vars_to_correct = c( "Index_cyl", "Index_ctl" )
    if(missing(knot_method)) knot_method = "grid"
    if(missing(max_cells)) max_cells = max( 2000, n_x*10 )
    if(missing(mesh_package)) mesh_package = "INLA"
  }

  ###################
  # Other settings and formatting
  ###################

  # Check for bad input
  if( purpose_found==FALSE ){
    stop("'purpose' is currently set up only for index-standardization models, correlations between condition and density, MICE-in-space models, ordination, or empirical-orthogonal-function analysis")
  }

  # Other defaults
  grid_size_km = 25
  Method = "Mesh"
  if(missing(use_anisotropy)) use_anisotropy = TRUE

  # Default naming
  names(RhoConfig) = c("Beta1","Beta2","Epsilon1","Epsilon2")

  # Bundle and export
  settings = list("Version" = Version,
             "n_x" = n_x,
             "Region" = Region,
             "strata.limits" = strata.limits,
             "zone" = zone,
             "FieldConfig" = FieldConfig,
             "RhoConfig" = RhoConfig,
             "VamConfig" = VamConfig,
             "OverdispersionConfig" = OverdispersionConfig,
             "ObsModel" = ObsModel,
             "vars_to_correct" = vars_to_correct,
             "Options" = Options,
             "grid_size_km" = grid_size_km,
             "max_cells" = max_cells,
             "knot_method" = knot_method,
             "Method" = Method,
             "use_anisotropy" = use_anisotropy,
             "fine_scale" = fine_scale,
             "bias.correct" = bias.correct,
             "mesh_package" = mesh_package )
  return(settings)
}
