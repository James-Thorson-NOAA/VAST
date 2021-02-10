
#' Build data input for VAST model
#'
#' \code{make_data} builds a tagged list of data inputs used by TMB for running the model
#'
#' Specification of \code{FieldConfig} can be seen by calling \code{\link[FishStatsUtils]{make_settings}},
#'   which is the recommended way of generating this input for beginning users.
#'
#' Argument \code{FieldConfig} is a matrix of form
#'   \code{FieldConfig = matrix( c(0,10,"IID","Identity", "AR1",10,"IID","Identity"), ncol=2, nrow=4, dimnames=list(c("Omega","Epsilon","Beta","Epsilon_year"),c("Component_1","Component_2"))}.
#'   However, for backwards compatibility, \code{FieldConfig} can instead be specified as
#'   a vector of format \code{FieldConfig = c("Omega1"=0, "Epsilon1"=10, "Omega2"="AR1", "Epsilon2"=10)},
#'   which generates the same settings as the matrix specification, given the values of each shown in this example.
#'
#' Both vector (simplified) and matrix (full) specification of \code{FieldConfig} involve named elements
#'   using the following naming conventions:
#' \describe{
#'    \item{Omega}{specifies whether spatial variation is present and/or correlated among variables}
#'    \item{Epsilon}{specifies whether spatio-temporal variation is present and/or correlated among variables}
#'    \item{Beta}{specifies whether temporal variation (a.k.a. "intercepts") is present and/or correlated among variables}
#'    \item{Epsilon_year}{specifies whether spatio-temporal variation is correlated among years}
#' }
#' Meanwhile, \code{Component_1} (or the numeral "1" after a component name) refers to the 1st lienar predictor (e.g., of a delta-model),
#'   while \code{Component_2} (or the numeral "2" after a component name) refers to the 2nd linear predictor.
#'   The simplified vector-specification does not include slots for \code{Beta} or \code{Epsilon_year} and therefore
#'   is not as general.
#'
#' In each slot of \code{FieldConfig}, the user can specify various options:
#' \describe{
#'   \item{0}{turns off a given model component}
#'   \item{integer greater than zero}{specifies the rank (number of factors) in a factor-analysis covariance matrix}
#'   \item{"AR1"}{specifies that a model component is correlated following an first-order autoregressive process}
#'   \item{"IID"}{specifies that a given model component is a random effect that is independent for every level}
#'   \item{"Identity"}{specifies that a given model component has covariance of an identity-matrix;
#'     this is only useful for \code{Epsilon_year} to "turn off" covariance among years while still including spatio-temporal variation}
#' }
#'
#' \code{make_data} generates arrays of covariates \code{X_gtp} and \code{X_itp} using \code{\link[FishStatsUtils]{make_covariates}};
#' see that function for more details.
#'
#' @inheritParams FishStatsUtils::make_covariates
#' @param b_i Numeric vector, providing sampled value (biomass, counts, etc.) for each observation i
#' @param a_i Numeric vector containing values greater than zero, providing sampled area for each
#'        observation i;  use \code{a_i=1} for observations without a natural area measurement, while
#'        noting that resulting densities no longer have interpretable units in that case)
#' @param c_iz Vector of integers ranging from 0 to the number of variables minus 1, providing the
#'        category (e.g., species, length-bin) for each observation i
#' @param t_i Vector of integers, providing the time (e.g., calendar year) for each observation i
#' @param e_i Optional vector of integers ranging from 0 to the number of different error distributions,
#'        providing the error distribution to use for each observation i;
#'        by default \code{e_i=c_iz[,1]} such that each category has a unique estimated magnitude of measurement error
#' @param v_i Vector of integers ranging from 0 to the number of vessels minus 1,
#'        providing sampling category (e.g., vessel or tow) associated with overdispersed variation for each observation i
#'        (by default \code{v_i=0} for all samples, which will not affect things given the default values for \code{OverdispersionConfig})
#' @param Version Which CPP version to use.  If missing, defaults to latest version using \code{\link[FishStatsUtils]{get_latest_version(package="VAST")}}.
#'        Can be used to specify using an older CPP, to maintain backwards compatibility.
#' @param FieldConfig See Details section of \code{\link[VAST]{make_data}} for details
#' @param OverdispersionConfig a vector of format \code{c("eta1"=0, "eta2"="AR1")} governing any correlated overdispersion
#'        among categories for each level of \code{v_i}, where eta1 is for encounter probability, and eta2 is for positive catch rates,
#'        where \code{0} is off, code{"AR1"} is an AR1 process, and aninteger greater than zero (e.g., \code{2}) is the number of elements in a factor-analysis covariance (by default,
#'        \code{c("eta1"=0, "eta2"=0)} and this turns off overdispersion)
#' @param ObsModel_ez an optional matrix with two columns where the first column specifies the distribution for positive catch rates,
#'        and the second column specifies the functional form for encounter probabilities
#' \describe{
#'   \item{\code{ObsModel_ez[e,1]=0}}{Normal}
#'   \item{\code{ObsModel_ez[e,1]=1}}{Lognormal, using bias-corrected-mean and log-sd as parameters; I recommend using \code{ObsModel_ez[e,1]=4}
#'        instead for consistent interpretation of logSigmaM as log-CV as used for Gamma and inverse-Gaussian distribution}
#'   \item{\code{ObsModel_ez[e,1]=2}}{Gamma}
#'   \item{\code{ObsModel_ez[e,1]=4}}{Lornormal, using bias-corrected-mean and log-coefficient of variation (CV) as parameters;  see note for \code{ObsModel_ez[e,1]=1}}
#'   \item{\code{ObsModel_ez[e,1]=5}}{Negative binomial}
#'   \item{\code{ObsModel_ez[e,1]=7}}{Poisson (more numerically stable than negative-binomial)}
#'   \item{\code{ObsModel_ez[e,1]=10}}{Tweedie distribution, where expected biomass (lambda) is the product of 1st-component and 2nd-component,
#'          variance scalar (phi) is the 1st component, and logis-SigmaM is the power. This parameterization is fast (i.e., comparable to delta-models)
#'          as long as random effects are turned off for the 1st component, but is otherwise extremely slow.}
#'   \item{\code{ObsModel_ez[e,1]=11}}{Zero-inflated Poisson with additional normally-distributed variation overdispersion in the
#'        log-intensity of the Poisson distribution}
#'   \item{\code{ObsModel_ez[e,1]=12}}{Poisson distribution (not zero-inflated) with log-intensity from the 1st linear predictor,
#'        to be used in combination with the Poisson-link delta model for combining multiple data types}
#'   \item{\code{ObsModel_ez[e,1]=13}}{Bernoulli distribution using complementary log-log (cloglog) link from the 1st linear predictor,
#'        to be used in combination with the Poisson-link delta model for combining multiple data types}
#'   \item{\code{ObsModel_ez[e,1]=14}}{Similar to 12, but also including lognormal overdispersion}
#'   \item{\code{ObsModel_ez[e,2]=0}}{Conventional delta-model using logit-link for encounter probability and log-link for positive catch rates}
#'   \item{\code{ObsModel_ez[e,2]=1}}{Alternative "Poisson-link delta-model" using log-link for numbers-density and log-link for biomass per number}
#'   \item{\code{ObsModel_ez[e,2]=2}}{Link function for Tweedie distribution, necessary for \code{ObsModel_ez[e,1]=8} or \code{ObsModel_ez[e,1]=10}}
#'   \item{\code{ObsModel_ez[e,2]=3}}{Conventional delta-model, but fixing encounter probability=1 for any year where all samples encounter the species}
#'   \item{\code{ObsModel_ez[e,2]=4}}{Poisson-link delta-model, but fixing encounter probability=1 for any year where all samples encounter the
#'        species and encounter probability=0 for any year where no samples encounter the species}
#' }
#' @param RhoConfig vector of form \code{c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0)} specifying whether either intercepts (Beta1 and Beta2)
#'        or spatio-temporal variation (Epsilon1 and Epsilon2) is structured among time intervals, e.g.
#'        for component \code{Epsilon2} indicated in the 4rd slot:
#' \describe{
#'   \item{\code{RhoConfig[4]=0}}{Each year as fixed effect}
#'   \item{\code{RhoConfig[4]=1}}{Each year as an independent and identically distributed random effect, thus estimating the variance as fixed effect}
#'   \item{\code{RhoConfig[4]=2}}{Each year as a random effect following a random walk, thus estimating the variance as fixed effect}
#'   \item{\code{RhoConfig[4]=3}}{Constant among years as fixed effect}
#'   \item{\code{RhoConfig[4]=4}}{Each year as a random effect following a first-order autoregressive process, thus estimating the variance as fixed effects and a single first-order autoregression parameter}
#'   \item{\code{RhoConfig[4]=5}}{Each year as a random effect following a first-order autoregressive process, estimating the variance as fixed effects and a separate first-order autoregression parameter for each category}
#'   \item{\code{RhoConfig[4]=6}}{Only possible for \code{Epsilon2} or \code{Beta2}, and specifying that that associated hyperparameters parameters have identical values to the first component \code{Epsilon1} or \code{Beta1}}
#' }
#' If missing, the default is to assume a value of zero for each element (i.e., \code{RhoConfig[1:4]=0})
#' @param VamConfig Options to estimate interactions, containing three slots:
#' \describe{
#'   \item{\code{VamConfig[0]}}{selects method for forming interaction matrix; Turn off feature using 0, or I recommend using 2 by default}
#'   \item{\code{VamConfig[1]}}{indicates the rank of the interaction matrix, indicating the number of community axes that have regulated dynamics}
#'   \item{\code{VamConfig[2]}}{Indicates whether interactions occur before spatio-temporal variation (\code{VamConfig[2]=0}) or after \code{VamConfig[2]=1}}
#' }
#' @param X1_formula right-sided formula affecting the 1st linear predictor for density, e.g.,
#'        \code{X1_formula=~BOT_DEPTH+BOT_DEPTH^2} for a quadratic effect of variable \code{BOT_DEPTH}
#' @param X2_formula same as \code{X1_formula} but affecting the 2nd linear predictor for density
#' @param covariate_data data-frame of covariates for use when specifying \code{X1_formula} and \code{X2_formula}
#' @param X1config_cp matrix of settings for each density covariate for the 1st lienar predictor,
#'        where the row corresponds to model category, and column corresponds to each density covariate
#' \describe{
#'   \item{\code{X1config_cp[c,p]=0}}{\code{X_ip[,p]} has no effect on the 1st linear predictor for category c}
#'   \item{\code{X1config_cp[c,p]=1}}{\code{X_ip[,p]} has a linear effect on 1st linear predictor for category c}
#'   \item{\code{X1config_cp[c,p]=2}}{\code{X_ip[,p]} has a spatially varying, zero-centered linear effect on 1st linear predictor for category c}
#'   \item{\code{X1config_cp[c,p]=3}}{\code{X_ip[,p]} has a spatially varying linear effect on 1st linear predictor for category c}
#'   \item{\code{X1config_cp[c,p]=-1}}{\code{X1config_cp[c,p]=-1} is the same as \code{X1config_cp[c,p]=2}, but without including the log-likelihood term; this is useful in special cases when carefully mirroring spatially varying coefficients, e.g., to use cohort effects in a age-structured spatio-temporal model}
#' }
#' @param X2config_cp Same as argument \code{X1config_cp} but for 2nd linear predictor
#' @param Q1_formula same as \code{X1_formula} but affecting the 1st linear predictor for catchability (a.k.a. detectability)
#' @param Q2_formula same as \code{X1_formula} but affecting the 2nd linear predictor for catchability (a.k.a. detectability)
#' @param catchability_data data-frame of covariates for use when specifying \code{Q1_formula} and \code{Q2_formula}
#' @param Q1config_k Same as argument \code{X1config_cp} but affecting affecting the 1st linear predictor for catchability,
#'        and note that it is a vector (instead of matrix) given that catchability responses do not vary among variables \code{c}
#'        by default (but can be specified to do so when an appropriate `formula` is supplied)
#' @param Q2config_k Same as argument \code{Q1config_cp} but affecting affecting the 2nd linear predictor for catchability
#' @param spatial_list tagged list of locational information from , i.e., from \code{\link[FishStatsUtils]{make_spatial_info}}
#' @param PredTF_i OPTIONAL, whether each observation i is included in the likelihood, \code{PredTF_i[i]=0}, or in the predictive probability, \code{PredTF_i[i]=1}
#' @param Aniso whether to assume isotropy, \code{Aniso=0}, or geometric anisotropy, \code{Aniso=1}
#' @param Z_gm matrix specifying coordinates to use when calculating center-of-gravity and range-edge statistics.
#'        Defaults to eastings and northings for each knots or extrapolation-grid cell.
#' @param Expansion_cz matrix specifying how densities are expanded when calculating annual indices, with a row for each category \code{c} and two columns.
#'        The first column specifies whether to calculate annual index for category \code{c} as the weighted-sum across density estimates,
#'        where density is weighted by area ("area-weighted expansion", \code{Expansion[c,1]=0}, the default),
#'        where density is weighted by the expanded value for another category ("abundance weighted expansion" \code{Expansion[c1,1]=1}),
#'        or the index is calculated as the weighted average of density weighted by the expanded value for another category
#'        ("abundance weighted-average expansion" \code{Expansion[c1,1]=2}).  The 2nd column is only used when \code{Expansion[c1,1]=1} or \code{Expansion[c1,1]=2},
#'        and specifies the category to use for abundance-weighted expansion, where \code{Expansion[c1,2]=c2} and \code{c2} must be lower than \code{c1}.
#' @param F_ct matrix of fishing mortality for each category c and year t (only feasible when using a Poisson-link delta model
#'        and specifying temporal structure on intercepts, and mainly interpretable when species interactions via VamConfig)
#' @param Options a tagged-vector that is empty by default, \code{Options=c()}, but where the following slots might be helpful to add,
#'        either by passing \code{Options} to \code{\link[FishStatsUtils]{make_settings}}, or editing after a call to that function:
#' \describe{
#'   \item{\code{Options["SD_site_logdensity"]=TRUE}}{Turns on standard error calculation for local log-density (which is very slow to calculate!)}
#'   \item{\code{Options["Calculate_Range"]=TRUE}}{Turns on internal calculation and SE for center-of-gravity}
#'   \item{\code{Options["Calculate_effective_area"]=TRUE}}{Turns on internal calculation and SE for effective area occupied measuring range expansion/contraction}
#'   \item{\code{Options["Calculate_Cov_SE"]=TRUE}}{Turns on internal calculation and SE for covariance among categories (i.e. in factor model)}
#'   \item{\code{Options["Calculate_proportion"]=TRUE}}{Turns on internal calculation and SE for proportion of response within each category (e.g., for calculating proportion-at-age or species turnover)}
#'   \item{\code{Options["Calculate_Synchrony"]=TRUE}}{Turns on internal calculation and SE for Loreau metric of synchrony (a.k.a. portfolio effects)}
#'   \item{\code{Options["report_additional_variables"]=TRUE}}{Export additional variables to \code{Report} object, to use for diagnostics or additional exploration}
#' }
#' @param yearbounds_zz matrix with two columns, giving first and last years for defining one or more periods (rows) used to
#'        calculate changes in synchrony over time (only used if \code{Options['Calculate_Synchrony']=1})
#' @param CheckForErrors whether to check for errors in input (NOTE: when \code{CheckForErrors=TRUE}, the function will throw an error if
#'        it detects a problem with inputs.  However, failing to throw an error is no guarantee that the inputs are all correct)
#' @param ... interface to pass deprecated inputs, included for backwards compatibility with previous versions which, e.g., specified elements of \code{spatial_list}
#'        individually instead of as a single object

#' @return Object of class \code{make_data}, containing inputs to function \code{\link{make_model}}

#' @export
make_data <-
function( b_i,
          a_i,
          t_i,
          c_iz = rep(0,length(b_i)),
          e_i = c_iz[,1],
          v_i = rep(0,length(b_i)),
          FieldConfig = c("Omega1" = "IID","Epsilon1" = "IID","Omega2" = "IID","Epsilon2" = "IID"),
          ObsModel_ez = c("PosDist" = 1,"Link" = 0),
          OverdispersionConfig = c("eta1" = 0,"eta2" = 0),
          RhoConfig = c("Beta1" = 0,"Beta2" = 0,"Epsilon1" = 0,"Epsilon2" = 0),
          VamConfig = c("Method" = 0,"Rank" = 0,"Timing" = 0),
          Aniso = TRUE,
          PredTF_i = rep(0,length(b_i)),
          covariate_data = NULL,
          X1_formula = ~0,
          X2_formula = ~0,
          X1config_cp = NULL,
          X2config_cp = NULL,
          catchability_data = NULL,
          Q1_formula = ~0,
          Q2_formula = ~0,
          Q1config_k = NULL,
          Q2config_k = NULL,
          spatial_list,
          Network_sz = NULL,
          F_ct = NULL,
          F_init = 1,
          CheckForErrors = TRUE,
          yearbounds_zz = NULL,
          Options = c(),
          Expansion_cz = NULL,
          Z_gm = NULL,
          Version = FishStatsUtils::get_latest_version(package = "VAST"),
          overlap_zz = matrix(ncol = 7,nrow = 0),
          ... ){

  # Deprecated inputs for backwards compatibility e.g., in transition from CPP Version < 8.0.0 to >= 8.0.0
  alternate_inputs = list( ... )
  if( missing(spatial_list) ){
    warning("Consider changing use of `make_data` to include `spatial_list` as input")
    a_xl = a_gl = alternate_inputs[["a_xl"]]
    MeshList = alternate_inputs[["MeshList"]]
    GridList = alternate_inputs[["GridList"]]
    Method = alternate_inputs[["Method"]]
    s_i = alternate_inputs[["s_i"]]
  }else{
    MeshList = spatial_list[["MeshList"]]
    GridList = spatial_list[["GridList"]]
    Method = spatial_list[["Method"]]
    a_xl = a_gl = spatial_list[["a_gl"]]
    s_i = spatial_list[["knot_i"]] - 1
  }
  # Maintain interface for t_iz
  if(missing(t_i) & !is.null(alternate_inputs[["t_iz"]])){
    message( "Detecting deprecated input `t_iz` and coercing this to expected input `t_i`" )
    t_i = alternate_inputs$t_iz[,1]
  }
  if( "X_xj" %in% names(alternate_inputs) ) stop("`X_xj` is fully deprecated; please check inputs")

  # Specify default values for `Options`
  Options2use = c('SD_site_density'=FALSE, 'SD_site_logdensity'=FALSE, 'Calculate_Range'=FALSE, 'SD_observation_density'=FALSE, 'Calculate_effective_area'=FALSE,
    'Calculate_Cov_SE'=FALSE, 'Calculate_Synchrony'=FALSE, 'Calculate_Coherence'=FALSE, 'Calculate_proportion'=FALSE, 'normalize_GMRF_in_CPP'=TRUE,
    'Calculate_Fratio'=FALSE, 'Estimate_B0'=FALSE, 'Project_factors'=FALSE, 'treat_nonencounter_as_zero'=FALSE, 'simulate_random_effects'=TRUE,
    'observation_error_as_CV'=TRUE, 'report_additional_variables'=FALSE, 'zerosum_penalty'=0 )

  # Replace defaults for `Options` with provided values (if any)
  for( i in seq_along(Options) ){
    if(tolower(names(Options)[i]) %in% tolower(names(Options2use))){
      Options2use[[match(tolower(names(Options)[i]),tolower(names(Options2use)))]] = Options[[i]]
    }
  }

  #### Deals with backwards compatibility for FieldConfig
  # Converts from 4-vector to 3-by-2 matrix
  if( is.vector(FieldConfig) && length(FieldConfig)==4 ){
    FieldConfig = rbind( matrix(FieldConfig,ncol=2,dimnames=list(c("Omega","Epsilon"),c("Component_1","Component_2"))), "Beta"=c("IID","IID") )
  }
  # Converts from 3-by-2 matrix to 4-by-2 matrix
  if( is.matrix(FieldConfig) & all(dim(FieldConfig)==c(3,2)) ){
    FieldConfig = rbind( FieldConfig, "Epsilon_time"=c("Identity","Identity") )
  }
  # Checks for errors
  if( !is.matrix(FieldConfig) || !all(dim(FieldConfig)==c(4,2)) ){
    stop("`FieldConfig` has the wrong dimensions in `make_data`")
  }
  # Renames
  dimnames(FieldConfig) = list( c("Omega","Epsilon","Beta","Epsilon_time"), c("Component_1","Component_2") )

  # Rescale tprime_iz to start at 0
  tprime_i = t_i - min(t_i,na.rm=TRUE)

  # Increment first tprime_iz if t=0 corresponds to B0
  if( Options2use[12]==1 ){
    tprime_i = tprime_i + 1
    F_ct = cbind( 0, F_ct )
  }

  # Coerce c_iz to be a matrix
  if( !is.matrix(c_iz) ) c_iz = matrix(c_iz,ncol=1)

  # Determine dimensions
  n_t = max(tprime_i,na.rm=TRUE) + 1
  n_c = max(c_iz,na.rm=TRUE) + 1
  n_e = max(e_i) + 1
  n_v = length(unique(v_i))   # If n_v=1, then turn off overdispersion later
  n_i = length(b_i)
  n_x = nrow(a_gl)
  n_l = ncol(a_gl)
  n_g = ifelse( is.null(spatial_list), 1, spatial_list$n_g )

  # Coerce ObsModel_ez to be a matrix
  if( !is.matrix(ObsModel_ez) ) ObsModel_ez = matrix( ObsModel_ez, ncol=2, nrow=n_e, byrow=TRUE )

  # if Options2use['treat_nonencounter_as_zero']==TRUE, then replace b_i for all year-categories with all zeros with NAs
  if( Options2use['treat_nonencounter_as_zero']==TRUE ){
    # Determine year-category pairs with no data
    Index = list( factor(c_iz[,1],levels=0:max(c_iz[,1])), factor(tprime_i,levels=0:max(tprime_i)) )
    Num_ct = tapply( b_i, INDEX=Index, FUN=function(vec){sum(vec>0,na.rm=TRUE)} )
    Num_ct = ifelse( is.na(Num_ct), 0, Num_ct )
    b_i = ifelse( Num_ct[cbind(as.numeric(Index[[1]]),as.numeric(Index[[2]]))]==0, NA, b_i )
  }

  ####################
  # Density and catchability covariates
  #   Note:  Section is complicated to maintain backwards compatibility across many different specifications
  # User can specify density covariates in multiple ways (where the first possible will be preferentially used, and order ensures backwards compatibility):
  #   1: X1_ip OR X1_itp, AND X1_gctp OR X1_gtp, AND X2_ip OR X2_itp, AND X2_gctp OR X2_gtp
  #   2: X_ip OR X_itp, AND X_gctp OR X_gtp
  #   3: formula AND covariate_data
  #   4: covariate_data, and then using X1_formula and X2_formula defaults
  #   5: covariate_data AND X_gctp AND X_itp are all NULL
  # Code creates and uses objects internally
  #   X1_itp , X2_itp (where X_itp = X1_itp for backwards compatibility)
  #   X1_gctp , X2_gctp (where X_gctp = X1_gctp for backwards compatibility)
  #   X1config_cp , X2config_cp
  #   X_xtp (for backwards compatibility, which hasn't yet been entirely deprecated)
  ####################

  # Inputs that are deprecated or avoid user-interface
  X_xtp = alternate_inputs[["X_xtp"]]
  X_gctp = alternate_inputs[["X_gctp"]]
  X1_gctp = alternate_inputs[["X1_gctp"]]
  X2_gctp = alternate_inputs[["X2_gctp"]]
  X_itp = alternate_inputs[["X_itp"]]
  X1_itp = alternate_inputs[["X1_itp"]]
  X2_itp = alternate_inputs[["X2_itp"]]

  # Backwards compatibility for covariate interface from VAST <= 3.5.0 to >= 3.6.0
  if( missing(X1_formula) & missing(X2_formula) & "formula" %in% names(alternate_inputs) ){
    X1_formula = X2_formula = alternate_inputs[["formula"]]
  }
  if( is.null(X1config_cp) & is.null(X2config_cp) & "Xconfig_zcp" %in% names(alternate_inputs) ){
    X1config_cp = array( alternate_inputs[["Xconfig_zcp"]][1,,], dim=dim(alternate_inputs[["Xconfig_zcp"]])[2:3] )
    X2config_cp = array( alternate_inputs[["Xconfig_zcp"]][2,,], dim=dim(alternate_inputs[["Xconfig_zcp"]])[2:3] )
  }

  # Add options for supplying X_ip but not X_itp
  if( is.null(X_itp) & "X_ip" %in% names(alternate_inputs) ){
    X_itp = aperm( outer(alternate_inputs[["X_ip"]],rep(1,n_t)), c(1,3,2) )
  }
  # Add options for supplying X_gtp but not X_gctp
  if( is.null(X_gctp) & "X_gtp" %in% names(alternate_inputs) ){
    X_gctp = aperm( outer(alternate_inputs[["X_gtp"]],rep(1,n_c)), c(1,4,2,3) )
  }

  # Check for backwards-compatibility issues
  if( FishStatsUtils::convert_version_name(Version) >= FishStatsUtils::convert_version_name("VAST_v8_0_0") ){
    if( !is.null(X_xtp) ){
      stop("`X_xtp` is not used in version >= 8.0.0. If you'd like to specify covariates using input `X_xtp` please use `Version='VAST_v7_0_0'`")
    }
  }
  if( FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v7_0_0") ){
    if( !is.null(X_gctp) | !is.null(X_itp) ){
      stop("`X_gctp` and `X_itp` are not used in version <= 7.0.0. If you'd like to specify covariates using input `X_gctp` and `X_itp` please use `Version='VAST_v8_0_0'` or higher")
    }
  }
  if( FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v11_0_0") ){
    if( X1_formula != ~0 | X1_formula != ~0 ){
      stop("`X1_formula` and `X2_formula`; to use these please use `Version='VAST_v12_0_0'` or higher")
    }
  }

  # Specify objects directly using names >= 10.0.0
  Covariates_created = FALSE
  if( Covariates_created==FALSE ){
    if( !is.null(X1_gctp) & !is.null(X1_itp) & !is.null(X2_gctp) & !is.null(X2_itp) ){
      Covariates_created = TRUE
      if( !is.array(X1_gctp) || !(all(dim(X1_gctp)[1:3]==c(n_g,n_c,n_t))) ){
        stop("`X1_gctp` has wrong dimensions")
      }
      if( !is.array(X1_itp) || !(all(dim(X1_itp)[1:2]==c(n_i,n_t))) ){
        stop("`X1_itp` has wrong dimensions")
      }
      if( !is.array(X2_gctp) || !(all(dim(X2_gctp)[1:3]==c(n_g,n_c,n_t))) ){
        stop("`X2_gctp` has wrong dimensions")
      }
      if( !is.array(X2_itp) || !(all(dim(X2_itp)[1:2]==c(n_i,n_t))) ){
        stop("`X2_itp` has wrong dimensions")
      }
    }
  }
  # Specify objects using names <= 9.4.0
  if( Covariates_created==FALSE ){
    if( !is.null(X_gctp) & !is.null(X_itp) ){
      Covariates_created = TRUE
      X1_gctp = X2_gctp = X_gctp
      X1_itp = X2_itp = X_itp
      if( !is.array(X_gctp) || !(all(dim(X_gctp)[1:3]==c(n_g,n_c,n_t))) ){
        stop("`X_gctp` has wrong dimensions")
      }
      if( !is.array(X_itp) || !(all(dim(X_itp)[1:2]==c(n_i,n_t))) ){
        stop("`X_itp` has wrong dimensions")
      }
    }
  }
  # Specify formula and data-frame using names <= 9.4.0
  if( Covariates_created==FALSE ){
    if( !is.null(covariate_data) ){
      if( "formula" %in% names(alternate_inputs) ){
        Covariates_created = TRUE
        warning("Using input `formula` to generate covariates. This interface is soft-deprecated but still available for backwards compatibility; please switch to using `X1_formula` and `X2_formula`")
        covariate_list = FishStatsUtils::make_covariates( formula=alternate_inputs[["formula"]], covariate_data=covariate_data, Year_i=t_i,
          spatial_list=spatial_list, extrapolation_list=extrapolation_list )
        X1_gtp = X2_gtp = covariate_list$X_gtp
        X1_itp = X2_itp = covariate_list$X_itp
        X1_gctp = X2_gctp = aperm( outer(X1_gtp,rep(1,n_c)), c(1,4,2,3) )
      }
    }
  }
  # Specify formula and data-frame using names >= 10.0.0
  if( Covariates_created==FALSE ){
    if( !is.null(covariate_data) ){
      Covariates_created = TRUE
      if( FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v9_4_0") ){
        stop("To use separate formula interface for linear predictors, please use version >= CPP 10.0.0")
      }

      # First linear predictor
      covariate_list = FishStatsUtils::make_covariates( formula=X1_formula, covariate_data=covariate_data, Year_i=t_i,
        spatial_list=spatial_list, extrapolation_list=extrapolation_list )
      X1_gtp = covariate_list$X_gtp
      X1_itp = covariate_list$X_itp
      X1_gctp = aperm( outer(X1_gtp,rep(1,n_c)), c(1,4,2,3) )

      # Second linear predictor
      covariate_list = FishStatsUtils::make_covariates( formula=X2_formula, covariate_data=covariate_data, Year_i=t_i,
        spatial_list=spatial_list, extrapolation_list=extrapolation_list )
      X2_gtp = covariate_list$X_gtp
      X2_itp = covariate_list$X_itp
      X2_gctp = aperm( outer(X2_gtp,rep(1,n_c)), c(1,4,2,3) )
    }
  }
  # Fill in if missing
  if( Covariates_created==FALSE ){
    if( is.null(covariate_data) ){
      Covariates_created = TRUE
      X1_gctp = X2_gctp = array(0, dim=c(n_g,n_c,n_t,0))
      X1_itp = X2_itp = array(0, dim=c(n_i,n_t,0))
    }
  }

  # X1config / X2config defaults
  create_Xconfig = function( Xconfig_cp, n_c, n_p ){
    if( is.null(Xconfig_cp) ){
      Xconfig_cp = array(1, dim=c(n_c,n_p))
    }else{
      if( !is.array(Xconfig_cp) || !(all(dim(Xconfig_cp)==c(n_c,n_p))) ){
        stop("`Xconfig_cp` has wrong dimensions")
      }
      if( !all(Xconfig_cp %in% c(-1,0,1,2,3)) ){
        stop("`Xconfig_cp` has some wrong element(s)")
      }
      if( any(Xconfig_cp %in% -1) ){
        warning( "Using `Xconfig_cp[] = -1` is unconventional and warrents caution" )
      }
    }
    return( Xconfig_cp )
  }
  X1config_cp = create_Xconfig( Xconfig_cp=X1config_cp, n_c=n_c, n_p=dim(X1_gctp)[4] )
  X2config_cp = create_Xconfig( Xconfig_cp=X2config_cp, n_c=n_c, n_p=dim(X2_gctp)[4] )

  ####################
  # Catchability covariates
  # User can specify catchability covariates in multiple ways (where the first possible will be preferentially used, and order ensures backwards compatibility):
  #   1: Q_ik
  #   2: catchability_data, and then using Q1_formula and Q2_formula defaults
  #   3: catchability_data AND Q_ik are NULL
  # Code creates and uses objects internally
  #   Q1_ik and Q2_ik (where Q_ik = Q1_ik for backwards compatibility)
  #   Q1config_k , Q2config_k
  ####################

  # Check for backwards-compatibility issues
  if( FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v11_0_0") ){
    if( Q1_formula != ~0 | Q1_formula != ~0 ){
      stop("`Q1_formula` and `Q1_formula`; to use these please use `Version='VAST_v12_0_0'` or higher")
    }
  }

  Catchability_created = FALSE
  if( !is.null(alternate_inputs[["Q_ik"]]) ){
    if( !is.array(alternate_inputs[["Q_ik"]]) || nrow(alternate_inputs[["Q_ik"]])!=n_i ){
      stop("`Q_ik` has wrong dimensions")
    }
    Q1_ik = Q2_ik = alternate_inputs[["Q_ik"]]
    Catchability_created = TRUE
  }
  if( Catchability_created==FALSE ){
    if( !is.null(catchability_data) ){
      if( nrow(catchability_data)!=n_i ) stop("`catchability_data` has the wrong number of rows; please supply one row for each observation `i`")
      Catchability_created = TRUE
      # First predictor
      Model_matrix1 = stats::model.matrix( stats::update.formula(Q1_formula, ~.+1), data=catchability_data )
      Columns_to_keep = which( attr(Model_matrix1,"assign") != 0 )
      coefficient_names_Q1 = attr(Model_matrix1,"dimnames")[[2]][Columns_to_keep]
      Q1_ik = Model_matrix1[,Columns_to_keep,drop=FALSE]
      # First predictor
      Model_matrix2 = stats::model.matrix( stats::update.formula(Q2_formula, ~.+1), data=catchability_data )
      Columns_to_keep = which( attr(Model_matrix2,"assign") != 0 )
      coefficient_names_Q2 = attr(Model_matrix2,"dimnames")[[2]][Columns_to_keep]
      Q2_ik = Model_matrix2[,Columns_to_keep,drop=FALSE]
    }
  }
  if( Catchability_created==FALSE ){
    #Q1_ik = Q2_ik = matrix(0, nrow=n_i, ncol=1)
    Q1_ik = Q2_ik = matrix(0, nrow=n_i, ncol=0)
  }

  # X1config / X2config defaults
  create_Qconfig = function( Qconfig_k, n_k ){
    if( is.null(Qconfig_k) ){
      Qconfig_k = array(1, dim=c(n_k))
    }else{
      if( !is.vector(Qconfig_k) || length(Qconfig_k)!=n_k ){
        stop("`Qconfig_k` has wrong length")
      }
      if( !all(Qconfig_k %in% c(0,1,2,3)) ){
        stop("`Qconfig_k` has some wrong element(s)")
      }
    }
    return( Qconfig_k )
  }
  Q1config_k = create_Qconfig( Qconfig_k=Q1config_k, n_k=ncol(Q1_ik) )
  Q2config_k = create_Qconfig( Qconfig_k=Q2config_k, n_k=ncol(Q2_ik) )

  ##################
  # more defaults
  ##################

  if( is.null(yearbounds_zz)){
    yearbounds_zz = matrix(c(0,n_t-1),nrow=1)
  }else{
    if( !is.array(yearbounds_zz) || !(all(dim(yearbounds_zz)[2]==2)) ){
      stop("`yearbounds_zz` has wrong dimensions")
    }
  }
  t_yz = matrix(0:max(tprime_i,na.rm=TRUE), ncol=1)
  n_y = nrow(t_yz)

  # Other defaults
  if( is.null(F_ct) ){
    F_ct = matrix(0, nrow=n_c, ncol=n_t)
  }else{
    if( !is.array(F_ct) || !(all(dim(F_ct)==c(n_c,n_t))) ){
      stop("`F_ct` has wrong dimensions")
    }
  }
  if( is.null(Expansion_cz) ){
    Expansion_cz = matrix( 0, nrow=n_c, ncol=2 )
  }else{
    if( !is.array(Expansion_cz) || !(all(dim(Expansion_cz)==c(n_c,2))) ){
      stop("`Expansion_cz` has wrong dimensions")
    }
  }

  # Translate FieldConfig from input formatting to CPP formatting
  FieldConfig_input = array(NA, dim=dim(FieldConfig), dimnames=dimnames(FieldConfig) )
  g = function(mat) suppressWarnings( array(as.numeric(mat),dim=dim(mat)) )
  FieldConfig_input[] = ifelse( tolower(FieldConfig)=="ar1", 0, FieldConfig_input)
  FieldConfig_input[] = ifelse( tolower(FieldConfig)=="iid", -2, FieldConfig_input)
  FieldConfig_input[] = ifelse( tolower(FieldConfig)=="identity", -3, FieldConfig_input)
  FieldConfig_input[] = ifelse( tolower(FieldConfig)=="full", n_c, FieldConfig_input)
  FieldConfig_input[1:3,] = ifelse( !is.na(g(FieldConfig[1:3,])) & g(FieldConfig[1:3,])>0 & g(FieldConfig[1:3,])<=n_c, g(FieldConfig[1:3,]), FieldConfig_input[1:3,])
  FieldConfig_input[4,] = ifelse( !is.na(g(FieldConfig[4,,drop=FALSE])) & g(FieldConfig[4,,drop=FALSE])>0 & g(FieldConfig[4,,drop=FALSE])<=n_t, g(FieldConfig[4,,drop=FALSE]), FieldConfig_input[4,,drop=FALSE])
  FieldConfig_input[] = ifelse( !is.na(g(FieldConfig)) & g(FieldConfig)==0, -1, FieldConfig_input)
  if( any(is.na(FieldConfig_input)) ) stop( "'FieldConfig' must be: 0 (turn off overdispersion); 'IID' (independent for each factor); 'AR1' (use AR1 structure); or 0<n_f<=n_c (factor structure)" )
  message( "FieldConfig_input is:" )
  print(FieldConfig_input)

  # Translate OverdispersionConfig from input formatting to CPP formatting
  OverdispersionConfig_input = rep(NA, length(OverdispersionConfig))
  names(OverdispersionConfig_input) = names(OverdispersionConfig)
  g = function(vec) suppressWarnings(as.numeric(vec))
  OverdispersionConfig_input[] = ifelse( OverdispersionConfig=="AR1", 0, OverdispersionConfig_input)
  OverdispersionConfig_input[] = ifelse( !is.na(g(OverdispersionConfig)) & g(OverdispersionConfig)>0 & g(OverdispersionConfig)<=n_c, g(OverdispersionConfig), OverdispersionConfig_input)
  OverdispersionConfig_input[] = ifelse( !is.na(g(OverdispersionConfig)) & g(OverdispersionConfig)==0, -1, OverdispersionConfig_input)
  if( all(OverdispersionConfig_input<0) ){
    v_i = rep(0,length(b_i))
    n_v = 1
  }
  if( any(is.na(OverdispersionConfig_input)) ) stop( "'OverdispersionConfig' must be: 0 (turn off overdispersion); 'AR1' (use AR1 structure); or 0<n_f<=n_c (factor structure)" )
  message( "OverdispersionConfig_input is:" )
  print(OverdispersionConfig_input)

  # by default, add nothing as Z_xl
  if( Options2use['Calculate_Range']==FALSE ){
    Z_gm = Z_xm = matrix(0, nrow=nrow(a_gl), ncol=ncol(a_gl) ) # Size so that it works for Version 3g-3j
  }else{
    if( is.null(spatial_list) ){
      # Maintain backwards compatibility with Version < 8.0.0
      Z_xm = Z_gm = MeshList$loc_x
    }else{
      # Improved interface used for Version >= 8.0.0
      Z_xm = spatial_list$loc_x
      if(is.null(Z_gm)){
        Z_gm = spatial_list$loc_g
      }else{
        if(nrow(Z_gm)!=nrow(spatial_list$loc_g)) stop("Check dimensions for input `Z_gm`")
      }
    }
    message( "Calculating range shift for stratum #1:",colnames(a_gl[1]))
  }

  ###################
  # Check for bad data entry
  ###################

  if( CheckForErrors==TRUE ){
    if( ncol(ObsModel_ez)!=2 | nrow(ObsModel_ez)!=n_e ) stop("Wrong dimensions for ObsModel_ez")
    if( !is.matrix(a_gl) ) stop("a_gl should be a matrix")
    if( !is.array(X1_gctp) | !is.array(X1_itp) ) stop( "X1_gctp and X1_itp should be arrays")
    if( (max(s_i)-1)>n_x | min(s_i)<0 ) stop("s_i exceeds bounds in MeshList")
    if( any(a_i<=0) ) stop("a_i must be greater than zero for all observations, and at least one value of a_i is not")
    # Logic-check that all categories have data
    if( !all(0:(n_c-1) %in% c_iz) ) stop("One or more categories has no associated observation: please check input `c_iz` for errors")
    # Warnings about all positive or zero
    Prop_nonzero = tapply( b_i, INDEX=list(tprime_i,c_iz[,1]), FUN=function(vec){mean(vec>0)} )
    if( Options2use[12]==1 ){
      Prop_nonzero = Prop_nonzero[-1,]
    }
    if( any(ObsModel_ez[,2] %in% c(0,1)) ){
      if( RhoConfig[1]==0 ){
        if( any(!is.na(Prop_nonzero) & (Prop_nonzero==1)) ){
          print( Prop_nonzero )
          stop("Some years and/or categories have 100% encounters, and this requires either temporal structure of a different link-function")
        }
      }
    }
    if( any(ObsModel_ez[,2] %in% c(0,1,3,4)) ){
      # Require Options2use['treat_nonencounter_as_zero']=TRUE to use any standard link function without temporal smoother given that there's 0% encounters
      if( RhoConfig[1]==0 ){
        if( any(!is.na(Prop_nonzero) & (Prop_nonzero==0)) & Options2use['treat_nonencounter_as_zero']==FALSE ){
          print( Prop_nonzero )
          stop("Some years and/or categories have 0% encounters, and this requires either temporal structure of a different link-function")
        }
      }
    }
    if( length(OverdispersionConfig)!=2 ) stop("length(OverdispersionConfig)!=2")
    if( ncol(yearbounds_zz)!=2 ) stop("yearbounds_zz must have two columns")
    if( Options2use['Calculate_Coherence']==1 & any(ObsModel_ez[,2]==0) ) stop("Calculating coherence only makes sense when 'ObsModel_ez[,2]=1'")
    if( any(yearbounds_zz<0) | any(yearbounds_zz>=max(n_t)) ) stop("yearbounds_zz exceeds bounds for modeled years")
    if( n_c!=length(unique(na.omit(as.vector(c_iz)))) ) stop("n_c doesn't equal the number of levels in c_i")
    if( any(ObsModel_ez[,1]==9) & !all(b_i%in%0:3) ) stop("If using 'ObsModel_ez[e,1]=9', all 'b_i' must be in {0,1,2,3}")
    if( length(unique(ObsModel_ez[,2]))>1 ) stop("All `ObsModel_ez[,2]` must have the same value")
    if( any(OverdispersionConfig>0) & length(unique(v_i))==1 ) stop("It doesn't make sense to use use `OverdispersionConfig` when using only one level of `v_i`")
    if( any(ObsModel_ez[,1] %in% c(12,13,14)) ){
      if( any(ObsModel_ez[,2] != 1) ) stop("If using `ObsModel_ez[e,1]` in {12,13,14} then must use `ObsModel_ez[e,2]=1`")
      if( !any(ObsModel_ez[,1] %in% c(0,1,2,3,4)) ) stop("Using `ObsModel_ez[e,1]` in {12,13,14} is only intended when combining data with biomass-sampling data")
    }
    if( all(b_i>0) & all(ObsModel_ez[,1]==0) & !all(FieldConfig_input[1:2,1]==-1) ) stop("All data are positive and using a conventional delta-model, so please turn off `Omega1` and `Epsilon1` terms")
    if( !(all(ObsModel_ez[,1] %in% c(0,1,2,4,5,7,10,11,12,13,14))) ) stop("Please check `ObsModel_ez[,1]` input")
    if( !(all(ObsModel_ez[,2] %in% c(0,1,2,3,4))) ) stop("Please check `ObsModel_ez[,2]` input")
    if( !all(RhoConfig[1]%in%c(0,1,2,3,4)) | !all(RhoConfig[2]%in%c(0,1,2,3,4,6)) | !all(RhoConfig[3]%in%c(0,1,2,4,5)) | !all(RhoConfig[4]%in%c(0,1,2,4,5,6)) ) stop("Check `RhoConfig` inputs")
    if( any(is.na(X_xtp)) ) stop("Some `X_xtp` is NA, and this is not allowed")
    if( any(is.na(X1_gctp)) ) stop("Some `X1_gctp` is NA, and this is not allowed")
    if( any(is.na(X1_itp)) ) stop("Some `X1_itp` is NA, and this is not allowed")
    if( n_c==1 && !all(FieldConfig_input[1:3,] %in% c(-3,-2,-1,1)) ) stop("If using a univariate model, `FieldConfig` must be 0, 1, or `IID` for all entries")
  }

  # Check for wrong dimensions
  if( CheckForErrors==TRUE ){
    if( any(c(length(b_i),length(a_i),nrow(c_iz),length(tprime_i),length(v_i),length(PredTF_i))!=n_i) ) stop("b_i, a_i, c_i, s_i, v_i, or tprime_i doesn't have length n_i")
    if( nrow(a_gl)!=n_x | ncol(a_gl)!=n_l ) stop("a_xl has wrong dimensions")
    if( FishStatsUtils::convert_version_name(Version) >= FishStatsUtils::convert_version_name("VAST_v8_0_0") ){
      if( any(dim(X1_gctp)[1:3] != c(n_g,n_c,n_t)) ) stop("X1_gctp has wrong dimensions")
      if( any(dim(X1_itp)[1:2] != c(n_i,n_t)) ) stop("X1_itp has wrong dimensions")
    }else{
      if( any(dim(X_xtp)[1:2] != c(n_x,n_t)) ) stop("X_xtp has wrong dimensions")
    }
    if( ncol(c_iz)>1 & any(ObsModel_ez[,2]!=1) ) stop("Using multiple columnns in `c_iz` only makes sense using a Poisson-link delta model via `ObsModel[2]=1`")
    if( nrow(F_ct)!=n_c | ncol(F_ct)!=n_t ) stop("F_ct has wrong dimensions")
    if( ncol(overlap_zz) != 7 ) stop("Input `overlap_zz` must contain 7 columns but doesn't")
    if( any(overlap_zz[,c(1,3)] >= n_c) ) stop("Check `overlap_zz[,c(1,3)]` entries")
    if( any(overlap_zz[,c(2,4)] >= n_t) ) stop("Check `overlap_zz[,c(2,4)]` entries")
  }

  ###################
  # Check for incompatibilities amongst versions
  ###################

  if( FishStatsUtils::convert_version_name(Version) >= FishStatsUtils::convert_version_name("VAST_v9_4_0") ){
    if(VamConfig[1]==3) stop("`VamConfig[1]=3` feature causes compile issues on macOS, and has been removed from the CPP; please contact package author if interested in using it.")
  }
  if( FishStatsUtils::convert_version_name(Version) >= FishStatsUtils::convert_version_name("VAST_v8_0_0") ){
    if( is.null(spatial_list) ) stop("Must provide `spatial_list` for Version >= 8.0.0")
  }
  if( any(ObsModel_ez[,2]==1) ){
    # Versions 1.6.0 through 2.2.0 used a previous interpretation of area-swept for Poisson-link model and are not consistent with Q-Q plotting
    if( FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v2_2_0") ){
      stop("Problem with Poisson-link model for VAST versions 1.6.0 through 2.2.0")
    }
    # Versions 1.0.0 through 1.5.0 don't have Poisson-link model
    if( FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v1_5_0") ){
      stop("Poisson-link model for VAST versions 1.0.0 through 1.5.0")
    }
    # Can't use multiple error distributions prior to version 3.0.0
    if( FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v2_8_0") ){
      if( length(unique(ObsModel_ez[,1]))>1 ) stop("Can't use multiple error distributions prior to version 3.0.0")
    }
  }
  if( any(FieldConfig_input == -2) ){
    # Versions 2.6.0 was the first to allow "IID" setting for FieldConfig_input elements
    if( FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v2_5_0") ){
      stop("Problem with 'IID' option for factors for VAST versions 1.0.0 through 2.6.0")
    }
  }
  if( ncol(c_iz)>1 ){
    if( FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v3_0_0") ){
      stop("Can't have observations assigned to more than one category prior to version 4.0.0")
    }
  }
  if( any(RhoConfig[1:2]==3) ){
    if( FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v4_1_0") ){
      stop("There was bug in specifying fixed intercepts among years for versions prior to V4.2.0")
    }
  }
  if( Options2use['SD_observation_density']==1 ){
    if( FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v4_1_0") ){
      stop("Calculating 'SD_observation_density' is not possible prior to V4.2.0")
    }
  }
  if( any(ObsModel_ez[,1]==3) ){
    if( FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v8_2_0") ){
      stop("Inverse-gaussian distribution only available for CPP version >= 8_3_0")
    }
  }

  ###################
  # Check for incompatible settings
  ###################

  # intercepts
  if( (FieldConfig_input[2,1]==(-1) & RhoConfig[3]!=0) | (FieldConfig_input[2,2]==(-1) & RhoConfig[4]!=0) ){
    stop("Spatio-temporal variation is turned off for a component with temporal structure, and this combination doesn't make sense")
  }

  # Prohibitively slow
  if( !is.null(spatial_list$fine_scale) && spatial_list$fine_scale==TRUE ){
    if( Options2use['SD_site_density']==TRUE | Options2use['SD_site_logdensity']==TRUE ){
      warning("'SD_site_density' and 'SD_site_logdensity' are very slow when using `fine_scale=TRUE`")
    }
  }

  # Recommend against dlnorm with mean+sd parameterization
  if( any(ObsModel_ez[,1]==1) ){
    #warning("the package author recommends using the lognormal mean-CV parameterization, `ObsModel_ez[,1]=4` instead of the mean-SD parameterization, `ObsModel_ez[,1]=1`")
  }
  if( any(ObsModel_ez[,1]==3) ){
    if( Version=="VAST_v8_3_0" ){
      warning("Version `VAST_v8_3_0` used a mean-sd parameterization for the Inverse-Gaussian, while later versions use a mean-CV parameterization to match the Gamma distribution")
    }
  }

  # Interactions
  if( VamConfig[1]!=0 ){
    if( any(ObsModel_ez[,2]!=1) ){
      stop("Must use Poisson-link delta model when estimating interactions")
    }
    if( any(RhoConfig[1:2]!=3) ){
      #stop("Must use constant intercepts when estimating interactions")
    }
    if( !(FieldConfig_input[2,1] %in% c(-2,-1,0,n_c)) & VamConfig[3]==1 ){
      stop("Spatio-temporal variation must either have full rank covariance or be turned off or IID for the 1st linear predictor for when using interactions and when VamConfig[`Timing`]==1")
    }
    if( VamConfig[2]>n_c | VamConfig[2]<0 ){
      stop("Rank for interactions must be an integer between 0 and `n_c`, inclusive")
    }
    if( VamConfig[2]==FieldConfig_input[2,1] & RhoConfig[3]!=0 ){
      stop("Can't simultaneously identify full-rank interactions and temporal correlation on spatio-temporal component for 1st linear predictor")
    }
  }

  # Mirroring
  if( RhoConfig[4]==6 ){
    if( FieldConfig_input[2,1]!=FieldConfig_input[2,2] ){
      stop("To fix 'Epsilon_rho2_f` equal to 'Epsilon_rho2_f`, you must specify the same rank using `FieldConfig_input[2,1]` equal to `FieldConfig_input[2,2]`")
    }
  }

  # Bratio reporting
  if( Options2use[12]!=0 ){
    if( FieldConfig_input[2,1]!=n_c ){
      stop("Must have full-rank spatio-temporal component matrix to estimate B0 using `Options2use[12]=1`")
    }
    if( !(Options2use[12] %in% c(0,1)) ){
      stop("`Options2use[12]` must be either 0 or 1")
    }
    if( any(ObsModel_ez[,2]!=1) ){
      stop("Must use Poisson-link delta model when estimating interactions")
    }
  }

  # Fratio reporting
  if( Options2use[11]!=0 ){
    if( FieldConfig_input[2,2]!=0 & !(RhoConfig[4] %in% c(6)) ){
      stop("To estimate Fratio, either Epsilon2 must be turned off (i.e., `FieldConfig_input[2,2]=0`) or B2 must equal B1_cc (i.e., `RhoConfig[4]=6`)")
    }
  }

  # Fishing mortality
  if( any(F_ct!=0) ){
    if( any(ObsModel_ez[,2]!=1) ){
      stop("Must use Poisson-link delta model when estimating the impact of fishing mortality")
    }
    if( any(RhoConfig[1:2]!=3) ){
      stop("Must use constant intercepts when estimating the impact of fishing mortality")
    }
    if( !(F_init %in% c(1,2)) ){
      stop("`F_init` must be either 1 or 2")
    }
    if( F_init==2 ){
      if( FieldConfig_input[2,2]!=0 & !(RhoConfig[4] %in% c(6)) ){
        stop("To estimate stationary distribution for initial F, either Epsilon2 must be turned off (i.e., `FieldConfig_input[2,2]=0`) or B2 must equal B1_cc (i.e., `RhoConfig[4]=6`)")
      }
    }
  }

  # Stream network stuff
  if( Method=="Stream_network" ){
    if(is.null(Network_sz)){
      stop("Must specify 'Network_sz' when using Method=='Stream_network'")
    }
    if( ncol(Network_sz)!=3 | !all(c("parent_s","child_s","dist_s") %in% colnames(Network_sz)) ){
      stop("'Network_sz' must have three columns, 'parent_s', 'child_s', and 'dist_s'")
    }
  }else{
    Network_sz = matrix( c(1,1,1), nrow=1, dimnames=list(NULL,c("parent_s","child_s","dist_s")) )
  }

  # Check incompatibility of constant intercept + 0/100% encounter option
  if( CheckForErrors==TRUE ){
    if( any(RhoConfig[1:2]!=0) & any(ObsModel_ez[,2]==3) ){
      stop( "RhoConfig[1:2] must be 0 when using ObsModel[2]=3:  Other options are not coded to work together" )
    }
    if( any(RhoConfig[1:2]!=0) & any(ObsModel_ez[,2]==4) ){
      stop( "RhoConfig[1:2] must be 0 when using ObsModel[2]=4:  Other options are not coded to work together" )
    }
  }

  # Factor model for intercepts + 0% or 100% encounter rate options doesn't make sense
  if( any(FieldConfig_input[3,1:2] != -2) & any(ObsModel_ez[,2] %in% c(3,4)) ){
    stop( "Factor model for intercepts is incompatible  with ObsModel_ez[,2] being 3 or 4")
  }

  # Rank-reduced factor model for intercepts + fixed intercepts doesn't make sense
  if( RhoConfig[1] == 0 ){
    if( FieldConfig_input[3,1] != -2 ) stop("Using a factor model doesn't make sense using fixed-effect intercepts.  If you want to use a factor model without temporal structure, please change `RhoConfig[1]=1` for covariance that is independent in each year, or use some other temporal structure on intercepts.")
  }
  if( RhoConfig[2] == 0 ){
    if( FieldConfig_input[3,2] != -2 ) stop("Using a factor model doesn't make sense using fixed-effect intercepts.  If you want to use a factor model without temporal structure, please change `RhoConfig[2]=1` for covariance that is independent in each year, or use some other temporal structure on intercepts.")
  }

  # Check for CV vs. variance parameterization for positive catch rates
  if( Options2use["observation_error_as_CV"]==FALSE ){
    if( !all(ObsModel_ez[,1] %in% c(2,3,4)) ){
      stop("Options['observation_error_as_CV']==FALSE only works with some observation distributions")
    }else{
      warning("Options['observation_error_as_CV']==FALSE is experimental; please use at your own risk")
    }
  }


  ###################
  # switch defaults if necessary
  ###################

  if( Method=="Grid" ){
    Aniso = 0
    message("Using isotropic 2D AR1 hyperdistribution, so switching to Aniso=0")
  }
  if( Method=="Spherical_mesh" ){
    Aniso = 0
    message("Using spherical projection for SPDE approximation, so switching to Aniso=0")
  }
  if( Method=="Stream_network" ){
    Aniso = 0
    message("Using stream network correlations, so switching to Aniso=0")
  }
  if( Method=="Barrier" ){
    Aniso = 0
    message("Using SPDE approach using geographical barriers, so switching to Aniso=0")
  }
  if( VamConfig[2]==0 & VamConfig[1]!=0 ){
    VamConfig[1] = 0
    message("Using interactions with zero rank (`VamConfig[2]==0`), so turning off interactions (`VamConfig[1]=0`)")
  }

  # Warning messages
  if( n_c>1 & any(FieldConfig_input==1)){
    warning( "Using 1 factor for more than one category:  Please note that this is non-standard, and it is more common to use multiple factors (often as many as the number of categories)" )
  }
  #SD_p = apply( X_xtp, MARGIN=3, FUN=sd )
  #if( any(SD_p>3) ){
  #  warning( "I highly recommend that you standardize each density covariate `X_xtp` to have a low standard deviation, to avoid numerical under/over-flow" )
  #}

  # Tweedie bug
  if( any(ObsModel_ez[,1]==10) ){
    if( FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v8_3_0") ){
      warning("CPP versions prior to 8.4.0 had a bug in dtweedie, where the power parameter xi was mistakenly constrained to range from 1.5 to 2.0 (rather than 1.0 to 2.0). This bug either had no effect, or resulted in the parameter hitting a bound and impaired model fit.")
    }
  }

  ###################
  # Output tagged list
  ###################

  # CMP_xmax should be >100 and CMP_breakpoint should be 1 for Tweedie model
  Options_vec = c("Aniso"=Aniso, "R2_interpretation"=0, "Rho_beta1_TF"=ifelse(RhoConfig[1]%in%c(1,2,4),1,0), "Rho_beta2_TF"=ifelse(RhoConfig[2]%in%c(1,2,4),1,0), "AreaAbundanceCurveTF"=0, "CMP_xmax"=200, "CMP_breakpoint"=1, "Method"=switch(Method,"Mesh"=0,"Grid"=1,"Spherical_mesh"=0,"Stream_network"=2,"Barrier"=3), "Include_F"=ifelse(all(F_ct==0),0,F_init) )
  Return = NULL
  if(Version%in%c("VAST_v1_1_0","VAST_v1_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=dim(X_xtp)[3], "n_k"=ncol(Q1_ik), "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=tprime_i, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
  }
  if(Version%in%c("VAST_v1_4_0","VAST_v1_3_0","VAST_v1_2_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=dim(X_xtp)[3], "n_k"=ncol(Q1_ik), "n_v"=n_v, "n_f_input"=OverdispersionConfig_input[1], "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=tprime_i, "v_i"=match(v_i,sort(unique(v_i)))-1, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
  }
  if(Version%in%c("VAST_v1_6_0","VAST_v1_5_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=dim(X_xtp)[3], "n_k"=ncol(Q1_ik), "n_v"=n_v, "n_f_input"=OverdispersionConfig_input[1], "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=tprime_i, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
  }
  if(Version%in%c("VAST_v1_7_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=dim(X_xtp)[3], "n_k"=ncol(Q1_ik), "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=tprime_i, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
  }
  if(Version%in%c("VAST_v1_8_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=dim(X_xtp)[3], "n_k"=ncol(Q1_ik), "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=tprime_i, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v1_9_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=dim(X_xtp)[3], "n_k"=ncol(Q1_ik), "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=tprime_i, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v2_8_0","VAST_v2_7_0","VAST_v2_6_0","VAST_v2_5_0","VAST_v2_4_0","VAST_v2_3_0","VAST_v2_2_0","VAST_v2_1_0","VAST_v2_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=dim(X_xtp)[3], "n_k"=ncol(Q1_ik), "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_iz"=matrix(tprime_i,ncol=1), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v3_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_j"=n_j, "n_p"=dim(X_xtp)[3], "n_k"=ncol(Q1_ik), "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "e_i"=e_i, "s_i"=s_i, "t_iz"=matrix(tprime_i,ncol=1), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v4_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_j"=n_j, "n_p"=dim(X_xtp)[3], "n_k"=ncol(Q1_ik), "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=matrix(tprime_i,ncol=1), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v4_4_0","VAST_v4_3_0","VAST_v4_2_0","VAST_v4_1_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=dim(X_xtp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "include_data"=TRUE, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=matrix(tprime_i,ncol=1), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v5_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=dim(X_xtp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "include_data"=TRUE, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=matrix(tprime_i,ncol=1), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v5_1_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=dim(X_xtp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "include_data"=TRUE, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=matrix(tprime_i,ncol=1), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v5_2_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=dim(X_xtp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "include_data"=TRUE, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=matrix(tprime_i,ncol=1), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v5_4_0","VAST_v5_3_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=dim(X_xtp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "include_data"=TRUE, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=matrix(tprime_i,ncol=1), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v5_5_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=dim(X_xtp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_list"=list("Options_vec"=Options_vec,"Options"=Options2use,"yearbounds_zz"=yearbounds_zz,"Expansion_cz"=Expansion_cz), "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "include_data"=TRUE, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=matrix(tprime_i,ncol=1), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=matrix(0,nrow=n_x,ncol=1), "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v6_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=dim(X_xtp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_list"=list("Options_vec"=Options_vec,"Options"=Options2use,"yearbounds_zz"=yearbounds_zz,"Expansion_cz"=Expansion_cz), "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "Xconfig_zcp"=abind::abind(X1config_cp,X2config_cp,along=3), "include_data"=TRUE, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=matrix(tprime_i,ncol=1), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v7_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=dim(X_xtp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_list"=list("Options_vec"=Options_vec,"Options"=Options2use,"yearbounds_zz"=yearbounds_zz,"Expansion_cz"=Expansion_cz), "FieldConfig"=FieldConfig_input, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "Xconfig_zcp"=abind::abind(X1config_cp,X2config_cp,along=3), "include_data"=TRUE, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=matrix(tprime_i,ncol=1), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xtp"=X_xtp, "Q_ik"=Q1_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v9_1_0","VAST_v9_0_0","VAST_v8_6_0","VAST_v8_5_0","VAST_v8_4_0","VAST_v8_3_0","VAST_v8_2_0","VAST_v8_1_0","VAST_v8_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_g"=n_g, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=dim(X1_itp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_gm), "Options_list"=list("Options_vec"=Options_vec,"Options"=Options2use,"yearbounds_zz"=yearbounds_zz,"Expansion_cz"=Expansion_cz), "FieldConfig"=FieldConfig_input, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "Xconfig_zcp"=abind::abind(X1config_cp,X2config_cp,along=3), "include_data"=TRUE, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "t_iz"=matrix(tprime_i,ncol=1), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_gl"=a_gl, "X_itp"=X1_itp, "X_gtp"=array(X1_gctp[,1,,],dim(X1_gctp)[c(1,3,4)]), "Q_ik"=Q1_ik, "t_yz"=t_yz, "Z_gm"=Z_gm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2, "Ais_ij"=cbind(spatial_list$A_is@i,spatial_list$A_is@j), "Ais_x"=spatial_list$A_is@x, "Ags_ij"=cbind(spatial_list$A_gs@i,spatial_list$A_gs@j), "Ags_x"=spatial_list$A_gs@x )
  }
  if(Version%in%c("VAST_v9_2_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_g"=n_g, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=dim(X1_itp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_gm), "Options_list"=list("Options_vec"=Options_vec,"Options"=Options2use,"yearbounds_zz"=yearbounds_zz,"Expansion_cz"=Expansion_cz,"overlap_zz"=overlap_zz), "FieldConfig"=FieldConfig_input, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "Xconfig_zcp"=abind::abind(X1config_cp,X2config_cp,along=3), "include_data"=TRUE, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "t_iz"=matrix(tprime_i,ncol=1), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_gl"=a_gl, "X_itp"=X1_itp, "X_gtp"=array(X1_gctp[,1,,],dim(X1_gctp)[c(1,3,4)]), "Q_ik"=Q1_ik, "t_yz"=t_yz, "Z_gm"=Z_gm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2, "Ais_ij"=cbind(spatial_list$A_is@i,spatial_list$A_is@j), "Ais_x"=spatial_list$A_is@x, "Ags_ij"=cbind(spatial_list$A_gs@i,spatial_list$A_gs@j), "Ags_x"=spatial_list$A_gs@x )
  }
  if(Version%in%c("VAST_v9_3_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x,MeshList$anisotropic_spde$n.spde)[Options_vec['Method']+1], "n_g"=n_g, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=dim(X1_itp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_gm), "Options_list"=list("Options_vec"=Options_vec,"Options"=Options2use,"yearbounds_zz"=yearbounds_zz,"Expansion_cz"=Expansion_cz,"overlap_zz"=overlap_zz), "FieldConfig"=FieldConfig_input, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "Xconfig_zcp"=abind::abind(X1config_cp,X2config_cp,along=3), "include_data"=TRUE, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "t_iz"=matrix(tprime_i,ncol=1), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_gl"=a_gl, "X_itp"=X1_itp, "X_gtp"=array(X1_gctp[,1,,],dim(X1_gctp)[c(1,3,4)]), "Q_ik"=Q1_ik, "t_yz"=t_yz, "Z_gm"=Z_gm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "spdeMatricesBarrier"=list(), "Barrier_scaling"=NA, "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2, "Ais_ij"=cbind(spatial_list$A_is@i,spatial_list$A_is@j), "Ais_x"=spatial_list$A_is@x, "Ags_ij"=cbind(spatial_list$A_gs@i,spatial_list$A_gs@j), "Ags_x"=spatial_list$A_gs@x )
  }
  if(Version%in%c("VAST_v9_4_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x,MeshList$anisotropic_spde$n.spde)[Options_vec['Method']+1], "n_g"=n_g, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=dim(X1_itp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_gm), "Options_list"=list("Options_vec"=Options_vec,"Options"=Options2use,"yearbounds_zz"=yearbounds_zz,"Expansion_cz"=Expansion_cz,"overlap_zz"=overlap_zz), "FieldConfig"=FieldConfig_input, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "Xconfig_zcp"=abind::abind(X1config_cp,X2config_cp,along=3), "include_data"=TRUE, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "t_i"=tprime_i, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_gl"=a_gl, "X_ip"=array(X1_itp[,1,],dim=dim(X1_itp)[c(1,3)]), "X_gctp"=X1_gctp, "Q_ik"=Q1_ik, "Z_gm"=Z_gm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "spdeMatricesBarrier"=list(), "Barrier_scaling"=NA, "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2, "Ais_ij"=cbind(spatial_list$A_is@i,spatial_list$A_is@j), "Ais_x"=spatial_list$A_is@x, "Ags_ij"=cbind(spatial_list$A_gs@i,spatial_list$A_gs@j), "Ags_x"=spatial_list$A_gs@x )
  }
  if(Version%in%c("VAST_v10_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x,MeshList$anisotropic_spde$n.spde)[Options_vec['Method']+1], "n_g"=n_g, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p1"=dim(X1_itp)[3], "n_p2"=dim(X2_itp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_gm), "Options_list"=list("Options_vec"=Options_vec,"Options"=Options2use,"yearbounds_zz"=yearbounds_zz,"Expansion_cz"=Expansion_cz,"overlap_zz"=overlap_zz), "FieldConfig"=FieldConfig_input, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "X1config_cp"=X1config_cp, "X2config_cp"=X2config_cp, "include_data"=TRUE, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "t_i"=tprime_i, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_gl"=a_gl, "X1_ip"=array(X1_itp[,1,],dim=dim(X1_itp)[c(1,3)]), "X1_gctp"=X1_gctp, "X2_ip"=array(X2_itp[,1,],dim=dim(X2_itp)[c(1,3)]), "X2_gctp"=X2_gctp, "Q_ik"=Q1_ik, "Z_gm"=Z_gm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "spdeMatricesBarrier"=list(), "Barrier_scaling"=NA, "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2, "Ais_ij"=cbind(spatial_list$A_is@i,spatial_list$A_is@j), "Ais_x"=spatial_list$A_is@x, "Ags_ij"=cbind(spatial_list$A_gs@i,spatial_list$A_gs@j), "Ags_x"=spatial_list$A_gs@x )
  }
  if(Version%in%c("VAST_v11_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x,MeshList$anisotropic_spde$n.spde)[Options_vec['Method']+1], "n_g"=n_g, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p1"=dim(X1_itp)[3], "n_p2"=dim(X2_itp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_gm), "Options_list"=list("Options_vec"=Options_vec,"Options"=Options2use,"yearbounds_zz"=yearbounds_zz,"Expansion_cz"=Expansion_cz,"overlap_zz"=overlap_zz), "FieldConfig"=FieldConfig_input, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "X1config_cp"=X1config_cp, "X2config_cp"=X2config_cp, "include_data"=TRUE, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "t_i"=tprime_i, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_gl"=a_gl, "X1_ip"=array(X1_itp[,1,],dim=dim(X1_itp)[c(1,3)]), "X1_gctp"=X1_gctp, "X2_ip"=array(X2_itp[,1,],dim=dim(X2_itp)[c(1,3)]), "X2_gctp"=X2_gctp, "Q1_ik"=Q1_ik, "Q2_ik"=Q2_ik, "Z_gm"=Z_gm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "spdeMatricesBarrier"=list(), "Barrier_scaling"=NA, "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2, "Ais_ij"=cbind(spatial_list$A_is@i,spatial_list$A_is@j), "Ais_x"=spatial_list$A_is@x, "Ags_ij"=cbind(spatial_list$A_gs@i,spatial_list$A_gs@j), "Ags_x"=spatial_list$A_gs@x )
  }
  if(Version%in%c("VAST_v13_0_0","VAST_v12_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x,MeshList$anisotropic_spde$n.spde)[Options_vec['Method']+1], "n_g"=n_g, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p1"=dim(X1_itp)[3], "n_p2"=dim(X2_itp)[3], "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_gm), "Options_list"=list("Options_vec"=Options_vec,"Options"=Options2use,"yearbounds_zz"=yearbounds_zz,"Expansion_cz"=Expansion_cz,"overlap_zz"=overlap_zz,"zerosum_penalty"=matrix(Options2use['zerosum_penalty'],1,1)), "FieldConfig"=FieldConfig_input, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "X1config_cp"=X1config_cp, "X2config_cp"=X2config_cp, "Q1config_k"=Q1config_k, "Q2config_k"=Q2config_k, "include_data"=TRUE, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "t_i"=tprime_i, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_gl"=a_gl, "X1_ip"=array(X1_itp[,1,],dim=dim(X1_itp)[c(1,3)]), "X1_gctp"=X1_gctp, "X2_ip"=array(X2_itp[,1,],dim=dim(X2_itp)[c(1,3)]), "X2_gctp"=X2_gctp, "Q1_ik"=Q1_ik, "Q2_ik"=Q2_ik, "Z_gm"=Z_gm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "spdeMatricesBarrier"=list(), "Barrier_scaling"=NA, "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2, "Ais_ij"=cbind(spatial_list$A_is@i,spatial_list$A_is@j), "Ais_x"=spatial_list$A_is@x, "Ags_ij"=cbind(spatial_list$A_gs@i,spatial_list$A_gs@j), "Ags_x"=spatial_list$A_gs@x )
  }
  if( is.null(Return) ) stop("`Version` provided does not match the list of possible values")
  if( "spde" %in% names(Return) ) Return[['spde']] = MeshList$isotropic_spde$param.inla[c("M0","M1","M2")]
  if( "spde_aniso" %in% names(Return) ) Return[['spde_aniso']] = list("n_s"=MeshList$anisotropic_spde$n.spde, "n_tri"=nrow(MeshList$anisotropic_mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$anisotropic_spde$param.inla$M0, "G0_inv"=INLA::inla.as.dgTMatrix(solve(MeshList$anisotropic_spde$param.inla$M0)) )
  if( "spdeMatricesBarrier" %in% names(Return) ){
    Return[['spdeMatricesBarrier']] = MeshList$barrier_list
    Return[['Barrier_scaling']] = c(0.2, 1)
  }

  #assign(x="X1_itp", value=X1_itp, envir=.GlobalEnv)
  #assign(x="X1_ip", value=array(X1_itp[,1,],dim=dim(X1_itp)[c(1,3)]), envir=.GlobalEnv)

  # Check for NAs
  if( CheckForErrors==TRUE ){
    NoNAs = setdiff( names(Return), c("t_iz","t_yz","c_iz","Network_sz","b_i") )
    if( any(sapply(Return[NoNAs], FUN=function(Array){any(is.na(Array))})==TRUE) ) stop("Please find and eliminate the NA from your inputs")
  }

  # Return
  class(Return) = "make_data"
  return( Return )
}

#' Print data fitted by \code{\link{VAST}}
#'
#' @title Print data
#' @param x Output from \code{\link{make_data}}
#' @param ... Not used
#' @return NULL
#' @method print make_data
#' @export
print.make_data <- function(x, ...)
{
  Samples_iz = data.frame( "b_i"=x$b_i, "a_i"=x$a_i, "c_iz"=x$c_iz, "t_i"=x$t_i, "e_i"=x$e_i, "v_i"=x$v_i, "PredTF_i"=x$PredTF_i )
  cat("make_data(.) result\n")
  cat( paste0("`n_i = `", x$n_i, " samples\n") )
  cat( "`summary(.)` of sampling data\n" )
  print( summary(Samples_iz) )

  invisible(Samples_iz)
}



