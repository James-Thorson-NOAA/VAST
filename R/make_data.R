
#' Build data input for VAST model
#'
#' \code{make_data} builds a tagged list of data inputs used by TMB for running the model
#'
#' @param b_i Sampled biomass for each observation i
#' @param a_i Sampled area for each observation i
#' @param c_iz Category (e.g., species, length-bin) for each observation i
#' @param t_iz Matrix where each row species the time for each observation i (if t_iz is a vector, it is coerced to a matrix with one column; if it is a matrix with two or more columns, it specifies multiple times for each observation, e.g., both year and season)
#' @param e_i Error distribution for each observation i (by default \code{e_i=c_i})
#' @param v_i sampling category (e.g., vessel or tow) associated with overdispersed variation for each observation i (by default \code{v_i=0} for all samples, which will not affect things given the default values for \code{OverdispersionConfig})
#' @param Version a version number;  If missing, defaults to latest version using \code{FishStatsUtils::get_latest_version(package="VAST")}
#' @param FieldConfig a vector of format c("Omega1"=0, "Epsilon1"=10, "Omega2"="AR1", "Epsilon2"=10), where Omega refers to spatial variation, Epsilon refers to spatio-temporal variation, Omega1 refers to variation in encounter probability, and Omega2 refers to variation in positive catch rates, where 0 is off, "AR1" is an AR1 process, and >0 is the number of elements in a factor-analysis covariance
#' @param OverdispersionConfig a vector of format c("eta1"=0, "eta2"="AR1") governing any correlated overdispersion among categories for each level of v_i, where eta1 is for encounter probability, and eta2 is for positive catch rates, where 0 is off, "AR1" is an AR1 process, and >0 is the number of elements in a factor-analysis covariance (by default, c("eta1"=0, "eta2"=0) and this turns off overdispersion)
#' @param ObsModel_ez an optional matrix with two columns where first column specifies the distribution for positive catch rates, and second element specifies the functional form for encounter probabilities
#' \describe{
#'   \item{ObsModel_ez[e,1]=0}{Normal}
#'   \item{ObsModel_ez[e,1]=1}{Lognormal}
#'   \item{ObsModel_ez[e,1]=2}{Gamma}
#'   \item{ObsModel_ez[e,1]=5}{Negative binomial}
#'   \item{ObsModel_ez[e,1]=6}{Conway-Maxwell-Poisson (likely to be very slow)}
#'   \item{ObsModel_ez[e,1]=7}{Poisson (more numerically stable than negative-binomial)}
#'   \item{ObsModel_ez[e,1]=8}{Compound-Poisson-Gamma, where the expected number of individuals is the 1st-component, the expected biomass per individual is the 2nd-component, and SigmaM is the variance in positive catches (likely to be very slow)}
#'   \item{ObsModel_ez[e,1]=9}{Binned-Poisson (for use with REEF data, where 0=0 individual; 1=1 individual; 2=2:10 individuals; 3=>10 individuals)}
#'   \item{ObsModel_ez[e,1]=10}{Tweedie distribution, where epected biomass (lambda) is the product of 1st-component and 2nd-component, variance scalar (phi) is the 1st component, and logis-SigmaM is the power}
#'   \item{ObsModel_ez[e,1]=11}{Zero-inflated Poisson with additional normally-distributed variation overdispersion in the log-intensity of the Poisson distribution}
#'   \item{ObsModel_ez[e,1]=12}{Poisson distribution (not zero-inflated) with log-intensity from the 1st linear predictor, to be used in combination with the Poisson-link delta model for combining multiple data types}
#'   \item{ObsModel_ez[e,1]=13}{Bernoilli distribution using complementary log-log (cloglog) link from the 1st linear predictor, to be used in combination with the Poisson-link delta model for combining multiple data types}
#'   \item{ObsModel_ez[e,1]=14}{Similar to 12, but also including lognormal overdispersion}
#'   \item{ObsModel_ez[e,2]=0}{Conventional delta-model using logit-link for encounter probability and log-link for positive catch rates}
#'   \item{ObsModel_ez[e,2]=1}{Alternative "Poisson-link delta-model" using log-link for numbers-density and log-link for biomass per number}
#'   \item{ObsModel_ez[e,2]=2}{Link function for Tweedie distribution, necessary for \code{ObsModel_ez[e,1]=8} or \code{ObsModel_ez[e,1]=10}}
#'   \item{ObsModel_ez[e,2]=3}{Conventional delta-model, but fixing encounter probability=1 for any year where all samples encounter the species}
#'   \item{ObsModel_ez[e,2]=4}{Poisson-link delta-model, but fixing encounter probability=1 for any year where all samples encounter the species and encounter probability=0 for any year where no samples encounter the species}
#' }
#' @param RhoConfig vector of form c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0) specifying whether either intercepts (Beta1 and Beta2) or spatio-temporal variation (Epsilon1 and Epsilon2) is structured among time intervals (0: each year as fixed effect; 1: each year as random following IID distribution; 2: each year as random following a random walk; 3: constant among years as fixed effect; 4: each year as random following AR1 process);  If missing, assumed to be zero for each element
#' @param VamConfig Options to estimate interactions, containing three slots:
#' \describe{
#'   \item{VamConfig[0]}{selects method for forming interaction matrix; Turn off feature using 0, or I recommend using 2 by default}
#'   \item{VamConfig[1]}{indicates the rank of the interaction matrix, indicating the number of community axes that have regulated dynamics}
#'   \item{VamConfig[2]}{Indicates whether interactions occur before spatio-temporal variation (\code{VamConfig[2]=0}) or after \code{VamConfig[2]=1}}
#' }
#' @param spatial_list tagged list of locatoinal information from , i.e., from \code{FishStatsUtils::make_spatial_info}
#' @param PredTF_i OPTIONAL, whether each observation i is included in the likelihood (PredTF_i[i]=0) or in the predictive probability (PredTF_i[i]=1)
#' @param X_gtp array of density covariates for each extrapolation-grid cell g, time t, and covariate p; if missing, assumed to not include covariates
#' @param X_itp array of density covariates for each observation i, time t, and covariate p
#' @param Xconfig_zcp OPTIONAL, 3D array of settings for each dynamic density covariate, where the first dimension corresponds to 1st or 2nd linear predictors, second dimension corresponds to model category, and third dimension corresponds to each density covariate
#' \describe{
#'   \item{Xconfig_zcp[z,c,p]=0}{\code{X_itp[,,p]} has no effect on linear predictor z for category c}
#'   \item{Xconfig_zcp[z,c,p]=1}{\code{X_itp[,,p]} has a linear effect on linear predictor z for category c}
#'   \item{Xconfig_zcp[z,c,p]=2}{\code{X_itp[,,p]} has a spatially varying, zero-centered linear effect on linear predictor z for category c}
#'   \item{Xconfig_zcp[z,c,p]=3}{\code{X_itp[,,p]} has a spatially varying linear effect on linear predictor z for category c}
#' }
#' @param Q_ik matrix of catchability covariates (e.g., measured variables affecting catch rates but not caused by variation in species density) for each observation i
#' @param Aniso whether to assume isotropy (Aniso=0) or geometric anisotropy (Aniso=1)
#' @param Expansion_cz matrix specifying how densities are expanded when calculating annual indices, with a row for each category \code{c} and two columns.  The first column specifies whether to calculate annual index for category \code{c} as the weighted-sum across density estimates, where density is weighted by area ("area-weighted expansion", \code{Expansion[c,1]=0}, the default) or where density is weighted by the expanded value for another category ("abundance weighted expansion" \code{Expansion[c1,1]=1}).  The 2nd column is only used when \code{Expansion[c1,1]=1}, and specifies the category to use for abundance-weighted expansion, where \code{Expansion[c1,2]=c2} and \code{c2} must be lower than \code{c1}.
#' @param F_ct matrix of fishing mortality for each category c and year t (only feasible when using a Poisson-link delta model and specifying temporal structure on intercepts, and mainly interpretable when species interactions via VamConfig)
#' @param t_yz matrix specifying combination of levels of \code{t_iz} to use when calculating different indices of abundance or range shifts
#' @param Options a vector of form c('SD_site_logdensity'=FALSE,'Calculate_Range'=FALSE,'Calculate_effective_area'=FALSE,'Calculate_Cov_SE'=FALSE,'Calculate_Synchrony'=FALSE,'Calculate_proportion'=FALSE), where Calculate_Range=1 turns on calculation of center of gravity, and Calculate_effective_area=1 turns on calculation of effective area occupied
#' @param yearbounds_zz matrix with two columns, giving first and last years for defining one or more periods (rows) used to calculate changes in synchrony over time (only used if \code{Options['Calculate_Synchrony']=1})
#' @param CheckForErrors whether to check for errors in input (NOTE: when CheckForErrors=TRUE, the function will throw an error if it detects a problem with inputs.  However, failing to throw an error is no guaruntee that the inputs are all correct)
#' @param ... interface to pass deprecated inputs, included for backwards compatibility with previous versions which specified elements of \code{spatial_list} individually instead of as a single object

#' @return Object of class \code{make_data}, containing inputs to function \code{VAST::Build_TMB_Fn()}

#' @export
make_data <-
function( b_i, a_i, c_iz, t_iz, e_i=c_iz[,1], v_i=rep(0,length(b_i)),
  FieldConfig, spatial_list, ObsModel_ez=c("PosDist"=1,"Link"=0),
  OverdispersionConfig=c("eta1"=0,"eta2"=0), RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0),
  VamConfig=c("Method"=0,"Rank"=0,"Timing"=0), Aniso=TRUE, PredTF_i=rep(0,length(b_i)),
  Xconfig_zcp=NULL, X_gtp=NULL, X_itp=NULL,
  Q_ik=NULL, Network_sz=NULL, F_ct=NULL, F_init=1,
  t_yz=NULL, CheckForErrors=TRUE, yearbounds_zz=NULL,
  Options=c(), Expansion_cz=NULL,
  Version=FishStatsUtils::get_latest_version(package="VAST"), ... ){

  # Deprecated inputs for backwards compatibility in transition from Version < 8.0.0 to >= 8.0.0
  deprecated_inputs = list( ... )
  X_xtp = deprecated_inputs[["X_xtp"]]
  if( missing(spatial_list) ){
    warning("Consider changing use of `make_data` to include `spatial_list` as input")
    a_xl = a_gl = deprecated_inputs[["a_xl"]]
    MeshList = deprecated_inputs[["MeshList"]]
    GridList = deprecated_inputs[["GridList"]]
    Method = deprecated_inputs[["Method"]]
    s_i = deprecated_inputs[["s_i"]]
  }else{
    MeshList = spatial_list[["MeshList"]]
    GridList = spatial_list[["GridList"]]
    Method = spatial_list[["Method"]]
    a_xl = a_gl = spatial_list[["a_gl"]]
    s_i = spatial_list[["knot_i"]] - 1
  }

  # Deprecated inputs for backwards compatibility with earlier versions
  X_xj = deprecated_inputs[["X_xj"]]

  # Specify default values for `Options`
  Options2use = c('SD_site_density'=FALSE, 'SD_site_logdensity'=FALSE, 'Calculate_Range'=FALSE, 'SD_observation_density'=FALSE, 'Calculate_effective_area'=FALSE,
    'Calculate_Cov_SE'=FALSE, 'Calculate_Synchrony'=FALSE, 'Calculate_Coherence'=FALSE, 'Calculate_proportion'=FALSE, 'normalize_GMRF_in_CPP'=TRUE,
    'Calculate_Fratio'=FALSE, 'Estimate_B0'=FALSE, 'Project_factors'=FALSE, 'treat_nonencounter_as_zero'=FALSE, 'simulate_random_effects'=TRUE )

  # Replace defaults for `Options` with provided values (if any)
  for( i in seq_along(Options) ){
    if(tolower(names(Options)[i]) %in% tolower(names(Options2use))){
      Options2use[[match(tolower(names(Options)[i]),tolower(names(Options2use)))]] = Options[[i]]
    }
  }

  # Check for backwards-compatibility issues
  if( FishStatsUtils::convert_version_name(Version) >= FishStatsUtils::convert_version_name("VAST_v8_0_0") ){
    if( !is.null(X_xtp) ){
      stop("`X_xtp` is not used in version >= 8.0.0. If you'd like to specify covariates using input `X_xtp` please use `Version='VAST_v7_0_0'`")
    }
  }
  if( FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v7_0_0") ){
    if( !is.null(X_gtp) | !is.null(X_itp) ){
      stop("`X_gtp` and `X_itp` are not used in version <= 7.0.0. If you'd like to specify covariates using input `X_gtp` and `X_itp` please use `Version='VAST_v8_0_0'` or higher")
    }
  }

  # Adds intercept defaults to FieldConfig if missing
  if( is.vector(FieldConfig) && length(FieldConfig)==4 ){
    FieldConfig = rbind( matrix(FieldConfig,ncol=2,dimnames=list(c("Omega","Epsilon"),c("Component_1","Component_2"))), "Beta"=c("Beta1"="IID","Beta2"="IID") )
  }else{
    if( !is.matrix(FieldConfig) || !all(dim(FieldConfig)==c(3,2)) ){
      stop("`FieldConfig` has the wrong dimensions in `make_data`")
    }else{
      dimnames(FieldConfig) = list( c("Omega","Epsilon","Beta"), c("Component_1","Component_2") )
    }
  }

  # Rescale tprime_iz to start at 0
  tprime_iz = t_iz - min(t_iz,na.rm=TRUE)

  # Coerce tprime_iz to be a matrix
  if( !is.matrix(tprime_iz) ) tprime_iz = matrix(tprime_iz,ncol=1)

  # Increment first tprime_iz if t=0 corresponds to B0
  if( Options2use[12]==1 ){
    tprime_iz = tprime_iz + 1
    F_ct = cbind( 0, F_ct )
  }

  # Coerce c_iz to be a matrix
  if( !is.matrix(c_iz) ) c_iz = matrix(c_iz,ncol=1)

  # Determine dimensions
  n_t = max(tprime_iz,na.rm=TRUE) + 1
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
    Index = list( factor(c_iz[,1],levels=0:max(c_iz[,1])), factor(tprime_iz[,1],levels=0:max(tprime_iz[,1])) )
    Num_ct = tapply( b_i, INDEX=Index, FUN=function(vec){sum(vec>0)} )
    Num_ct = ifelse( is.na(Num_ct), 0, Num_ct )
    b_i = ifelse( Num_ct[cbind(as.numeric(Index[[1]]),as.numeric(Index[[2]]))]==0, NA, b_i )
  }

  # Covariates and defaults
  if( is.null(X_xj) ){
    X_xj = matrix(0, nrow=n_x, ncol=1)
  }else{
    if( !is.array(X_xj) || !(dim(X_xj)[1]==n_x) ){
      stop("`X_xj` has wrong dimensions")
    }
  }
  if( is.null(X_xtp) ){
    X_xtp = array(0, dim=c(n_x,n_t,1))
  }else{
    if( !is.array(X_xtp) || !(all(dim(X_xtp)[1:2]==c(n_x,n_t))) ){
      stop("`X_xtp` has wrong dimensions")
    }
  }
  if( is.null(X_gtp) ){
    X_gtp = array(0, dim=c(n_g,n_t,1))
  }else{
    if( !is.array(X_gtp) || !(all(dim(X_gtp)[1:2]==c(n_g,n_t))) ){
      stop("`X_gtp` has wrong dimensions")
    }
  }
  if( is.null(X_itp) ){
    X_itp = array(0, dim=c(n_i,n_t,1))
  }else{
    if( !is.array(X_itp) || !(all(dim(X_itp)[1:2]==c(n_i,n_t))) ){
      stop("`X_itp` has wrong dimensions")
    }
  }
  if( is.null(Q_ik) ){
    Q_ik = matrix(0, nrow=n_i, ncol=1)
  }else{
    if( !is.array(Q_ik) || !(all(dim(Q_ik)[1]==c(n_i))) ){
      stop("`Q_ik` has wrong dimensions")
    }
  }
  if( is.null(yearbounds_zz)){
    yearbounds_zz = matrix(c(0,n_t-1),nrow=1)
  }else{
    if( !is.array(yearbounds_zz) || !(all(dim(yearbounds_zz)[2]==2)) ){
      stop("`yearbounds_zz` has wrong dimensions")
    }
  }
  if( is.null(t_yz) ){
    t_yz = matrix(0:max(tprime_iz[,1],na.rm=TRUE), ncol=1)
    for( cI in seq(2,ncol(tprime_iz),length=ncol(tprime_iz)-1)) t_yz = cbind(t_yz, min(tprime_iz[,cI],na.rm=TRUE))
  }
  n_j = ncol(X_xj)
  if( FishStatsUtils::convert_version_name(Version) >= FishStatsUtils::convert_version_name("VAST_v8_0_0") ){
    n_p = dim(X_gtp)[3]
  }else{
    n_p = dim(X_xtp)[3]
  }
  n_k = ncol(Q_ik)
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
  if( is.null(Xconfig_zcp) ){
    Xconfig_zcp = array(1, dim=c(2,n_c,n_p))
  }else{
    if( !is.array(Xconfig_zcp) || !(all(dim(Xconfig_zcp)==c(2,n_c,n_p))) ){
      stop("`Xconfig_zcp` has wrong dimensions")
    }
    if( !all(Xconfig_zcp %in% c(0,1,2,3)) ){
      stop("`Xconfig_zcp` has some wrong element(s)")
    }
  }

  # Translate FieldConfig from input formatting to CPP formatting
  FieldConfig_input = array(NA, dim=dim(FieldConfig), dimnames=dimnames(FieldConfig) )
  g = function(mat) suppressWarnings( array(as.numeric(mat),dim=dim(mat)) )
  FieldConfig_input[] = ifelse( FieldConfig=="AR1", 0, FieldConfig_input)
  FieldConfig_input[] = ifelse( FieldConfig=="IID", -2, FieldConfig_input)
  FieldConfig_input[] = ifelse( !is.na(g(FieldConfig)) & g(FieldConfig)>0 & g(FieldConfig)<=n_c, g(FieldConfig), FieldConfig_input)
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
      Z_gm = spatial_list$loc_g
    }
    message( "Calculating range shift for stratum #1:",colnames(a_gl[1]))
  }

  ###################
  # Check for bad data entry
  ###################

  if( CheckForErrors==TRUE ){
    if( ncol(ObsModel_ez)!=2 | nrow(ObsModel_ez)!=n_e ) stop("Wrong dimensions for ObsModel_ez")
    if( !is.matrix(a_gl) | !is.matrix(X_xj) | !is.matrix(Q_ik) ) stop("a_gl, X_xj, and Q_ik should be matrices")
    if( !is.array(X_xtp) | !is.array(X_gtp) | !is.array(X_itp) ) stop( "X_xtp, X_gtp, and X_itp should be arrays")
    if( (max(s_i)-1)>n_x | min(s_i)<0 ) stop("s_i exceeds bounds in MeshList")
    if( any(a_i<=0) ) stop("a_i must be greater than zero for all observations, and at least one value of a_i is not")
    # Warnings about all positive or zero
    Prop_nonzero = tapply( b_i, INDEX=list(tprime_iz[,1],c_iz[,1]), FUN=function(vec){mean(vec>0)} )
    if( Options2use[12]==1 ){
      Prop_nonzero = Prop_nonzero[-1,]
    }
    if( any(ObsModel_ez[,2] %in% c(0,1)) ){
      if( RhoConfig[1]==0 ){
        if( any(!is.na(Prop_nonzero) & (Prop_nonzero==1)) ){
          print( Prop_nonzero )
          stop("Some years and/or categories have 100% encounters, and this requires either temporal structure of a different link-function")
        }
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
    if( ncol(t_yz)!=ncol(tprime_iz) ) stop("t_yz and tprime_iz must have same number of columns")
    if( n_c!=length(unique(na.omit(as.vector(c_iz)))) ) stop("n_c doesn't equal the number of levels in c_i")
    if( any(X_xj!=0) ) stop("X_xj is deprecated, please use X_xtp to specify static or dynamic density covariates (which by default have constant effect among years but differ among categories)")
    if( any(ObsModel_ez[,1]==9) & !all(b_i%in%0:3) ) stop("If using 'ObsModel_ez[e,1]=9', all 'b_i' must be in {0,1,2,3}")
    if( length(unique(ObsModel_ez[,2]))>1 ) stop("All `ObsModel_ez[,2]` must have the same value")
    if( any(OverdispersionConfig>0) & length(unique(v_i))==1 ) stop("It doesn't make sense to use use `OverdispersionConfig` when using only one level of `v_i`")
    if( any(ObsModel_ez[,1] %in% c(12,13,14)) ){
      if( any(ObsModel_ez[,2] != 1) ) stop("If using `ObsModel_ez[e,1]` in {12,13,14} then must use `ObsModel_ez[e,2]=1`")
      if( !any(ObsModel_ez[,1] %in% c(0,1,2,3)) ) stop("Using `ObsModel_ez[e,1]` in {12,13,14} is only intended when combining data with biomass-sampling data")
    }
    if( all(b_i>0) & all(ObsModel_ez[,1]==0) & !all(FieldConfig_input[1:2,1]==-1) ) stop("All data are positive and using a conventional delta-model, so please turn off `Omega1` and `Epsilon1` terms")
    if( !(all(ObsModel_ez[,1] %in% c(0,1,2,5,6,7,8,9,10,11,12,13,14))) ) stop("Please check `ObsModel_ez[,1]` input")
    if( !(all(ObsModel_ez[,2] %in% c(0,1,2,3,4))) ) stop("Please check `ObsModel_ez[,2]` input")
    if( !all(RhoConfig[1]%in%c(0,1,2,3,4)) | !all(RhoConfig[2]%in%c(0,1,2,3,4,6)) | !all(RhoConfig[3]%in%c(0,1,2,4,5)) | !all(RhoConfig[4]%in%c(0,1,2,4,5,6)) ) stop("Check `RhoConfig` inputs")
    if( any(is.na(X_xtp)) ) stop("Some `X_xtp` is NA, and this is not allowed")
    if( any(is.na(X_gtp)) ) stop("Some `X_gtp` is NA, and this is not allowed")
    if( any(is.na(X_itp)) ) stop("Some `X_itp` is NA, and this is not allowed")
    if( n_c==1 && !all(FieldConfig_input %in% c(-2,-1,1)) ) stop("If using a univariate model, `FieldConfig` must be 0, 1, or `IID` for all entries")
  }

  # Check for wrong dimensions
  if( CheckForErrors==TRUE ){
    if( any(c(length(b_i),length(a_i),nrow(c_iz),length(s_i),nrow(tprime_iz),length(v_i),length(PredTF_i))!=n_i) ) stop("b_i, a_i, c_i, s_i, v_i, or tprime_i doesn't have length n_i")
    if( nrow(a_gl)!=n_x | ncol(a_gl)!=n_l ) stop("a_xl has wrong dimensions")
    if( nrow(X_xj)!=n_x | ncol(X_xj)!=n_j ) stop("X_xj has wrong dimensions")
    if( nrow(Q_ik)!=n_i | ncol(Q_ik)!=n_k ) stop("Q_ik has wrong dimensions")
    if( FishStatsUtils::convert_version_name(Version) >= FishStatsUtils::convert_version_name("VAST_v8_0_0") ){
      if( dim(X_gtp)[1]!=n_g | dim(X_gtp)[2]!=n_t | dim(X_gtp)[3]!=n_p ) stop("X_gtp has wrong dimensions")
      if( dim(X_itp)[1]!=n_i | dim(X_itp)[2]!=n_t | dim(X_itp)[3]!=n_p ) stop("X_itp has wrong dimensions")
    }else{
      if( dim(X_xtp)[1]!=n_x | dim(X_xtp)[2]!=n_t | dim(X_xtp)[3]!=n_p ) stop("X_xtp has wrong dimensions")
    }
    if( ncol(c_iz)>1 & any(ObsModel_ez[,2]!=1) ) stop("Using multiple columnns in `c_iz` only makes sense using a Poisson-link delta model via `ObsModel[2]=1`")
    if( nrow(F_ct)!=n_c | ncol(F_ct)!=n_t ) stop("F_ct has wrong dimensions")
  }

  ###################
  # Check for incompatibilities amongst versions
  ###################

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

  ###################
  # Check for incompatible settings
  ###################

  # Seasonal models and intercepts
  if( ncol(tprime_iz)>=2 & any(RhoConfig[1:2]!=0) ){
    stop("Temporal structure on intercepts is not implemented for seasonal models")
  }
  if( ncol(tprime_iz)>=2 & any(VamConfig[1]!=0) ){
    stop("Species interactions are not implemented for seasonal models")
  }
  if( (FieldConfig_input[2,1]==(-1) & RhoConfig[3]!=0) | (FieldConfig_input[2,2]==(-1) & RhoConfig[4]!=0) ){
    stop("Spatio-temporal variation is turned off for a component with temporal structure, and this combination doesn't make sense")
  }

  # Prohibitively slow
  if( !is.null(spatial_list$fine_scale) && spatial_list$fine_scale==TRUE ){
    if( Options2use['SD_site_density']==TRUE | Options2use['SD_site_logdensity']==TRUE ){
      warning("'SD_site_density' and 'SD_site_logdensity' are very slow when using `fine_scale=TRUE`")
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
  if( any(RhoConfig[1:2]!=0) & any(ObsModel_ez[,2]==3) ){
    stop( "RhoConfig[1:2] must be 0 when using ObsModel[2]=3:  Other options are not coded to work together" )
  }
  if( any(RhoConfig[1:2]!=0) & any(ObsModel_ez[,2]==4) ){
    stop( "RhoConfig[1:2] must be 0 when using ObsModel[2]=4:  Other options are not coded to work together" )
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
  if( VamConfig[2]==0 & VamConfig[1]!=0 ){
    VamConfig[1] = 0
    message("Using interactions with zero rank (`VamConfig[2]==0`), so turning off interactions (`VamConfig[1]=0`)")
  }

  # Warning messages
  if( n_c>1 & any(FieldConfig_input==1)){
    warning( "Using 1 factor for more than one category:  Please note that this is non-standard, and it is more common to use multiple factors (often as many as the number of categories)" )
  }
  SD_p = apply( X_xtp, MARGIN=3, FUN=sd )
  if( any(SD_p>3) ){
    warning( "I highly recommend that you standardize each density covariate `X_xtp` to have a low standard deviation, to avoid numerical under/over-flow" )
  }

  ###################
  # Output tagged list
  ###################

  # CMP_xmax should be >100 and CMP_breakpoint should be 1 for Tweedie model
  Options_vec = c("Aniso"=Aniso, "R2_interpretation"=0, "Rho_beta1_TF"=ifelse(RhoConfig[1]%in%c(1,2,4),1,0), "Rho_beta2_TF"=ifelse(RhoConfig[2]%in%c(1,2,4),1,0), "AreaAbundanceCurveTF"=0, "CMP_xmax"=200, "CMP_breakpoint"=1, "Method"=switch(Method,"Mesh"=0,"Grid"=1,"Spherical_mesh"=0,"Stream_network"=2), "Include_F"=ifelse(all(F_ct==0),0,F_init) )
  Return = NULL
  if(Version%in%c("VAST_v1_1_0","VAST_v1_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=tprime_iz[,1], "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
  }
  if(Version%in%c("VAST_v1_4_0","VAST_v1_3_0","VAST_v1_2_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_f_input"=OverdispersionConfig_input[1], "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=tprime_iz[,1], "v_i"=match(v_i,sort(unique(v_i)))-1, "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
  }
  if(Version%in%c("VAST_v1_6_0","VAST_v1_5_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_f_input"=OverdispersionConfig_input[1], "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=tprime_iz[,1], "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
  }
  if(Version%in%c("VAST_v1_7_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=tprime_iz[,1], "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
  }
  if(Version%in%c("VAST_v1_8_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=tprime_iz[,1], "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v1_9_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=tprime_iz[,1], "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v2_8_0","VAST_v2_7_0","VAST_v2_6_0","VAST_v2_5_0","VAST_v2_4_0","VAST_v2_3_0","VAST_v2_2_0","VAST_v2_1_0","VAST_v2_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_iz"=tprime_iz, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v3_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "e_i"=e_i, "s_i"=s_i, "t_iz"=tprime_iz, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v4_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=tprime_iz, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v4_4_0","VAST_v4_3_0","VAST_v4_2_0","VAST_v4_1_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=n_p, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "include_data"=TRUE, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=tprime_iz, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v5_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=n_p, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "include_data"=TRUE, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=tprime_iz, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v5_1_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=n_p, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "include_data"=TRUE, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=tprime_iz, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v5_2_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=n_p, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "include_data"=TRUE, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=tprime_iz, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v5_4_0","VAST_v5_3_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=n_p, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "include_data"=TRUE, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=tprime_iz, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v5_5_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=n_p, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_list"=list("Options_vec"=Options_vec,"Options"=Options2use,"yearbounds_zz"=yearbounds_zz,"Expansion_cz"=Expansion_cz), "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "include_data"=TRUE, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=tprime_iz, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v6_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=n_p, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_list"=list("Options_vec"=Options_vec,"Options"=Options2use,"yearbounds_zz"=yearbounds_zz,"Expansion_cz"=Expansion_cz), "FieldConfig"=as.vector(FieldConfig_input[1:2,]), "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "Xconfig_zcp"=Xconfig_zcp, "include_data"=TRUE, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=tprime_iz, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v7_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=n_p, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_list"=list("Options_vec"=Options_vec,"Options"=Options2use,"yearbounds_zz"=yearbounds_zz,"Expansion_cz"=Expansion_cz), "FieldConfig"=FieldConfig_input, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "Xconfig_zcp"=Xconfig_zcp, "include_data"=TRUE, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=tprime_iz, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_gl, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v8_2_0","VAST_v8_1_0","VAST_v8_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x,n_x)[Options_vec['Method']+1], "n_g"=n_g, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=n_p, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_gm), "Options_list"=list("Options_vec"=Options_vec,"Options"=Options2use,"yearbounds_zz"=yearbounds_zz,"Expansion_cz"=Expansion_cz), "FieldConfig"=FieldConfig_input, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "Xconfig_zcp"=Xconfig_zcp, "include_data"=TRUE, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "t_iz"=tprime_iz, "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_gl"=a_gl, "X_itp"=X_itp, "X_gtp"=X_gtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_gm"=Z_gm, "F_ct"=F_ct, "parent_s"=Network_sz[,'parent_s']-1, "child_s"=Network_sz[,'child_s']-1, "dist_s"=Network_sz[,'dist_s'], "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2, "Ais_ij"=cbind(spatial_list$A_is@i,spatial_list$A_is@j), "Ais_x"=spatial_list$A_is@x, "Ags_ij"=cbind(spatial_list$A_gs@i,spatial_list$A_gs@j), "Ags_x"=spatial_list$A_gs@x )
  }
  if( is.null(Return) ) stop("`Version` provided does not match the list of possible values")
  if( "spde" %in% names(Return) ) Return[['spde']] = MeshList$isotropic_spde$param.inla[c("M0","M1","M2")]
  if( "spde_aniso" %in% names(Return) ) Return[['spde_aniso']] = list("n_s"=MeshList$anisotropic_spde$n.spde, "n_tri"=nrow(MeshList$anisotropic_mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$anisotropic_spde$param.inla$M0, "G0_inv"=INLA::inla.as.dgTMatrix(solve(MeshList$anisotropic_spde$param.inla$M0)) )

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
  Samples_iz = data.frame( "b_i"=x$b_i, "a_i"=x$a_i, "c_iz"=x$c_iz, "t_iz"=x$t_iz, "e_i"=x$e_i, "v_i"=x$v_i, "PredTF_i"=x$PredTF_i )
  cat("make_data(.) result\n")
  cat( paste0("`n_i = `", x$n_i, " samples\n") )
  cat( "`summary(.)` of sampling data\n" )
  print( summary(Samples_iz) )

  invisible(Samples_iz)
}



