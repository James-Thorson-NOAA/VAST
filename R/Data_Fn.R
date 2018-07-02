
#' Build data input for VAST model
#'
#' \code{Data_Fn} builds a tagged list of data inputs used by TMB for running the model
#'
#' @param Version a version number (see example for current default).
#' @param FieldConfig a vector of format c("Omega1"=0, "Epsilon1"=10, "Omega2"="AR1", "Epsilon2"=10), where Omega refers to spatial variation, Epsilon refers to spatio-temporal variation, Omega1 refers to variation in encounter probability, and Omega2 refers to variation in positive catch rates, where 0 is off, "AR1" is an AR1 process, and >0 is the number of elements in a factor-analysis covariance
#' @param OverdispersionConfig OPTIONAL, a vector of format c("eta1"=0, "eta2"="AR1") governing any correlated overdispersion among categories for each level of v_i, where eta1 is for encounter probability, and eta2 is for positive catch rates, where 0 is off, "AR1" is an AR1 process, and >0 is the number of elements in a factor-analysis covariance
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
#'   \item{ObsModel_ez[e,2]=1}{Alternative delta-model using log-link for numbers-density and log-link for biomass per number}
#'   \item{ObsModel_ez[e,2]=2}{Link function for Tweedie distribution, necessary for \code{ObsModel_ez[e,1]=8} or \code{ObsModel_ez[e,1]=10}}
#'   \item{ObsModel_ez[e,2]=3}{Conventional delta-model, but fixing encounter probability=1 for any year where all samples encounter the species}
#' }
#' @param VamConfig Options to estimate interactions, where first slot selects method for forming interaction matrix and second indicates rank
#' @param b_i Sampled biomass for each observation i
#' @param a_i Sampled area for each observation i
#' @param c_iz Category (e.g., species, length-bin) for each observation i
#' @param s_i Spatial knot (e.g., grid cell) for each observation i
#' @param t_iz Matrix where each row species the time for each observation i (if t_iz is a vector, it is coerced to a matrix with one column; if it is a matrix with two or more columns, it specifies multiple times for each observation, e.g., both year and season)
#' @param a_xl Area associated with each knot
#' @param MeshList, tagged list representing location information for the SPDE mesh hyperdistribution, i.e., from \code{SpatialDeltaGLMM::Spatial_Information_Fn}
#' @param GridList, tagged list representing location information for the 2D AR1 grid hyperdistribution, i.e., from \code{SpatialDeltaGLMM::Spatial_Information_Fn}
#' @param Method, character (either "Mesh" or "Grid") specifying hyperdistribution (Default="Mesh")
#' @param v_i OPTIONAL, sampling category (e.g., vessel or tow) associated with overdispersed variation for each observation i
#' @param e_i Error distribution for each observation i (by default \code{e_i=c_i})
#' @param PredTF_i OPTIONAL, whether each observation i is included in the likelihood (PredTF_i[i]=0) or in the predictive probability (PredTF_i[i]=1)
#' @param X_xj OPTIONAL, matrix of static density covariates (e.g., measured variables affecting density, as used when interpolating density for calculating an index of abundance)
#' @param X_xtp OPTIONAL, array of dynamic (varying among time intervals) density covariates
#' @param Q_ik OPTIONAL, matrix of catchability covariates (e.g., measured variables affecting catch rates but not caused by variation in species density) for each observation i
#' @param Aniso OPTIONAL, whether to assume isotropy (Aniso=0) or geometric anisotropy (Aniso=1)
#' @param RhoConfig OPTIONAL, vector of form c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0) specifying whether either intercepts (Beta1 and Beta2) or spatio-temporal variation (Epsilon1 and Epsilon2) is structured among time intervals (0: each year as fixed effect; 1: each year as random following IID distribution; 2: each year as random following a random walk; 3: constant among years as fixed effect; 4: each year as random following AR1 process)
#' @param t_yz OPTIONAL, matrix specifying combination of levels of \code{t_iz} to use when calculating different indices of abundance or range shifts
#' @param Options OPTIONAL, a vector of form c('SD_site_logdensity'=0,'Calculate_Range'=0,'Calculate_effective_area'=0,'Calculate_Cov_SE'=0,'Calculate_Synchrony'=0,'Calculate_proportion'=0), where Calculate_Range=1 turns on calculation of center of gravity, and Calculate_effective_area=1 turns on calculation of effective area occupied
#' @param yearbounds_zz OPTIONAL, matrix with two columns, giving first and last years for defining one or more periods (rows) used to calculate changes in synchrony over time (only used if \code{Options['Calculate_Synchrony']=1})
#' @param CheckForErrors OPTIONAL, whether to check for errors in input (NOTE: when CheckForErrors=TRUE, the function will throw an error if it detects a problem with inputs.  However, failing to throw an error is no guaruntee that the inputs are all correct)

#' @return Tagged list containing inputs to function \code{VAST::Build_TMB_Fn()}

#' @export
Data_Fn <-
function( Version, FieldConfig, OverdispersionConfig=c("eta1"=0,"eta2"=0), ObsModel_ez=c("PosDist"=1,"Link"=0), VamConfig=c("Method"=0,"Rank"=0,"Timing"=0),
  b_i, a_i, c_iz, s_i, t_iz, a_xl, MeshList, GridList, Method, v_i=rep(0,length(b_i)), e_i=c_iz[,1],
  PredTF_i=rep(0,length(b_i)), X_xj=NULL, X_xtp=NULL, Q_ik=NULL, Aniso=1,
  RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0), t_yz=NULL, CheckForErrors=TRUE, yearbounds_zz=NULL,
  Options=c('SD_site_logdensity'=0,'Calculate_Range'=0,'Calculate_effective_area'=0,'Calculate_Cov_SE'=0,'Calculate_Synchrony'=0,'Calculate_proportion'=0) ){

  # Specify default values for `Options`
  Options2use = c('SD_site_density'=0, 'SD_site_logdensity'=0, 'Calculate_Range'=0, 'SD_observation_density'=0, 'Calculate_effective_area'=0,
    'Calculate_Cov_SE'=0, 'Calculate_Synchrony'=0, 'Calculate_Coherence'=0, 'Calculate_proportion'=0, 'normalize_GMRF_in_CPP'=TRUE)

  # Replace defaults for `Options` with provided values (if any)
  for( i in 1:length(Options)){
    if(tolower(names(Options)[i]) %in% tolower(names(Options2use))){
      Options2use[[match(tolower(names(Options)[i]),tolower(names(Options2use)))]] = Options[[i]]
    }
  }

  # Coerce t_iz to be a matrix
  if( !is.matrix(t_iz) ) t_iz = matrix(t_iz,ncol=1)

  # Coerce t_iz to be a matrix
  if( !is.matrix(c_iz) ) c_iz = matrix(c_iz,ncol=1)

  # Determine dimensions
  n_t = max(t_iz,na.rm=TRUE) - min(t_iz,na.rm=TRUE) + 1
  n_c = max(c_iz,na.rm=TRUE) + 1
  n_e = max(e_i) + 1
  n_v = length(unique(v_i))   # If n_v=1, then turn off overdispersion later
  n_i = length(b_i)
  n_x = nrow(a_xl)
  n_l = ncol(a_xl)

  # Coerce ObsModel_ez to be a matrix
  if( !is.matrix(ObsModel_ez) ) ObsModel_ez = matrix( ObsModel_ez, ncol=2, nrow=n_e, byrow=TRUE )

  # Covariates and defaults
  if( is.null(X_xj) ) X_xj = matrix(0, nrow=n_x, ncol=1)
  if( is.null(X_xtp) ) X_xtp = array(0, dim=c(n_x,n_t,1))
  if( is.null(Q_ik) ) Q_ik = matrix(0, nrow=n_i, ncol=1)
  if( is.null(yearbounds_zz)) yearbounds_zz = matrix(c(0,n_t-1),nrow=1)
  if( is.null(t_yz) ){
    t_yz = matrix(0:(max(t_iz[,1])-min(t_iz[,1])), ncol=1)
    for( cI in seq(2,ncol(t_iz),length=ncol(t_iz)-1)) t_yz = cbind(t_yz, min(t_iz[,cI],na.rm=TRUE)-min(t_iz[,1]))
  }
  n_j = ncol(X_xj)
  n_p = dim(X_xtp)[3]
  n_k = ncol(Q_ik)
  n_y = nrow(t_yz)

  # Translate FieldConfig from input formatting to CPP formatting
  FieldConfig_input = rep(NA, length(FieldConfig))
  names(FieldConfig_input) = names(FieldConfig)
  g = function(vec) suppressWarnings(as.numeric(vec))
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
    Z_xm = matrix(0, nrow=nrow(a_xl), ncol=ncol(a_xl) ) # Size so that it works for Version 3g-3j
  }else{
    Z_xm = MeshList$loc_x
    message( "Calculating range shift for stratum #1:",colnames(a_xl[1]))
  }

  # Check for bad data entry
  if( CheckForErrors==TRUE ){
    if( ncol(ObsModel_ez)!=2 | nrow(ObsModel_ez)!=n_e ) stop("Wrong dimensions for ObsModel_ez")
    if( !is.matrix(a_xl) | !is.matrix(X_xj) | !is.matrix(Q_ik) ) stop("a_xl, X_xj, and Q_ik should be matrices")
    if( (max(s_i)-1)>n_x | min(s_i)<0 ) stop("s_i exceeds bounds in MeshList")
    if( any(a_i<=0) ) stop("a_i must be greater than zero for all observations, and at least one value of a_i is not")
    # Warnings about all positive or zero
    Prop_nonzero = tapply( b_i, INDEX=list(t_iz[,1],c_iz[,1]), FUN=function(vec){mean(vec>0)} )
    if( any(Prop_nonzero==0|Prop_nonzero==1) & any(ObsModel_ez[1,]==0) ){
      print( Prop_nonzero )
      stop("Some years and/or categories have either all or no encounters, and this is not permissible when ObsModel_ez[,'Link']=0")
    }
    if( length(OverdispersionConfig)!=2 ) stop("length(OverdispersionConfig)!=2")
    if( ncol(yearbounds_zz)!=2 ) stop("yearbounds_zz must have two columns")
    if( Options2use['Calculate_Coherence']==1 & any(ObsModel_ez[,2]==0) ) stop("Calculating coherence only makes sense when 'ObsModel_ez[,2]=1'")
    if( any(yearbounds_zz<0) | any(yearbounds_zz>=max(n_t)) ) stop("yearbounds_zz exceeds bounds for modeled years")
    if( ncol(t_yz)!=ncol(t_iz) ) stop("t_yz and t_iz must have same number of columns")
    if( n_c!=length(unique(na.omit(as.vector(c_iz)))) ) stop("n_c doesn't equal the number of levels in c_i")
    if( any(X_xj!=0) ) stop("X_xj is deprecated, please use X_xtp to specify static or dynamic density covariates (which by default have constant effect among years but differ among categories)")
    if( any(ObsModel_ez[,1]==9) & !all(b_i%in%0:3) ) stop("If using 'ObsModel_ez[e,1]=9', all 'b_i' must be in {0,1,2,3}")
    if( length(unique(ObsModel_ez[,2]))>1 ) stop("All `ObsModel_ez[,2]` must have the same value")
    if( any(OverdispersionConfig>0) & length(unique(v_i))==1 ) stop("It doesn't make sense to use use `OverdispersionConfig` when using only one level of `v_i`")
    if( any(ObsModel_ez[,1] %in% c(12,13,14)) ){
      if( any(ObsModel_ez[,2] != 1) ) stop("If using `ObsModel_ez[e,1]` in {12,13,14} then must use `ObsModel_ez[e,2]=1`")
      if( !any(ObsModel_ez[,1] %in% c(0,1,2,3)) ) stop("Using `ObsModel_ez[e,1]` in {12,13,14} is only intended when combining data with biomass-sampling data")
    }
    if( all(b_i>0) & any(FieldConfig_input[1:2]>0) & any(ObsModel_ez[1,]==0) ) stop("All data are positive, so please turn off `Omega1` and `Epsilon1` terms")
  }

  # Check for wrong dimensions
  if( CheckForErrors==TRUE ){
    if( any(c(length(b_i),length(a_i),nrow(c_iz),length(s_i),nrow(t_iz),length(v_i),length(PredTF_i))!=n_i) ) stop("b_i, a_i, c_i, s_i, v_i, or t_i doesn't have length n_i")
    if( nrow(a_xl)!=n_x | ncol(a_xl)!=n_l ) stop("a_xl has wrong dimensions")
    if( nrow(X_xj)!=n_x | ncol(X_xj)!=n_j ) stop("X_xj has wrong dimensions")
    if( nrow(Q_ik)!=n_i | ncol(Q_ik)!=n_k ) stop("Q_ik has wrong dimensions")
    if( dim(X_xtp)[1]!=n_x | dim(X_xtp)[2]!=n_t | dim(X_xtp)[3]!=n_p ) stop("X_xtp has wrong dimensions")
    if( ncol(c_iz)>1 & any(ObsModel_ez[,2]!=1) ) stop("Using multiple columnns in `c_iz` only makes sense using a Poisson-link delta model via `ObsModel[2]=1`")
  }

  # Check for incompatibilities amongst versions
  if( any(ObsModel_ez[,2]==1) ){
    # Versions 1.6.0 through 2.2.0 used a previous interpretation of area-swept for Poisson-link model and are not consistent with Q-Q plotting
    if( Version%in%c("VAST_v2_2_0","VAST_v2_1_0","VAST_v2_0_0","VAST_v1_9_0","VAST_v1_8_0","VAST_v1_7_0","VAST_v1_6_0") ){
      stop("Problem with Poisson-link model for VAST versions 1.6.0 through 2.2.0")
    }
    # Versions 1.0.0 through 1.5.0 don't have Poisson-link model
    if( Version%in%c("VAST_v1_5_0","VAST_v1_4_0","VAST_v1_3_0","VAST_v1_2_0","VAST_v1_1_0","VAST_v1_0_0") ){
      stop("Poisson-link model for VAST versions 1.0.0 through 1.5.0")
    }
    # Can't use multiple error distributions prior to version 3.0.0
    if( Version%in%c("VAST_v2_8_0","VAST_v2_7_0","VAST_v2_6_0","VAST_v2_5_0","VAST_v2_4_0","VAST_v2_3_0","VAST_v2_2_0","VAST_v2_1_0","VAST_v2_0_0","VAST_v1_9_0","VAST_v1_8_0","VAST_v1_7_0","VAST_v1_6_0","VAST_v1_5_0","VAST_v1_4_0","VAST_v1_3_0","VAST_v1_2_0","VAST_v1_1_0","VAST_v1_0_0") ){
      if( length(unique(ObsModel_ez[,1]))>1 ) stop("Can't use multiple error distributions prior to version 3.0.0")
    }
  }
  if( any(FieldConfig_input == -2) ){
    # Versions 2.6.0 was the first to allow "IID" setting for FieldConfig elements
    if( Version%in%c("VAST_v2_5_0","VAST_v2_4_0","VAST_v2_3_0","VAST_v2_2_0","VAST_v2_1_0","VAST_v2_0_0","VAST_v1_9_0","VAST_v1_8_0","VAST_v1_7_0","VAST_v1_6_0","VAST_v1_5_0","VAST_v1_4_0","VAST_v1_3_0","VAST_v1_2_0","VAST_v1_1_0","VAST_v1_0_0") ){
      stop("Problem with 'IID' option for factors for VAST versions 1.0.0 through 2.6.0")
    }
  }
  if( ncol(c_iz)>1 ){
    if( Version%in%c("VAST_v3_0_0","VAST_v2_8_0","VAST_v2_7_0","VAST_v2_6_0","VAST_v2_5_0","VAST_v2_4_0","VAST_v2_3_0","VAST_v2_2_0","VAST_v2_1_0","VAST_v2_0_0","VAST_v1_9_0","VAST_v1_8_0","VAST_v1_7_0","VAST_v1_6_0","VAST_v1_5_0","VAST_v1_4_0","VAST_v1_3_0","VAST_v1_2_0","VAST_v1_1_0","VAST_v1_0_0") ){
      stop("Can't have observations assigned to more than one category prior to version 4.0.0")
    }
  }
  if( any(RhoConfig[1:2]==3) ){
    if( Version%in%c("VAST_v4_1_0","VAST_v4_0_0","VAST_v3_0_0","VAST_v2_8_0","VAST_v2_7_0","VAST_v2_6_0","VAST_v2_5_0","VAST_v2_4_0","VAST_v2_3_0","VAST_v2_2_0","VAST_v2_1_0","VAST_v2_0_0","VAST_v1_9_0","VAST_v1_8_0","VAST_v1_7_0","VAST_v1_6_0","VAST_v1_5_0","VAST_v1_4_0","VAST_v1_3_0","VAST_v1_2_0","VAST_v1_1_0","VAST_v1_0_0") ){
      stop("There was bug in specifying fixed intercepts among years for versions prior to V4.2.0")
    }
  }
  if( Options2use['SD_observation_density']==1 ){
    if( Version%in%c("VAST_v4_1_0","VAST_v4_0_0","VAST_v3_0_0","VAST_v2_8_0","VAST_v2_7_0","VAST_v2_6_0","VAST_v2_5_0","VAST_v2_4_0","VAST_v2_3_0","VAST_v2_2_0","VAST_v2_1_0","VAST_v2_0_0","VAST_v1_9_0","VAST_v1_8_0","VAST_v1_7_0","VAST_v1_6_0","VAST_v1_5_0","VAST_v1_4_0","VAST_v1_3_0","VAST_v1_2_0","VAST_v1_1_0","VAST_v1_0_0") ){
      stop("Calculating 'SD_observation_density' is not possible prior to V4.2.0")
    }
  }


  ### Check for incompatible settings
  # Seasonal models and intercepts
  if( ncol(t_iz)>=2 & any(RhoConfig[1:2]!=0) ){
    stop("Temporal structure on intercepts is not implemented for seasonal models")
  }
  if( (FieldConfig_input[2]==(-1) & RhoConfig[3]!=0) | (FieldConfig_input[4]==(-1) & RhoConfig[4]!=0) ){
    stop("Spatio-temporal variation is turned off for a component with temporal structure, and this combination doesn't make sense")
  }
  # Interactions
  if( VamConfig[1]!=0 ){
    if( any(RhoConfig[1:2]!=3) | any(ObsModel_ez[,2]!=1) ){
      stop("Bad settings when turning on `VamConfig`")
    }
    if( !all(FieldConfig_input[c(2,4)] %in% c(-2,-1,0,n_c)) ){
      stop("Spatio-temporal variation must either have full rank covariance or be turned off or IID for each component for when using interactions")
    }
    if( VamConfig[2]>n_c | VamConfig[2]<0 ){
      stop("Rank for interactions must be an integer between 0 and `n_c`, inclusive")
    }
    if( VamConfig[2]==n_c & any(RhoConfig[3:4]!=0) ){
      stop("Can't simultaneously identify full-rank interactions and temporal correlation on spatio-temporal components")
    }
  }

  # switch defaults if necessary
  if( Method=="Grid" ){
    Aniso = 0
    message("Using isotropic 2D AR1 hyperdistribution, so switching to Aniso=0")
  }
  if( Method=="Spherical_mesh" ){
    Aniso = 0
    message("Using spherical projection for SPDE approximation, so switching to Aniso=0")
  }

  # Warning messages
  if( n_c>1 & any(FieldConfig_input==1)){
    warning( "Using 1 factor for more than one category:  Please note that this is non-standard, and it is more common to use multiple factors (often as many as the number of categories)" )
  }

  # Output tagged list
  # CMP_xmax should be >100 and CMP_breakpoint should be 1 for Tweedie model
  Options_vec = c("Aniso"=Aniso, "R2_interpretation"=0, "Rho_beta1_TF"=ifelse(RhoConfig[1]%in%c(1,2,4),1,0), "Rho_beta2_TF"=ifelse(RhoConfig[2]%in%c(1,2,4),1,0), "AreaAbundanceCurveTF"=0, "CMP_xmax"=200, "CMP_breakpoint"=1, "Method"=switch(Method,"Mesh"=0,"Grid"=1,"Spherical_mesh"=0) )
  if(Version%in%c("VAST_v1_1_0","VAST_v1_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=t_iz[,1]-min(t_iz[,1]), "a_xl"=a_xl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
  }
  if(Version%in%c("VAST_v1_4_0","VAST_v1_3_0","VAST_v1_2_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_f_input"=OverdispersionConfig_input[1], "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=t_iz[,1]-min(t_iz[,1]), "v_i"=match(v_i,sort(unique(v_i)))-1, "a_xl"=a_xl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
  }
  if(Version%in%c("VAST_v1_6_0","VAST_v1_5_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_f_input"=OverdispersionConfig_input[1], "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=t_iz[,1]-min(t_iz[,1]), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_xl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
  }
  if(Version%in%c("VAST_v1_7_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig_input, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=t_iz[,1]-min(t_iz[,1]), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_xl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
  }
  if(Version%in%c("VAST_v1_8_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig_input, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=t_iz[,1]-min(t_iz[,1]), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_xl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v1_9_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig_input, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_i"=t_iz[,1]-min(t_iz[,1]), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_xl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v2_8_0","VAST_v2_7_0","VAST_v2_6_0","VAST_v2_5_0","VAST_v2_4_0","VAST_v2_3_0","VAST_v2_2_0","VAST_v2_1_0","VAST_v2_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig_input, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel"=ObsModel_ez[1,], "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "s_i"=s_i, "t_iz"=t_iz-min(t_iz[,1],na.rm=TRUE), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_xl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v3_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig_input, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_i"=c_iz[,1], "e_i"=e_i, "s_i"=s_i, "t_iz"=t_iz-min(t_iz,na.rm=TRUE), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_xl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v4_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig_input, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=t_iz-min(t_iz,na.rm=TRUE), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_xl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v4_4_0","VAST_v4_3_0","VAST_v4_2_0","VAST_v4_1_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=n_p, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig_input, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "include_data"=TRUE, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=t_iz-min(t_iz,na.rm=TRUE), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_xl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if(Version%in%c("VAST_v5_0_0")){
    Return = list( "n_i"=n_i, "n_s"=c(MeshList$anisotropic_spde$n.spde,n_x)[Options_vec['Method']+1], "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_e"=n_e, "n_p"=n_p, "n_v"=n_v, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig_input, "OverdispersionConfig"=OverdispersionConfig_input, "ObsModel_ez"=ObsModel_ez, "VamConfig"=VamConfig, "include_data"=TRUE, "Options"=Options2use, "yearbounds_zz"=yearbounds_zz, "b_i"=b_i, "a_i"=a_i, "c_iz"=c_iz, "e_i"=e_i, "s_i"=s_i, "t_iz"=t_iz-min(t_iz,na.rm=TRUE), "v_i"=match(v_i,sort(unique(v_i)))-1, "PredTF_i"=PredTF_i, "a_xl"=a_xl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "t_yz"=t_yz, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list(), "M0"=GridList$M0, "M1"=GridList$M1, "M2"=GridList$M2 )
  }
  if( "spde" %in% names(Return) ) Return[['spde']] = MeshList$isotropic_spde$param.inla[c("M0","M1","M2")]
  if( "spde_aniso" %in% names(Return) ) Return[['spde_aniso']] = list("n_s"=MeshList$anisotropic_spde$n.spde, "n_tri"=nrow(MeshList$anisotropic_mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$anisotropic_spde$param.inla$M0, "G0_inv"=INLA::inla.as.dgTMatrix(solve(MeshList$anisotropic_spde$param.inla$M0)) )

  # Check for NAs
  if( CheckForErrors==TRUE ){
    NoNAs = setdiff( names(Return), c("t_iz","t_yz","c_iz") )
    if( any(sapply(Return[NoNAs], FUN=function(Array){any(is.na(Array))})==TRUE) ) stop("Please find and eliminate the NA from your inputs")
  }

  # Return
  return( Return )
}
