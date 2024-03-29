#' Calculate parameter inputs for TMB
#'
#' \code{make_parameters} generates the \code{parameters} input for \code{TMB::MakeADFun}
#'
#' @param DataList list outputted from \code{make_data}
#' @inheritParams make_data

#' @return Tagged list containing starting values for all fixed effects (parameters) and random effects (coefficients)
#' \describe{
#'   \item{ln_H_input}{two parameters governing geometric anisotropy (rotation matrix H)}
#'   \item{beta1_ct}{intercepts for encounter probability}
#'   \item{gamma1_cp}{effect of density covariates on encounter prob}
#'   \item{lambda1_k}{effect of catchability covariates on encounter prob}
#'   \item{L1_z}{trimmed Cholesky (i.e., parameters for square-root of covariance) for overdispersion in encounter prob}
#'   \item{L_omega1_z}{trimmed Cholesky pointwise variance in spatial variation in encounter prob}
#'   \item{L_epsilon1_z}{trimmed Cholesky pointwise variance in spatio-temporal variation in encounter prob}
#'   \item{logkappa1}{governs decorrelation distance in encounter prob}
#'   \item{Beta_mean1}{average intercept of encounter prob (used with RhoConfig options)}
#'   \item{logsigmaB1}{SD for intercept of encounter prob (used with RhoConfig options)}
#'   \item{Beta_rho1}{first-order autoregressive coefficient for intercept of encounter prob (used with RhoConfig options)}
#'   \item{Epsilon_rho1}{first-order autoregressive coefficient for spatio-temporal variation of encounter prob (used with RhoConfig options)}
#'   \item{eta1_vf}{overdispersion parameters (i.e., vessel or tow-level effects) on encounter prob}
#'   \item{Omegainput1_sf}{Spatial variation in encounter prob}
#'   \item{Epsiloninput1_sft}{Spatio-temporal variation in encounter prob}
#'   \item{beta2_ct}{intercepts for positive catch rates}
#'   \item{gamma2_cp}{effect of density covariates on positive catch rates}
#'   \item{lambda2_k}{effect of catchability covariates on positive catch rates}
#'   \item{L2_z}{trimmed Cholesky (i.e., parameters for square-root of covariance) for overdispersion in positive catch rates}
#'   \item{L_omega2_z}{trimmed Cholesky pointwise variance in spatial variation in positive catch rates}
#'   \item{L_epsilon2_z}{trimmed Cholesk pointwise variance in spatio-temporal variation in positive catch rates}
#'   \item{logkappa2}{governs decorrelation distance in positive catch rates}
#'   \item{Beta_mean2}{average intercept of positive catch rates (used with RhoConfig options)}
#'   \item{logsigmaB2}{SD for intercept of positive catch rates (used with RhoConfig options)}
#'   \item{Beta_rho2}{first-order autoregressive coefficient for intercept of positive catch rates (used with RhoConfig options)}
#'   \item{Epsilon_rho2}{first-order autoregressive coefficient for spatio-temporal variation of positive catch rates (used with RhoConfig options)}
#'   \item{eta2_vf}{overdispersion parameters (i.e., vessel or tow-level effects) on positive catch rates}
#'   \item{Omegainput2_sf}{Spatial variation in positive catch rates}
#'   \item{Epsiloninput2_sft}{Spatio-temporal variation in positive catch rates}
#'   \item{logSigmaM}{variance parameters for positive catch rates}
#' }

#' @export
make_parameters <-
function( Version,
          DataList,
          RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0) ){

  # Local function to make a random array
  rarray = function( dim, mean=0, sd=0.01 ) array( rnorm(prod(dim),mean=mean,sd=sd), dim=dim)

  # Local funciton to add L_z and eta_vf to Parameters
  Add_factor = function( List, n_c, n_f, n_i, n_t=NA, list_names, sd=0.01 ){
    # n_f factors
    if( n_f>0 ){
      List[[which(names(List)==list_names[1])]] = rnorm(n_f*n_c - n_f*(n_f-1)/2)
      List[[which(names(List)==list_names[2])]] = rarray(dim=as.vector(na.omit(c(n_i,n_f,n_t))), sd=sd)
    }
    # AR1 process
    if( n_f==0 ){
      List[[which(names(List)==list_names[1])]] = c(1,0.5) # Pointwise SD / Correlation
      List[[which(names(List)==list_names[2])]] = rarray(dim=as.vector(na.omit(c(n_i,n_c,n_t))), sd=sd)
    }
    # Turned off
    if( n_f== -1 ){
      if( FishStatsUtils::convert_version_name(Version) >= FishStatsUtils::convert_version_name("VAST_v12_0_0") ){
        List[[which(names(List)==list_names[1])]] = vector()  # Turn off SD when zero factors, i.e., n_f = -1, MUST BE 1.0 BY DEFAULT
        List[[which(names(List)==list_names[2])]] = rarray(dim=as.vector(na.omit(c(n_i,0,n_t))), sd=sd)
      }else{
        List[[which(names(List)==list_names[1])]] = 1  # Turn off SD when zero factors, i.e., n_f = -1, MUST BE 1.0 BY DEFAULT
        List[[which(names(List)==list_names[2])]] = rarray(dim=as.vector(na.omit(c(n_i,abs(n_f),n_t))), sd=sd)
      }
    }
    # IID process
    if( n_f== -2 ){
      List[[which(names(List)==list_names[1])]] = rep(1,n_c) # Pointwise SD / Correlation
      List[[which(names(List)==list_names[2])]] = rarray(dim=as.vector(na.omit(c(n_i,n_c,n_t))), sd=sd)
    }
    # Identity matrix (1.0 on diagonal; 0.0 off diagonal)
    if( n_f== -3 ){
      List[[which(names(List)==list_names[1])]] = vector() # Pointwise SD / Correlation
      List[[which(names(List)==list_names[2])]] = rarray(dim=as.vector(na.omit(c(n_i,n_c,n_t))), sd=sd)
    }
    return( List )
  }

  #### Deals with backwards compatibility for FieldConfig
  # Converts from 4-vector to 3-by-2 matrix
  if( is.vector(DataList$FieldConfig) && length(DataList$FieldConfig)==4 ){
    DataList$FieldConfig = rbind( matrix(DataList$FieldConfig,ncol=2,dimnames=list(c("Omega","Epsilon"),c("Component_1","Component_2"))), "Beta"=c("IID","IID") )
  }
  # Converts from 3-by-2 matrix to 4-by-2 matrix
  if( is.matrix(DataList$FieldConfig) & all(dim(DataList$FieldConfig)==c(3,2)) ){
    DataList$FieldConfig = rbind( DataList$FieldConfig, "Epsilon_time"=c("IID","IID") )
  }
  # Checks for errors
  if( !is.matrix(DataList$FieldConfig) || !all(dim(DataList$FieldConfig)==c(4,2)) ){
    stop("`FieldConfig` has the wrong dimensions in `make_data`")
  }
  # Renames
  dimnames(DataList$FieldConfig) = list( c("Omega","Epsilon","Beta","Epsilon_time"), c("Component_1","Component_2") )

  #######################
  # Make Parameters for each version
  #######################

  if(Version%in%c("VAST_v8_4_0","VAST_v8_3_0","VAST_v8_2_0","VAST_v8_1_0","VAST_v8_0_0")){
    Return = list("ln_H_input"=c(0,0), "Chi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "Psi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "beta1_ft"=NA, "gamma1_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda1_k"=rep(0,ncol(DataList$Q_ik)), "L1_z"=NA, "L_omega1_z"=NA, "L_epsilon1_z"=NA, "L_beta1_z"=NA, "logkappa1"=log(0.9), "Beta_mean1_c"=rep(0,DataList$n_c), "Beta_rho1_f"=NA, "Epsilon_rho1_f"=NA, "log_sigmaXi1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p)), "log_sigmaratio1_z"=0, "eta1_vf"=NA, "Xiinput1_scp"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p)), "Omegainput1_sf"=NA, "Epsiloninput1_sft"=NA, "beta2_ft"=NA, "gamma2_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda2_k"=rep(0,ncol(DataList$Q_ik)), "L2_z"=NA, "L_omega2_z"=NA, "L_epsilon2_z"=NA, "L_beta2_z"=NA, "logkappa2"=log(0.9), "Beta_mean2_c"=rep(0,DataList$n_c), "Beta_rho2_f"=NA, "Epsilon_rho2_f"=NA, "log_sigmaXi2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p)), "log_sigmaratio2_z"=0, "logSigmaM"=rep(1,DataList$n_e)%o%c(log(5),log(2),log(1)), "delta_i"=rnorm(n=ifelse(any(DataList$ObsModel_ez[,1]%in%c(11,14)),DataList$n_i,0),sd=0.1), "eta2_vf"=NA, "Xiinput2_scp"=rarray(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p)), "Omegainput2_sf"=NA, "Epsiloninput2_sft"=NA )
  }
  if(Version%in%c("VAST_v8_6_0","VAST_v8_5_0")){
    Return = list("ln_H_input"=c(0,0), "Chi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "Psi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "beta1_ft"=NA, "gamma1_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda1_k"=rep(0,ncol(DataList$Q_ik)), "L1_z"=NA, "L_omega1_z"=NA, "L_epsilon1_z"=NA, "L_beta1_z"=NA, "logkappa1"=log(0.9), "Beta_mean1_c"=rep(0,DataList$n_c), "Beta_mean1_t"=rep(0,DataList$n_t), "Beta_rho1_f"=NA, "Epsilon_rho1_f"=NA, "log_sigmaXi1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p)), "log_sigmaratio1_z"=0, "eta1_vf"=NA, "Xiinput1_scp"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p)), "Omegainput1_sf"=NA, "Epsiloninput1_sft"=NA, "beta2_ft"=NA, "gamma2_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda2_k"=rep(0,ncol(DataList$Q_ik)), "L2_z"=NA, "L_omega2_z"=NA, "L_epsilon2_z"=NA, "L_beta2_z"=NA, "logkappa2"=log(0.9), "Beta_mean2_c"=rep(0,DataList$n_c), "Beta_mean2_t"=rep(0,DataList$n_t), "Beta_rho2_f"=NA, "Epsilon_rho2_f"=NA, "log_sigmaXi2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p)), "log_sigmaratio2_z"=0, "logSigmaM"=rep(1,DataList$n_e)%o%c(log(5),log(2),log(1)), "delta_i"=rnorm(n=ifelse(any(DataList$ObsModel_ez[,1]%in%c(11,14)),DataList$n_i,0),sd=0.1), "eta2_vf"=NA, "Xiinput2_scp"=rarray(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p)), "Omegainput2_sf"=NA, "Epsiloninput2_sft"=NA )
  }
  if(Version%in%c("VAST_v9_0_0")){
    Return = list("ln_H_input"=c(0,0), "Chi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "Psi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "beta1_ft"=NA, "gamma1_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda1_k"=rep(0,ncol(DataList$Q_ik)), "L1_z"=NA, "L_omega1_z"=NA, "L_epsilon1_z"=NA, "L_beta1_z"=NA, "Ltime_epsilon1_z"=NA, "logkappa1"=log(0.9), "Beta_mean1_c"=rep(0,DataList$n_c), "Beta_mean1_t"=rep(0,DataList$n_t), "Beta_rho1_f"=NA, "Epsilon_rho1_f"=NA, "log_sigmaXi1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p)), "log_sigmaratio1_z"=0, "eta1_vf"=NA, "Xiinput1_scp"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p)), "Omegainput1_sf"=NA, "Epsiloninput1_sff"=NA, "beta2_ft"=NA, "gamma2_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda2_k"=rep(0,ncol(DataList$Q_ik)), "L2_z"=NA, "L_omega2_z"=NA, "L_epsilon2_z"=NA, "L_beta2_z"=NA, "Ltime_epsilon2_z"=NA, "logkappa2"=log(0.9), "Beta_mean2_c"=rep(0,DataList$n_c), "Beta_mean2_t"=rep(0,DataList$n_t), "Beta_rho2_f"=NA, "Epsilon_rho2_f"=NA, "log_sigmaXi2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p)), "log_sigmaratio2_z"=0, "logSigmaM"=rep(1,DataList$n_e)%o%c(log(5),log(2),log(1)), "delta_i"=rnorm(n=ifelse(any(DataList$ObsModel_ez[,1]%in%c(11,14)),DataList$n_i,0),sd=0.1), "eta2_vf"=NA, "Xiinput2_scp"=rarray(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p)), "Omegainput2_sf"=NA, "Epsiloninput2_sff"=NA )
  }
  if(Version%in%c("VAST_v9_3_0","VAST_v9_2_0","VAST_v9_1_0")){
    Return = list("ln_H_input"=c(0,0), "Chi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "Psi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "beta1_ft"=NA, "gamma1_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda1_k"=rep(0,ncol(DataList$Q_ik)), "L_eta1_z"=NA, "L_omega1_z"=NA, "L_epsilon1_z"=NA, "L_beta1_z"=NA, "Ltime_epsilon1_z"=NA, "logkappa1"=log(0.9), "Beta_mean1_c"=rep(0,DataList$n_c), "Beta_mean1_t"=rep(0,DataList$n_t), "Beta_rho1_f"=NA, "Epsilon_rho1_f"=NA, "log_sigmaXi1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p)), "log_sigmaratio1_z"=0, "eta1_vf"=NA, "Xiinput1_scp"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p)), "Omegainput1_sf"=NA, "Epsiloninput1_sff"=NA, "beta2_ft"=NA, "gamma2_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda2_k"=rep(0,ncol(DataList$Q_ik)), "L_eta2_z"=NA, "L_omega2_z"=NA, "L_epsilon2_z"=NA, "L_beta2_z"=NA, "Ltime_epsilon2_z"=NA, "logkappa2"=log(0.9), "Beta_mean2_c"=rep(0,DataList$n_c), "Beta_mean2_t"=rep(0,DataList$n_t), "Beta_rho2_f"=NA, "Epsilon_rho2_f"=NA, "log_sigmaXi2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p)), "log_sigmaratio2_z"=0, "logSigmaM"=rep(1,DataList$n_e)%o%c(log(5),log(2),log(1)), "delta_i"=rnorm(n=ifelse(any(DataList$ObsModel_ez[,1]%in%c(11,14)),DataList$n_i,0),sd=0.1), "eta2_vf"=NA, "Xiinput2_scp"=rarray(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p)), "Omegainput2_sf"=NA, "Epsiloninput2_sff"=NA )
  }
  if(Version%in%c("VAST_v9_4_0")){
    Return = list("ln_H_input"=c(0,0), "Chi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "Psi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "beta1_ft"=NA, "gamma1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p)), "lambda1_k"=rep(0,ncol(DataList$Q_ik)), "L_eta1_z"=NA, "L_omega1_z"=NA, "L_epsilon1_z"=NA, "L_beta1_z"=NA, "Ltime_epsilon1_z"=NA, "logkappa1"=log(0.9), "Beta_mean1_c"=rep(0,DataList$n_c), "Beta_mean1_t"=rep(0,DataList$n_t), "Beta_rho1_f"=NA, "Epsilon_rho1_f"=NA, "log_sigmaXi1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p)), "eta1_vf"=NA, "Xiinput1_scp"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p)), "Omegainput1_sf"=NA, "Epsiloninput1_sff"=NA, "beta2_ft"=NA, "gamma2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p)), "lambda2_k"=rep(0,ncol(DataList$Q_ik)), "L_eta2_z"=NA, "L_omega2_z"=NA, "L_epsilon2_z"=NA, "L_beta2_z"=NA, "Ltime_epsilon2_z"=NA, "logkappa2"=log(0.9), "Beta_mean2_c"=rep(0,DataList$n_c), "Beta_mean2_t"=rep(0,DataList$n_t), "Beta_rho2_f"=NA, "Epsilon_rho2_f"=NA, "log_sigmaXi2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p)), "logSigmaM"=rep(1,DataList$n_e)%o%c(log(5),log(2),log(1)), "delta_i"=rnorm(n=ifelse(any(DataList$ObsModel_ez[,1]%in%c(11,14)),DataList$n_i,0),sd=0.1), "eta2_vf"=NA, "Xiinput2_scp"=rarray(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p)), "Omegainput2_sf"=NA, "Epsiloninput2_sff"=NA )
  }
  if(Version%in%c("VAST_v10_0_0")){
    Return = list("ln_H_input"=c(0,0), "Chi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "Psi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "beta1_ft"=NA, "gamma1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p1)), "lambda1_k"=rep(0,ncol(DataList$Q_ik)), "L_eta1_z"=NA, "L_omega1_z"=NA, "L_epsilon1_z"=NA, "L_beta1_z"=NA, "Ltime_epsilon1_z"=NA, "logkappa1"=log(0.9), "Beta_mean1_c"=rep(0,DataList$n_c), "Beta_mean1_t"=rep(0,DataList$n_t), "Beta_rho1_f"=NA, "Epsilon_rho1_f"=NA, "log_sigmaXi1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p1)), "eta1_vf"=NA, "Xiinput1_scp"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p1)), "Omegainput1_sf"=NA, "Epsiloninput1_sff"=NA, "beta2_ft"=NA, "gamma2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p2)), "lambda2_k"=rep(0,ncol(DataList$Q_ik)), "L_eta2_z"=NA, "L_omega2_z"=NA, "L_epsilon2_z"=NA, "L_beta2_z"=NA, "Ltime_epsilon2_z"=NA, "logkappa2"=log(0.9), "Beta_mean2_c"=rep(0,DataList$n_c), "Beta_mean2_t"=rep(0,DataList$n_t), "Beta_rho2_f"=NA, "Epsilon_rho2_f"=NA, "log_sigmaXi2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p2)), "logSigmaM"=rep(1,DataList$n_e)%o%c(log(5),log(2),log(1)), "delta_i"=rnorm(n=ifelse(any(DataList$ObsModel_ez[,1]%in%c(11,14)),DataList$n_i,0),sd=0.1), "eta2_vf"=NA, "Xiinput2_scp"=rarray(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p2)), "Omegainput2_sf"=NA, "Epsiloninput2_sff"=NA )
  }
  if(Version%in%c("VAST_v11_0_0")){
    Return = list("ln_H_input"=c(0,0), "Chi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "Psi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "beta1_ft"=NA, "gamma1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p1)), "lambda1_k"=rep(0,ncol(DataList$Q1_ik)), "L_eta1_z"=NA, "L_omega1_z"=NA, "L_epsilon1_z"=NA, "L_beta1_z"=NA, "Ltime_epsilon1_z"=NA, "logkappa1"=log(0.9), "Beta_mean1_c"=rep(0,DataList$n_c), "Beta_mean1_t"=rep(0,DataList$n_t), "Beta_rho1_f"=NA, "Epsilon_rho1_f"=NA, "log_sigmaXi1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p1)), "eta1_vf"=NA, "Xiinput1_scp"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p1)), "Omegainput1_sf"=NA, "Epsiloninput1_sff"=NA, "beta2_ft"=NA, "gamma2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p2)), "lambda2_k"=rep(0,ncol(DataList$Q2_ik)), "L_eta2_z"=NA, "L_omega2_z"=NA, "L_epsilon2_z"=NA, "L_beta2_z"=NA, "Ltime_epsilon2_z"=NA, "logkappa2"=log(0.9), "Beta_mean2_c"=rep(0,DataList$n_c), "Beta_mean2_t"=rep(0,DataList$n_t), "Beta_rho2_f"=NA, "Epsilon_rho2_f"=NA, "log_sigmaXi2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p2)), "logSigmaM"=rep(1,DataList$n_e)%o%c(log(5),log(2),log(1)), "delta_i"=rnorm(n=ifelse(any(DataList$ObsModel_ez[,1]%in%c(11,14)),DataList$n_i,0),sd=0.1), "eta2_vf"=NA, "Xiinput2_scp"=rarray(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p2)), "Omegainput2_sf"=NA, "Epsiloninput2_sff"=NA )
  }
  if(Version%in%c("VAST_v13_0_0","VAST_v12_0_0")){
    Return = list("ln_H_input"=c(0,0), "Chi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "Psi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "beta1_ft"=NA, "gamma1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p1)), "lambda1_k"=rep(0,ncol(DataList$Q1_ik)), "L_eta1_z"=NA, "L_omega1_z"=NA, "L_epsilon1_z"=NA, "L_beta1_z"=NA, "Ltime_epsilon1_z"=NA, "logkappa1"=log(0.9), "Beta_mean1_c"=rep(0,DataList$n_c), "Beta_mean1_t"=rep(0,DataList$n_t), "Beta_rho1_f"=NA, "Epsilon_rho1_f"=NA, "log_sigmaXi1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p1)), "log_sigmaPhi1_k"=rep(0,ncol(DataList$Q1_ik)), "eta1_vf"=NA, "Xiinput1_scp"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p1)), "Phiinput1_sk"=array(0,dim=c(DataList$n_s,ncol(DataList$Q1_ik))), "Omegainput1_sf"=NA, "Epsiloninput1_sff"=NA, "beta2_ft"=NA, "gamma2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p2)), "lambda2_k"=rep(0,ncol(DataList$Q2_ik)), "L_eta2_z"=NA, "L_omega2_z"=NA, "L_epsilon2_z"=NA, "L_beta2_z"=NA, "Ltime_epsilon2_z"=NA, "logkappa2"=log(0.9), "Beta_mean2_c"=rep(0,DataList$n_c), "Beta_mean2_t"=rep(0,DataList$n_t), "Beta_rho2_f"=NA, "Epsilon_rho2_f"=NA, "log_sigmaXi2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p2)), "log_sigmaPhi2_k"=rep(0,ncol(DataList$Q2_ik)), "logSigmaM"=rep(1,DataList$n_e)%o%c(log(5),log(2),log(1)), "delta_i"=rnorm(n=ifelse(any(DataList$ObsModel_ez[,1]%in%c(11,14)),DataList$n_i,0),sd=0.1), "eta2_vf"=NA, "Xiinput2_scp"=rarray(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p2)), "Phiinput2_sk"=array(0,dim=c(DataList$n_s,ncol(DataList$Q2_ik))), "Omegainput2_sf"=NA, "Epsiloninput2_sff"=NA )
  }
  if(Version%in%c("VAST_v13_1_0")){
    Return = list("ln_H_input"=c(0,0), "Chi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "Psi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "beta1_ft"=NA, "gamma1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p1),dimnames=dimnames(DataList$X1_gctp)[c(2,4)]), "lambda1_k"=rep(0,ncol(DataList$Q1_ik)), "L_eta1_z"=NA, "L_omega1_z"=NA, "L_epsilon1_z"=NA, "L_beta1_z"=NA, "Ltime_epsilon1_z"=NA, "logkappa1"=log(0.9), "Beta_mean1_c"=rep(0,DataList$n_c), "Beta_mean1_t"=rep(0,DataList$n_t), "Beta_rho1_f"=NA, "Epsilon_rho1_f"=NA, "log_sigmaXi1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p1)), "log_sigmaPhi1_k"=rep(0,ncol(DataList$Q1_ik)), "eta1_vf"=NA, "Xiinput1_scp"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p1)), "Phiinput1_sk"=array(0,dim=c(DataList$n_s,ncol(DataList$Q1_ik))), "Omegainput1_sf"=NA, "Epsiloninput1_sff"=NA, "beta2_ft"=NA, "gamma2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p2),dimnames=dimnames(DataList$X2_gctp)[c(2,4)]), "lambda2_k"=rep(0,ncol(DataList$Q2_ik)), "L_eta2_z"=NA, "L_omega2_z"=NA, "L_epsilon2_z"=NA, "L_beta2_z"=NA, "Ltime_epsilon2_z"=NA, "logkappa2"=log(0.9), "Beta_mean2_c"=rep(0,DataList$n_c), "Beta_mean2_t"=rep(0,DataList$n_t), "Beta_rho2_f"=NA, "Epsilon_rho2_f"=NA, "log_sigmaXi2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p2)), "log_sigmaPhi2_k"=rep(0,ncol(DataList$Q2_ik)), "logSigmaM"=rep(1,DataList$n_e)%o%c(log(5),log(2),log(1)), "lagrange_tc"=rarray(0,dim=c(DataList$n_t,DataList$n_c)), "delta_i"=rnorm(n=ifelse(any(DataList$ObsModel_ez[,1]%in%c(11,14)),DataList$n_i,0),sd=0.1), "eta2_vf"=NA, "Xiinput2_scp"=rarray(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p2)), "Phiinput2_sk"=array(0,dim=c(DataList$n_s,ncol(DataList$Q2_ik))), "Omegainput2_sf"=NA, "Epsiloninput2_sff"=NA )
  }
  if(Version%in%c("VAST_v14_0_1","VAST_v14_0_0")){
    Return = list("ln_H_input"=c(0,0), "Chi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "Psi_fr"=rarray(dim=c(max(DataList$FieldConfig[2,1],1),DataList$VamConfig[2])), "beta1_ft"=NA, "gamma1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p1),dimnames=dimnames(DataList$X1_gctp)[c(2,4)]), "lambda1_k"=rep(0,ncol(DataList$Q1_ik)), "L_eta1_z"=NA, "L_omega1_z"=NA, "L_epsilon1_z"=NA, "L_beta1_z"=NA, "Ltime_epsilon1_z"=NA, "logkappa1"=log(0.9), "Beta_mean1_c"=rep(0,DataList$n_c), "Beta_mean1_t"=rep(0,DataList$n_t), "Beta_rho1_f"=NA, "Epsilon_rho1_f"=NA, "log_sigmaXi1_cp"=array(0,dim=c(DataList$n_c,DataList$n_p1)), "log_sigmaPhi1_k"=rep(0,ncol(DataList$Q1_ik)), "eta1_vf"=NA, "Xiinput1_scp"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p1)), "Phiinput1_sk"=array(0,dim=c(DataList$n_s,ncol(DataList$Q1_ik))), "Omegainput1_sf"=NA, "Epsiloninput1_sff"=NA, "beta2_ft"=NA, "gamma2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p2),dimnames=dimnames(DataList$X2_gctp)[c(2,4)]), "lambda2_k"=rep(0,ncol(DataList$Q2_ik)), "L_eta2_z"=NA, "L_omega2_z"=NA, "L_epsilon2_z"=NA, "L_beta2_z"=NA, "Ltime_epsilon2_z"=NA, "logkappa2"=log(0.9), "Beta_mean2_c"=rep(0,DataList$n_c), "Beta_mean2_t"=rep(0,DataList$n_t), "Beta_rho2_f"=NA, "Epsilon_rho2_f"=NA, "log_sigmaXi2_cp"=array(0,dim=c(DataList$n_c,DataList$n_p2)), "log_sigmaPhi2_k"=rep(0,ncol(DataList$Q2_ik)), "logSigmaM"=rep(1,DataList$n_e)%o%c(log(5),log(2),log(1)), "lagrange_tc"=rarray(0,dim=c(DataList$n_t,DataList$n_c)), "delta_i"=rnorm(n=ifelse(any(DataList$ObsModel_ez[,1]%in%c(11,14)),DataList$n_i,0),sd=0.1), "eta2_vf"=NA, "Xiinput2_scp"=rarray(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_p2)), "Phiinput2_sk"=array(0,dim=c(DataList$n_s,ncol(DataList$Q2_ik))), "Omegainput2_sf"=NA, "Epsiloninput2_sff"=NA, "eps_Index_ctl"=array(0,dim=c(0,0,0)) )
  }

  #######################
  # Fill in values that are shared across versions
  #######################

  # Overdispersion
  if( "n_f_input" %in% names(DataList) ){
    if( "L1_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$n_f_input, n_i=DataList$n_v, list_names=c("L1_z","eta1_vf"), sd=0 )
    if( "L2_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$n_f_input, n_i=DataList$n_v, list_names=c("L2_z","eta2_vf"), sd=0 )
    if( "L_eta1_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$n_f_input, n_i=DataList$n_v, list_names=c("L_eta1_z","eta1_vf"), sd=0 )
    if( "L_eta2_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$n_f_input, n_i=DataList$n_v, list_names=c("L_eta2_z","eta2_vf"), sd=0 )
  }
  if( "OverdispersionConfig" %in% names(DataList) ){
    if( "L1_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$OverdispersionConfig[1], n_i=DataList$n_v, list_names=c("L1_z","eta1_vf"), sd=0 )
    if( "L2_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$OverdispersionConfig[2], n_i=DataList$n_v, list_names=c("L2_z","eta2_vf"), sd=0 )
    if( "L_eta1_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$OverdispersionConfig[1], n_i=DataList$n_v, list_names=c("L_eta1_z","eta1_vf"), sd=0 )
    if( "L_eta2_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$OverdispersionConfig[2], n_i=DataList$n_v, list_names=c("L_eta2_z","eta2_vf"), sd=0 )
  }

  # Fields
  if( "L_omega1_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$FieldConfig[1,1], n_i=DataList$n_s, list_names=c("L_omega1_z","Omegainput1_sf"), sd=0 )
  if( "Epsiloninput1_sft" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$FieldConfig[2,1], n_i=DataList$n_s, n_t=DataList$n_t, list_names=c("L_epsilon1_z","Epsiloninput1_sft"), sd=0 )
  if( "Epsiloninput1_sff" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$FieldConfig[2,1], n_i=DataList$n_s, n_t=DataList$n_t, list_names=c("L_epsilon1_z","Epsiloninput1_sff"), sd=0 )
  if( "L_beta1_z" %in% names(Return)){
    Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$FieldConfig[3,1], n_i=DataList$n_t, list_names=c("L_beta1_z","beta1_ft"), sd=0 )
    Return[["beta1_ft"]] = t(Return[["beta1_ft"]])
  }
  if( "Ltime_epsilon1_z" %in% names(Return)){
    if(DataList$FieldConfig[4,1] <= 0) Return[["Ltime_epsilon1_z"]] = vector()
    if(DataList$FieldConfig[4,1] >= 1){
      Return[["Ltime_epsilon1_z"]] = rnorm(DataList$FieldConfig[4,1]*DataList$n_t - DataList$FieldConfig[4,1]*(DataList$FieldConfig[4,1]-1)/2)
      Return[["Epsiloninput1_sff"]] = Return[["Epsiloninput1_sff"]][,,1:DataList$FieldConfig[4,1],drop=FALSE]
    }
  }
  if( "L_omega2_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$FieldConfig[1,2], n_i=DataList$n_s, list_names=c("L_omega2_z","Omegainput2_sf"), sd=0 )
  if( "Epsiloninput2_sft" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$FieldConfig[2,2], n_i=DataList$n_s, n_t=DataList$n_t, list_names=c("L_epsilon2_z","Epsiloninput2_sft"), sd=0 )
  if( "Epsiloninput2_sff" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$FieldConfig[2,2], n_i=DataList$n_s, n_t=DataList$n_t, list_names=c("L_epsilon2_z","Epsiloninput2_sff"), sd=0 )
  if( "L_beta2_z" %in% names(Return)){
    Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$FieldConfig[3,2], n_i=DataList$n_t, list_names=c("L_beta2_z","beta2_ft"), sd=0 )
    Return[["beta2_ft"]] = t(Return[["beta2_ft"]])
  }
  if( "Ltime_epsilon2_z" %in% names(Return)){
    if(DataList$FieldConfig[4,2] <= 0) Return[["Ltime_epsilon2_z"]] = vector()
    if(DataList$FieldConfig[4,2] >= 1){
      Return[["Ltime_epsilon2_z"]] = rnorm(DataList$FieldConfig[4,2]*DataList$n_t - DataList$FieldConfig[4,2]*(DataList$FieldConfig[4,2]-1)/2)
      Return[["Epsiloninput2_sff"]] = Return[["Epsiloninput2_sff"]][,,1:DataList$FieldConfig[4,2],drop=FALSE]
    }
  }

  # Autocorrelation (must be ZERO by default)
  if( "Epsilon_rho1_f" %in% names(Return)){
    if("Epsiloninput1_sft" %in% names(Return)) Return[["Epsilon_rho1_f"]] = rep(0, dim(Return[["Epsiloninput1_sft"]])[2])
    if("Epsiloninput1_sff" %in% names(Return)) Return[["Epsilon_rho1_f"]] = rep(0, dim(Return[["Epsiloninput1_sff"]])[2])
  }
  if( "Epsilon_rho2_f" %in% names(Return)){
    if("Epsiloninput2_sft" %in% names(Return)) Return[["Epsilon_rho2_f"]] = rep(0, dim(Return[["Epsiloninput2_sft"]])[2])
    if("Epsiloninput2_sff" %in% names(Return)) Return[["Epsilon_rho2_f"]] = rep(0, dim(Return[["Epsiloninput2_sff"]])[2])
  }
  if( "Beta_rho1_f" %in% names(Return)) Return[["Beta_rho1_f"]] = rep(0, nrow(Return[["beta1_ft"]]))
  if( "Beta_rho2_f" %in% names(Return)) Return[["Beta_rho2_f"]] = rep(0, nrow(Return[["beta2_ft"]]))

  # Informative starting values < 7.0.0 OR >= 7.0.0 when not using factor-model feature
  Use_informative_starts = FALSE
  if( all(c("beta1_ct","beta2_ct") %in% names(Return)) ){
    Use_informative_starts = TRUE
  }
  if( all(c("beta1_ft","beta2_ft") %in% names(Return)) ){
    if( all(DataList$FieldConfig[3,1:2] == -2) ){
      Use_informative_starts = TRUE
    }
  }
  if( Use_informative_starts==TRUE ){
    # Temporary object for mapping
    Params_tmp = list( "beta1_ct"=NA, "beta2_ct"=NA )
    # Starting values
    if( all(DataList$ObsModel_ez[,2] %in% c(0,3)) ){
      Prop_ct = abind::adrop(DataList$Options_list$metadata_ctz[,,'prop_nonzero',drop=FALSE], drop=3)
      Params_tmp[["beta1_ct"]] = qlogis( 0.01 + 0.98*Prop_ct )
      Posdens_ct = abind::adrop(DataList$Options_list$metadata_ctz[,,'mean_posdens',drop=FALSE], drop=3)
      Params_tmp[["beta2_ct"]] = log(Posdens_ct)
    }
    if( all(DataList$ObsModel_ez[,2] %in% c(1,2,4)) ){
      Params_tmp[["beta1_ct"]] = array(0, dim=c(DataList$n_c,DataList$n_t))
      Params_tmp[["beta2_ct"]] = array(0, dim=c(DataList$n_c,DataList$n_t))
    }
    # Over-ride starting values for 100% and 0% encounters
    if( all(DataList$ObsModel_ez[,2] %in% c(3)) ){
      Prop_ct = abind::adrop(DataList$Options_list$metadata_ctz[,,'prop_nonzero',drop=FALSE], drop=3)
      if( any(is.na(Prop_ct) | Prop_ct==1) ) Params_tmp[["beta1_ct"]][which(is.na(Prop_ct) | Prop_ct==1)] = 20
    }
    if( all(DataList$ObsModel_ez[,2] %in% c(4)) ){
      Prop_ct = abind::adrop(DataList$Options_list$metadata_ctz[,,'prop_nonzero',drop=FALSE], drop=3)
      if( any(is.na(Prop_ct) | Prop_ct==1) ) Params_tmp[["beta1_ct"]][which(is.na(Prop_ct) | Prop_ct==1)] = 20
      if( any(is.na(Prop_ct) | Prop_ct==0) ) Params_tmp[["beta1_ct"]][which(is.na(Prop_ct) | Prop_ct==0)] = -20
      if( any(is.na(Prop_ct) | Prop_ct==0) ) Params_tmp[["beta2_ct"]][which(is.na(Prop_ct) | Prop_ct==0)] = 0
    }
    # Deal with any potential problems
    if( any(is.na(Params_tmp[["beta1_ct"]])) ){
      Params_tmp[["beta1_ct"]] = array(0, dim=c(DataList$n_c,DataList$n_t))
    }
    if( any(is.na(Params_tmp[["beta2_ct"]])) ){
      Params_tmp[["beta2_ct"]] = array(0, dim=c(DataList$n_c,DataList$n_t))
    }
    # Insert with name appropriate for a given version
    if( all(c("beta1_ct","beta2_ct") %in% names(Return)) ){
      Return[["beta1_ct"]] = Params_tmp[["beta1_ct"]]
      Return[["beta2_ct"]] = Params_tmp[["beta2_ct"]]
    }
    if( all(c("beta1_ft","beta2_ft") %in% names(Return)) ){
      Return[["beta1_ft"]] = Params_tmp[["beta1_ct"]]
      Return[["beta2_ft"]] = Params_tmp[["beta2_ct"]]
    }
  }

  # If either beta or epsilon is a random-walk process, fix starting value at 1
  if( RhoConfig[["Beta1"]]==2 ){
    if( "Beta_rho1" %in% names(Return) ) Return[["Beta_rho1"]] = 1
    if( "Beta_rho1_f" %in% names(Return) ) Return[["Beta_rho1_f"]][] = 1
  }
  if( RhoConfig[["Beta2"]]==2 ){
    if( "Beta_rho2" %in% names(Return) ) Return[["Beta_rho2"]] = 1
    if( "Beta_rho2_f" %in% names(Return) ) Return[["Beta_rho2_f"]][] = 1
  }
  if( RhoConfig[["Epsilon1"]] %in% c(2) ){
    if( "Epsilon_rho1"%in%names(Return) ) Return[["Epsilon_rho1"]] = 1
    if( "Epsilon_rho1_f"%in%names(Return) ) Return[["Epsilon_rho1_f"]][] = 1
  }
  if( RhoConfig[["Epsilon2"]] %in% c(2) ){
    if( "Epsilon_rho2"%in%names(Return) ) Return[["Epsilon_rho2"]] = 1
    if( "Epsilon_rho2_f"%in%names(Return) ) Return[["Epsilon_rho2_f"]][] = 1
  }

  # If either beta or epsilon is a AR1 process, fix starting value at 0.01 to ensure a non-zero starting gradient
  if( RhoConfig[["Beta1"]]==4 ){
    if( "Beta_rho1" %in% names(Return) ) Return[["Beta_rho1"]] = 0.01
    if( "Beta_rho1_f" %in% names(Return) ) Return[["Beta_rho1_f"]][] = 0.01
  }
  if( RhoConfig[["Beta2"]]==4 ){
    if( "Beta_rho2"%in%names(Return) ) Return[["Beta_rho2"]] = 0.01
    if( "Beta_rho2_f"%in%names(Return) ) Return[["Beta_rho2_f"]][] = 0.01
  }
  if( RhoConfig[["Epsilon1"]] %in% c(4,5) ){
    if( "Epsilon_rho1"%in%names(Return) ) Return[["Epsilon_rho1"]] = 0.01
    if( "Epsilon_rho1_f"%in%names(Return) ) Return[["Epsilon_rho1_f"]][] = 0.01
  }
  if( RhoConfig[["Epsilon2"]] %in% c(4,5) ){
    if( "Epsilon_rho2"%in%names(Return) ) Return[["Epsilon_rho2"]] = 0.01
    if( "Epsilon_rho2_f"%in%names(Return) ) Return[["Epsilon_rho2_f"]][] = 0.01
  }

  # If estimating habitat-covariates, start coefficient at nonzero value
  # THIS HAS BEEN CHANGED BACK TO STARTING AT 0.0 TO AVOID ISSUES WHEN PEOPLE MAP THESE OFF BY HAND
  #if( "Xconfig_zcp" %in% names(DataList) ){
  #  for(cI in 1:DataList$n_c){
  #  for(pI in 1:DataList$n_p){
  #    if( DataList$Xconfig_zcp[1,cI,pI] %in% c(1,3) ){
  #      Return[["gamma1_ctp"]][cI,,pI] = 0.1
  #    }
  #    if( DataList$Xconfig_zcp[2,cI,pI] %in% c(1,3) ){
  #      Return[["gamma2_ctp"]][cI,,pI] = 0.1
  #    }
  #  }}
  #}

  # replace missing values function
  tmpfn = function( vec ){
    Return = ifelse( abs(vec)==Inf, NA, vec)
    return( ifelse(is.na(Return),mean(Return,na.rm=TRUE),Return) )
  }
  if( all(c("beta1_ct","beta2_ct") %in% names(Return)) ){
    Return[["beta1_ct"]] = tmpfn( Return[["beta1_ct"]] )
    Return[["beta2_ct"]] = tmpfn( Return[["beta2_ct"]] )
  }
  if( all(c("beta1_ft","beta2_ft") %in% names(Return)) ){
    Return[["beta1_ft"]] = tmpfn( Return[["beta1_ft"]] )
    Return[["beta2_ft"]] = tmpfn( Return[["beta2_ft"]] )
  }

  # Interactions
  if( "VamConfig"%in%names(DataList) & all(c("Chi_fr","Psi_fr")%in%names(Return)) ){
    Return[["Psi_fr"]][1:ncol(Return[["Psi_fr"]]),] = diag(nrow=ncol(Return[["Psi_fr"]]))
    if( DataList$VamConfig[1]==2 ){
      if( "Epsilon_rho1" %in% names(Return) ) Mean = Return[["Epsilon_rho1"]]
      if( "Epsilon_rho1_f" %in% names(Return) ) Mean = mean(Return[["Epsilon_rho1_f"]])
      Return[["Psi_fr"]][cbind(1:ncol(Return[["Psi_fr"]]),1:ncol(Return[["Psi_fr"]]))] = seq(0.2,0.9,length=ncol(Return[["Psi_fr"]])) - Mean  # "diag" threw an error when ncol(Return[["Psi_fr"]])=1
    }
  }

  # Error messages
  if( any(sapply(Return, FUN=function(num){any(is.na(num))})) ){
    stop("Some parameter is NA")
  }

  # Return tagged list
  return( Return )
}
