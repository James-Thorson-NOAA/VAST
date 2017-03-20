#' Calculate parameter inputs for TMB
#'
#' \code{Param_Fn} generates the \code{parameters} input for \code{TMB::MakeADFun}
#'
#' @param DataList list outputted from \code{VAST::Data_Fn}
#' @inheritParams Data_Fn

#' @return Tagged list containing starting values for all fixed effects (parameters) and random effects (coefficients)
#' \describe{
#'   \item{ln_H_input}{two parameters governing geometric anisotropy (rotation matrix H)}
#'   \item{beta1_ct}{intercepts for encounter probability}
#'   \item{gamma1_j}{effect of density covariates that are static over time on encounter prob}
#'   \item{gamma1_ctp}{effect of density covariates that change over time on encounter prob}
#'   \item{lambda1_k}{effect of catchability covariates on encounter prob}
#'   \item{L1_z}{trimmed cholesky (i.e., parameters for square-root of covariance) for overdispersion in encounter prob}
#'   \item{L_omega1_z}{trimmed cholesky pointwise variance in spatial variation in encounter prob}
#'   \item{L_epsilon1_z}{trimmed cholesky pointwise variance in spatio-temporal variation in encounter prob}
#'   \item{logkappa1}{governs decorrelation distance in encounter prob}
#'   \item{Beta_mean1}{average intercept of encounter prob (used with RhoConfig options)}
#'   \item{logsigmaB1}{SD for intercept of encounter prob (used with RhoConfig options)}
#'   \item{Beta_rho1}{first-order autoregressive coefficient for intercept of encounter prob (used with RhoConfig options)}
#'   \item{Epsilon_rho1}{first-order autoregressive coefficient for spatio-temporal variation of encounter prob (used with RhoConfig options)}
#'   \item{eta1_vf}{overdispersion parameters (i.e., vessel or tow-level effects) on encounter prob}
#'   \item{Omegainput1_sf}{Spatial variation in encounter prob}
#'   \item{Epsiloninput1_sft}{Spatio-temporal variation in encounter prob}
#'   \item{beta2_ct}{intercepts for positive catch rates}
#'   \item{gamma2_j}{effect of density covariates that are static over time on positive catch rates}
#'   \item{gamma2_ctp}{effect of density covariates that change over time on positive catch rates}
#'   \item{lambda2_k}{effect of catchability covariates on positive catch rates}
#'   \item{L2_z}{trimmed cholesky (i.e., parameters for square-root of covariance) for overdispersion in positive catch rates}
#'   \item{L_omega2_z}{trimmed cholesky pointwise variance in spatial variation in positive catch rates}
#'   \item{L_epsilon2_z}{trimmed cholesky pointwise variance in spatio-temporal variation in positive catch rates}
#'   \item{logkappa2}{governs decorrelation distance in positive catch rates}
#'   \item{Beta_mean2}{average intercept of positive catch rates (used with RhoConfig options)}
#'   \item{logsigmaB2}{SD for intercept of positive catch rates (used with RhoConfig options)}
#'   \item{Beta_rho2}{first-order autoregressive coefficient for intercept of positive catch rates (used with RhoConfig options)}
#'   \item{Epsilon_rho2}{first-order autoregressive coefficient for spatio-temporal variation of positive catch rates (used with RhoConfig options)}
#'   \item{eta2_vf}{overdispersion parameters (i.e., vessel or tow-level effects) on positive catch rates}
#'   \item{Omegainput2_sf}{Spatial variation in positive catch rates}
#'   \item{Epsiloninput2_sft}{Spatio-temporal variation in positive catch rates}
#'   \item{logSigmaM}{variance parameters for positive catch ratesA}
#' }

#' @export
Param_Fn <-
function( Version, DataList, RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0) ){
  #
  rarray = function( dim, mean=0, sd=0.01 ) array( rnorm(prod(dim),mean=mean,sd=sd), dim=dim)
  # Add L_z and eta_vf to Parameters
  Add_factor = function( List, n_c, n_f, n_i, n_t=NA, list_names, sd=0.01 ){
    if( n_f>0 ){
      List[[which(names(List)==list_names[1])]] = rnorm(n_f*n_c - n_f*(n_f-1)/2)
      List[[which(names(List)==list_names[2])]] = rarray(dim=as.vector(na.omit(c(n_i,n_f,n_t))), sd=sd)
    }
    if( n_f<=0 ){
      List[[which(names(List)==list_names[1])]] = c(1,0.5) # Pointwise SD / Correlation
      List[[which(names(List)==list_names[2])]] = rarray(dim=as.vector(na.omit(c(n_i,n_c,n_t))), sd=sd)
    }
    return( List )
  }
  #
  if(Version%in%c("VAST_v1_0_0")){
    Return = list("ln_H_input"=c(0,0), "beta1_ct"=NA, "gamma1_j"=rep(0,DataList$n_j), "gamma1_tp"=matrix(0,nrow=DataList$n_t,ncol=DataList$n_p), "lambda1_k"=rep(0,DataList$n_k), "logetaE1"=0, "logetaO1"=0, "logkappa1"=0, "Beta_mean1"=0, "logsigmaB1"=log(1), "Beta_rho1"=0, "Epsilon_rho1"=0, "rho_c1"=0, "Omegainput1_sc"=matrix(0,DataList$n_s,DataList$n_c), "Epsiloninput1_sct"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_t)), "beta2_ct"=NA, "gamma2_j"=rep(0,DataList$n_j), "gamma2_tp"=matrix(0,nrow=DataList$n_t,ncol=DataList$n_p), "lambda2_k"=rep(0,DataList$n_k), "logetaE2"=0, "logetaO2"=0, "logkappa2"=0, "Beta_mean2"=0, "logsigmaB2"=log(1), "Beta_rho2"=0, "Epsilon_rho2"=0, "rho_c2"=0, "logSigmaM"=c(log(5),log(2)), "Omegainput2_sc"=matrix(0,DataList$n_s,DataList$n_c), "Epsiloninput2_sct"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_t)) )
  }
  if(Version%in%c("VAST_v1_1_0")){
    Return = list("ln_H_input"=c(0,0), "beta1_ct"=NA, "gamma1_j"=rep(0,DataList$n_j), "gamma1_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda1_k"=rep(0,DataList$n_k), "logetaE1"=0, "logetaO1"=0, "logkappa1"=0, "Beta_mean1"=0, "logsigmaB1"=log(1), "Beta_rho1"=0, "Epsilon_rho1"=0, "rho_c1"=0, "Omegainput1_sc"=matrix(0,DataList$n_s,DataList$n_c), "Epsiloninput1_sct"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_t)), "beta2_ct"=NA, "gamma2_j"=rep(0,DataList$n_j), "gamma2_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda2_k"=rep(0,DataList$n_k), "logetaE2"=0, "logetaO2"=0, "logkappa2"=0, "Beta_mean2"=0, "logsigmaB2"=log(1), "Beta_rho2"=0, "Epsilon_rho2"=0, "rho_c2"=0, "logSigmaM"=c(log(5),log(2)), "Omegainput2_sc"=matrix(0,DataList$n_s,DataList$n_c), "Epsiloninput2_sct"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_t)) )
  }
  if(Version%in%c("VAST_v1_3_0","VAST_v1_2_0")){
    Return = list("ln_H_input"=c(0,0), "beta1_ct"=NA, "lambda1_k"=rep(0,DataList$n_k), "L1_z"=NA, "logetaE1"=0, "logetaO1"=0, "logkappa1"=0, "Beta_mean1"=0, "logsigmaB1"=log(1), "Beta_rho1"=0, "Epsilon_rho1"=0, "rho_c1"=0, "eta1_vf"=NA, "Omegainput1_sc"=matrix(0,DataList$n_s,DataList$n_c), "Epsiloninput1_sct"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_t)), "beta2_ct"=NA, "gamma2_j"=rep(0,DataList$n_j), "gamma2_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda2_k"=rep(0,DataList$n_k), "L2_z"=NA, "logetaE2"=0, "logetaO2"=0, "logkappa2"=0, "Beta_mean2"=0, "logsigmaB2"=log(1), "Beta_rho2"=0, "Epsilon_rho2"=0, "rho_c2"=0, "logSigmaM"=rep(1,DataList$n_c)%o%c(log(5),log(2)), "eta2_vf"=NA, "Omegainput2_sc"=matrix(0,DataList$n_s,DataList$n_c), "Epsiloninput2_sct"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_t)) )
  }
  if(Version%in%c("VAST_v1_5_0","VAST_v1_4_0")){
    Return = list("ln_H_input"=c(0,0), "beta1_ct"=NA, "gamma1_j"=rep(0,DataList$n_j), "gamma1_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda1_k"=rep(0,DataList$n_k), "L1_z"=NA, "L_omega1_z"=NA, "L_epsilon1_z"=NA, "logkappa1"=0, "Beta_mean1"=0, "logsigmaB1"=log(1), "Beta_rho1"=0, "Epsilon_rho1"=0, "eta1_vf"=NA, "Omegainput1_sf"=NA, "Epsiloninput1_sft"=NA, "beta2_ct"=NA, "gamma2_j"=rep(0,DataList$n_j), "gamma2_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda2_k"=rep(0,DataList$n_k), "L2_z"=NA, "L_omega2_z"=NA, "L_epsilon2_z"=NA, "logkappa2"=0, "Beta_mean2"=0, "logsigmaB2"=log(1), "Beta_rho2"=0, "Epsilon_rho2"=0, "logSigmaM"=rep(1,DataList$n_c)%o%c(log(5),log(2)), "eta2_vf"=NA, "Omegainput2_sf"=NA, "Epsiloninput2_sft"=NA )
  }
  if(Version%in%c("VAST_v1_9_0","VAST_v1_8_0","VAST_v1_7_0","VAST_v1_6_0")){
    Return = list("ln_H_input"=c(0,0), "beta1_ct"=NA, "gamma1_j"=rep(0,DataList$n_j), "gamma1_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda1_k"=rep(0,DataList$n_k), "L1_z"=NA, "L_omega1_z"=NA, "L_epsilon1_z"=NA, "logkappa1"=log(0.9), "Beta_mean1"=0, "logsigmaB1"=log(1), "Beta_rho1"=0, "Epsilon_rho1"=0, "eta1_vf"=NA, "Omegainput1_sf"=NA, "Epsiloninput1_sft"=NA, "beta2_ct"=NA, "gamma2_j"=rep(0,DataList$n_j), "gamma2_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda2_k"=rep(0,DataList$n_k), "L2_z"=NA, "L_omega2_z"=NA, "L_epsilon2_z"=NA, "logkappa2"=log(0.9), "Beta_mean2"=0, "logsigmaB2"=log(1), "Beta_rho2"=0, "Epsilon_rho2"=0, "logSigmaM"=rep(1,DataList$n_c)%o%c(log(5),log(2),log(5)), "eta2_vf"=NA, "Omegainput2_sf"=NA, "Epsiloninput2_sft"=NA )
  }
  if(Version%in%c("VAST_v2_4_0","VAST_v2_3_0","VAST_v2_2_0","VAST_v2_1_0","VAST_v2_0_0")){
    Return = list("ln_H_input"=c(0,0), "beta1_ct"=NA, "gamma1_j"=rep(0,DataList$n_j), "gamma1_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda1_k"=rep(0,DataList$n_k), "L1_z"=NA, "L_omega1_z"=NA, "L_epsilon1_z"=NA, "logkappa1"=log(0.9), "Beta_mean1"=0, "logsigmaB1"=log(1), "Beta_rho1"=0, "Epsilon_rho1"=0, "log_sigmaratio1_z"=rep(0,ncol(DataList$t_iz)), "eta1_vf"=NA, "Omegainput1_sf"=NA, "Epsiloninput1_sft"=NA, "beta2_ct"=NA, "gamma2_j"=rep(0,DataList$n_j), "gamma2_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda2_k"=rep(0,DataList$n_k), "L2_z"=NA, "L_omega2_z"=NA, "L_epsilon2_z"=NA, "logkappa2"=log(0.9), "Beta_mean2"=0, "logsigmaB2"=log(1), "Beta_rho2"=0, "Epsilon_rho2"=0, "log_sigmaratio2_z"=rep(0,ncol(DataList$t_iz)), "logSigmaM"=rep(1,DataList$n_c)%o%c(log(5),log(2),log(5)), "eta2_vf"=NA, "Omegainput2_sf"=NA, "Epsiloninput2_sft"=NA )
  }
  # Overdispersion
  if( "n_f_input" %in% names(DataList) ){
    if( "L1_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$n_f_input, n_i=DataList$n_v, list_names=c("L1_z","eta1_vf"), sd=0 )
    if( "L2_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$n_f_input, n_i=DataList$n_v, list_names=c("L2_z","eta2_vf"), sd=0 )
  }
  if( "OverdispersionConfig" %in% names(DataList) ){
    if( "L1_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$OverdispersionConfig[1], n_i=DataList$n_v, list_names=c("L1_z","eta1_vf"), sd=0 )
    if( "L2_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$OverdispersionConfig[2], n_i=DataList$n_v, list_names=c("L2_z","eta2_vf"), sd=0 )
  }
  # Fields
  if( "L_omega1_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$FieldConfig[1], n_i=DataList$n_s, list_names=c("L_omega1_z","Omegainput1_sf"), sd=0 )
  if( "L_epsilon1_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$FieldConfig[2], n_i=DataList$n_s, n_t=DataList$n_t, list_names=c("L_epsilon1_z","Epsiloninput1_sft"), sd=0 )
  if( "L_omega2_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$FieldConfig[3], n_i=DataList$n_s, list_names=c("L_omega2_z","Omegainput2_sf"), sd=0 )
  if( "L_epsilon2_z" %in% names(Return)) Return = Add_factor( List=Return, n_c=DataList$n_c, n_f=DataList$FieldConfig[4], n_i=DataList$n_s, n_t=DataList$n_t, list_names=c("L_epsilon2_z","Epsiloninput2_sft"), sd=0 )
  # Initial values
  if( length(DataList$ObsModel)==1 || DataList$ObsModel[2]==0 ){
    Return[["beta1_ct"]] = qlogis(0.01*0.99*tapply(ifelse(DataList$b_i>0,1,0),INDEX=factor(DataList$c_i,levels=sort(unique(DataList$c_i))),FUN=mean)) %o% rep(1,DataList$n_t)
    Return[["beta2_ct"]] = log(tapply(ifelse(DataList$b_i>0,DataList$b_i/DataList$a_i,NA),INDEX=factor(DataList$c_i,levels=sort(unique(DataList$c_i))),FUN=mean,na.rm=TRUE)) %o% rep(1,DataList$n_t)
  }
  if( length(DataList$ObsModel)==2 && DataList$ObsModel[2]%in%c(1,2) ){
    Return[["beta1_ct"]] = array(0, dim=c(DataList$n_c,DataList$n_t))
    Return[["beta2_ct"]] = array(0, dim=c(DataList$n_c,DataList$n_t))
  }

  # If either beta or epsilon is a random-walk process, fix starting value at 1
  if( "Beta_rho1"%in%names(Return) && RhoConfig[["Beta1"]]==2 ) Return[["Beta_rho1"]] = 1
  if( "Beta_rho2"%in%names(Return) && RhoConfig[["Beta2"]]==2 ) Return[["Beta_rho2"]] = 1
  if( "Epsilon_rho1"%in%names(Return) && RhoConfig[["Epsilon1"]]==2 ) Return[["Epsilon_rho1"]] = 1
  if( "Epsilon_rho2"%in%names(Return) && RhoConfig[["Epsilon2"]]==2 ) Return[["Epsilon_rho2"]] = 1
  # replace missing values function
  tmpfn = function( vec ){
    Return = ifelse( abs(vec)==Inf, NA, vec)
    return( ifelse(is.na(Return),mean(Return,na.rm=TRUE),Return) )
  }
  Return[["beta1_ct"]] = tmpfn( Return[["beta1_ct"]] )
  Return[["beta2_ct"]] = tmpfn( Return[["beta2_ct"]] )
  # Error messages
  if( any(sapply(Return, FUN=function(num){any(is.na(num))})) ) stop("Some parameter is NA")
  # Return tagged list
  return( Return )
}
