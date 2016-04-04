Param_Fn <-
function( Version, DataList, RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0) ){
  #
  rarray = function( dim, mean=0, sd=0.01 ) array( rnorm(prod(dim),mean=mean,sd=sd), dim=dim)
  # 
  if(Version%in%c("comp_index_v1a")){
    Return = list("ln_H_input"=c(0,0), "beta1_ct"=qlogis(0.01*0.99*tapply(ifelse(DataList$b_i>0,1,0),INDEX=list(factor(DataList$c_i,levels=sort(unique(DataList$c_i))),factor(DataList$t_i,levels=1:DataList$n_t)),FUN=mean)), "gamma1_j"=rep(0,DataList$n_j), "gamma1_tp"=matrix(0,nrow=DataList$n_t,ncol=DataList$n_p), "lambda1_k"=rep(0,DataList$n_k), "logetaE1"=0, "logetaO1"=0, "logkappa1"=0, "Beta_mean1"=0, "logsigmaB1"=log(1), "Beta_rho1"=0, "Epsilon_rho1"=0, "rho_c1"=0, "Omegainput1_sc"=matrix(0,DataList$n_s,DataList$n_c), "Epsiloninput1_sct"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_t)), "beta2_ct"=log(tapply(ifelse(DataList$b_i>0,DataList$b_i/DataList$a_i,NA),INDEX=list(factor(DataList$c_i,levels=sort(unique(DataList$c_i))),factor(DataList$t_i,levels=1:DataList$n_t)),FUN=mean,na.rm=TRUE)), "gamma2_j"=rep(0,DataList$n_j), "gamma2_tp"=matrix(0,nrow=DataList$n_t,ncol=DataList$n_p), "lambda2_k"=rep(0,DataList$n_k), "logetaE2"=0, "logetaO2"=0, "logkappa2"=0, "Beta_mean2"=0, "logsigmaB2"=log(1), "Beta_rho2"=0, "Epsilon_rho2"=0, "rho_c2"=0, "logSigmaM"=c(log(5),log(2)), "Omegainput2_sc"=matrix(0,DataList$n_s,DataList$n_c), "Epsiloninput2_sct"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_t)) )
  }
  if(Version%in%c("comp_index_v1b")){
    Return = list("ln_H_input"=c(0,0), "beta1_ct"=qlogis(0.01*0.99*tapply(ifelse(DataList$b_i>0,1,0),INDEX=list(factor(DataList$c_i,levels=sort(unique(DataList$c_i))),factor(DataList$t_i,levels=1:DataList$n_t)),FUN=mean)), "gamma1_j"=rep(0,DataList$n_j), "gamma1_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda1_k"=rep(0,DataList$n_k), "logetaE1"=0, "logetaO1"=0, "logkappa1"=0, "Beta_mean1"=0, "logsigmaB1"=log(1), "Beta_rho1"=0, "Epsilon_rho1"=0, "rho_c1"=0, "Omegainput1_sc"=matrix(0,DataList$n_s,DataList$n_c), "Epsiloninput1_sct"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_t)), "beta2_ct"=log(tapply(ifelse(DataList$b_i>0,DataList$b_i/DataList$a_i,NA),INDEX=list(factor(DataList$c_i,levels=sort(unique(DataList$c_i))),factor(DataList$t_i,levels=1:DataList$n_t)),FUN=mean,na.rm=TRUE)), "gamma2_j"=rep(0,DataList$n_j), "gamma2_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda2_k"=rep(0,DataList$n_k), "logetaE2"=0, "logetaO2"=0, "logkappa2"=0, "Beta_mean2"=0, "logsigmaB2"=log(1), "Beta_rho2"=0, "Epsilon_rho2"=0, "rho_c2"=0, "logSigmaM"=c(log(5),log(2)), "Omegainput2_sc"=matrix(0,DataList$n_s,DataList$n_c), "Epsiloninput2_sct"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_t)) )
  }
  if(Version%in%c("comp_index_v1c")){
    Return = list("ln_H_input"=c(0,0), "beta1_ct"=qlogis(0.01*0.99*tapply(ifelse(DataList$b_i>0,1,0),INDEX=list(factor(DataList$c_i,levels=sort(unique(DataList$c_i))),factor(DataList$t_i,levels=1:DataList$n_t)),FUN=mean)), "gamma1_j"=rep(0,DataList$n_j), "gamma1_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda1_k"=rep(0,DataList$n_k), "L1_z"=NA, "logetaE1"=0, "logetaO1"=0, "logkappa1"=0, "Beta_mean1"=0, "logsigmaB1"=log(1), "Beta_rho1"=0, "Epsilon_rho1"=0, "rho_c1"=0, "eta1_vf"=NA, "Omegainput1_sc"=matrix(0,DataList$n_s,DataList$n_c), "Epsiloninput1_sct"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_t)), "beta2_ct"=log(tapply(ifelse(DataList$b_i>0,DataList$b_i/DataList$a_i,NA),INDEX=list(factor(DataList$c_i,levels=sort(unique(DataList$c_i))),factor(DataList$t_i,levels=1:DataList$n_t)),FUN=mean,na.rm=TRUE)), "gamma2_j"=rep(0,DataList$n_j), "gamma2_ctp"=array(0,dim=c(DataList$n_c,DataList$n_t,DataList$n_p)), "lambda2_k"=rep(0,DataList$n_k), "L2_z"=NA, "logetaE2"=0, "logetaO2"=0, "logkappa2"=0, "Beta_mean2"=0, "logsigmaB2"=log(1), "Beta_rho2"=0, "Epsilon_rho2"=0, "rho_c2"=0, "logSigmaM"=rep(1,DataList$n_c)%o%c(log(5),log(2)), "eta2_vf"=NA, "Omegainput2_sc"=matrix(0,DataList$n_s,DataList$n_c), "Epsiloninput2_sct"=array(0,dim=c(DataList$n_s,DataList$n_c,DataList$n_t)) )
    if( TmbData[["n_f_input"]]>0 ){
      Return[["L1_z"]] = Return[["L2_z"]] = rnorm(TmbData$n_f_input*TmbData$n_c-TmbData$n_f_input*(TmbData$n_f_input-1)/2)
      Return[["eta1_vf"]] = Return[["eta2_vf"]] = rarray(dim=c(TmbData$n_v,TmbData$n_f_input))
    }
    if( TmbData[["n_f_input"]]<=0 ){
      Return[["L1_z"]] = Return[["L2_z"]] = c(1,0.5) # Pointwise SD / Correlation
      Return[["eta1_vf"]] = Return[["eta2_vf"]] = rarray(dim=c(TmbData$n_v,TmbData$n_c))
    }
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
