Make_Map <-
function( Version, TmbData, TmbParams, CovConfig=TRUE, DynCovConfig=TRUE, Q_Config=TRUE, RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0), Aniso=TRUE){
  # Local functions
  fixval_fn <- function( fixvalTF ){
    vec = rep(0,length(fixvalTF))
    vec[which(fixvalTF)] = NA
    vec[which(!is.na(vec))] = 1:sum(!is.na(vec))
    vec = factor( vec ) 
    return( vec )
  }
  
  # Create tagged-list in TMB format for fixing parameters
  Map = list()
  # Configurations of spatial and spatiotemporal error
  if(TmbData[["FieldConfig"]]['Omega1']==0){
    Map[["Omegainput1_sc"]] = factor(matrix(NA,nrow=TmbData$n_s,ncol=TmbData$n_c))
    Map[["logetaO1"]] = factor(NA) 
  }
  if(TmbData[["FieldConfig"]]['Epsilon1']==0){
    Map[["Epsiloninput1_sct"]] = factor(array(NA,dim=c(TmbData$n_s,TmbData$n_c,TmbData$n_t)))
    Map[["logetaE1"]] = factor(NA) 
  }
  if(TmbData[["FieldConfig"]]['Omega1']==0 & TmbData[["FieldConfig"]]['Epsilon1']==0){
    Map[["logkappa1"]] = factor(NA)
    Map[["rho_c1"]] = factor(NA)
  }
  if(TmbData[["FieldConfig"]]['Omega2']==0){
    Map[["Omegainput2_sc"]] = factor(matrix(NA,nrow=TmbData$n_s,ncol=TmbData$n_c))
    Map[["logetaO2"]] = factor(NA) 
  }
  if(TmbData[["FieldConfig"]]['Epsilon2']==0){
    Map[["Epsiloninput2_sct"]] = factor(array(NA,dim=c(TmbData$n_s,TmbData$n_c,TmbData$n_t)))
    Map[["logetaE2"]] = factor(NA) 
  }
  if(TmbData[["FieldConfig"]]['Omega2']==0 & TmbData[["FieldConfig"]]['Epsilon2']==0){
    Map[["logkappa2"]] = factor(NA)
    Map[["rho_c2"]] = factor(NA)
  }

  # Measurement error models
  if(TmbData[["ObsModel"]]%in%c(0,1,2)){
    Map[["logSigmaM"]] = factor( c(1,NA) )
  }
  if(TmbData[["ObsModel"]]%in%c(6)){
    Map[["logSigmaM"]] = factor( c(NA,NA) )
  }
  # Anisotropy
  if(Aniso==0 | all(TmbData[["FieldConfig"]]==0)) Map[['ln_H_input']] = factor( rep(NA,2) )
  
  # Beta1 -- Fixed
  if( RhoConfig["Beta1"]==0 && Version%in%c("comp_index_v1b","comp_index_v1a")){
    Map[["Beta_mean1"]] = factor( NA )
    Map[["Beta_rho1"]] = factor( NA )
    Map[["logsigmaB1"]] = factor( NA )
  }
  # Beta1 -- White-noise
  if( RhoConfig["Beta1"]==1 && Version%in%c("comp_index_v1b","comp_index_v1a")){
    Map[["Beta_rho1"]] = factor( NA )
  }
  # Beta1 -- Random-walk
  if( RhoConfig["Beta1"]==2 && Version%in%c("comp_index_v1b","comp_index_v1a")){
    Map[["Beta_mean1"]] = factor( NA )
    Map[["Beta_rho1"]] = factor( NA )
  }
  # Beta1 -- Constant
  if( RhoConfig["Beta1"]==3 && Version%in%c("comp_index_v1b","comp_index_v1a")){
    Map[["Beta_mean1"]] = factor( NA )
    Map[["Beta_rho1"]] = factor( NA )
    Map[["logsigmaB1"]] = factor( NA )
    Map[["beta1_ct"]] = factor( array(1,dim=c(TmbData$n_c,TmbData$n_t)) )
  }
  # Beta2 -- Fixed
  if( RhoConfig["Beta2"]==0 && Version%in%c("comp_index_v1b","comp_index_v1a")){
    Map[["Beta_mean2"]] = factor( NA )
    Map[["Beta_rho2"]] = factor( NA )
    Map[["logsigmaB2"]] = factor( NA )
  }
  # Beta2 -- White-noise
  if( RhoConfig["Beta2"]==1 && Version%in%c("comp_index_v1b","comp_index_v1a")){
    Map[["Beta_rho2"]] = factor( NA )
  }
  # Beta2 -- Random-walk
  if( RhoConfig["Beta2"]==2 && Version%in%c("comp_index_v1b","comp_index_v1a")){
    Map[["Beta_mean2"]] = factor( NA )
    Map[["Beta_rho2"]] = factor( NA )
  }
  # Beta2 -- Constant
  if( RhoConfig["Beta2"]==3 && Version%in%c("comp_index_v1b","comp_index_v1a")){
    Map[["Beta_mean2"]] = factor( NA )
    Map[["Beta_rho2"]] = factor( NA )
    Map[["logsigmaB2"]] = factor( NA )
    Map[["beta2_ct"]] = factor( array(1,dim=c(TmbData$n_c,TmbData$n_t)) )
  }
  # Epsilon1 -- Fixed OR White-noise OR Random walk
  if( RhoConfig["Epsilon1"]%in%c(0,1,2) && Version%in%c("comp_index_v1b","comp_index_v1a")){
    Map[["Epsilon_rho1"]] = factor( NA )
  }
  # Epsilon2 -- Fixed OR White-noise OR Random walk
  if( RhoConfig["Epsilon2"]%in%c(0,1,2) && Version%in%c("comp_index_v1b","comp_index_v1a")){
    Map[["Epsilon_rho2"]] = factor( NA )
  }
  # fix betas and/or epsilons for missing years if betas are fixed-effects
  YearNotInData = !( (1:TmbData$n_t) %in% (unique(TmbData$t_i)+1) ) 
  if( sum(YearNotInData)>0 ){
    stop("Haven't built missing years feature yet")
    # Beta1 -- Fixed
    if( RhoConfig["Beta1"]==0 && Version%in%c("comp_index_v1b","comp_index_v1a")){
      Map[["beta1_t"]] = fixval_fn( fixvalTF=YearNotInData )
    }
    # Beta1 -- White-noise
    if( RhoConfig["Beta1"]==1 && Version%in%c("comp_index_v1b","comp_index_v1a")){
      Map[["beta1_t"]] = fixval_fn( fixvalTF=YearNotInData )
    }
    # Beta2 -- Fixed
    if( RhoConfig["Beta2"]==0 && Version%in%c("comp_index_v1b","comp_index_v1a")){
      Map[["beta2_t"]] = fixval_fn( fixvalTF=YearNotInData ) 
    }
    # Beta2 -- White-noise
    if( RhoConfig["Beta2"]==1 && Version%in%c("comp_index_v1b","comp_index_v1a")){
      Map[["beta2_t"]] = fixval_fn( fixvalTF=YearNotInData ) 
    }
  }
  # fix AR across bins
  if( TmbData$n_c==1 ){
    Map[["rho_c1"]] = factor(NA)
    Map[["rho_c2"]] = factor(NA)
  }

  # Static covariates
  Var_j = apply( TmbData[["X_xj"]], MARGIN=2, FUN=var )
  Map[["gamma1_j"]] = Map[["gamma2_j"]] = 1:TmbData$n_j
  for(j in 1:length(Var_j)){
    if( Var_j[j]==0 || sum(CovConfig)==0 ){
      Map[["gamma1_j"]][j] = NA
      Map[["gamma2_j"]][j] = NA
    }
  }
  Map[["gamma1_j"]] = factor(Map[["gamma1_j"]])
  Map[["gamma2_j"]] = factor(Map[["gamma1_j"]])

  # Catchability variables
  Var_k = apply( TmbData[["Q_ik"]], MARGIN=2, FUN=var )
  Map[["lambda1_k"]] = Map[["lambda2_k"]] = 1:TmbData$n_k
  for(k in 1:length(Var_k)){
    if( Var_k[k]==0 || sum(Q_Config)==0 ){
      Map[["lambda1_k"]][k] = NA
      Map[["lambda2_k"]][k] = NA
    }
  }
  Map[["lambda1_k"]] = factor(Map[["lambda1_k"]])
  Map[["lambda2_k"]] = factor(Map[["lambda2_k"]])

  # Dynamic covariates
  if( "X_xtp" %in% names(TmbData) ){
    if( "gamma1_tp" %in% names(TmbParams) ){
      Var_tp = apply( TmbData[["X_xtp"]], MARGIN=2:3, FUN=var )
      Map[["gamma1_tp"]] = Map[["gamma2_tp"]] = matrix( 1:(TmbData$n_t*TmbData$n_p), nrow=TmbData$n_t, ncol=TmbData$n_p )
      for(t in 1:nrow(Var_tp)){
      for(p in 1:ncol(Var_tp)){
        if( Var_tp[t,p]==0 || sum(DynCovConfig)==0 ){
          Map[["gamma1_tp"]][t,p] = NA
          Map[["gamma2_tp"]][t,p] = NA
        }
      }}
      # By default, assume constant coefficient for all years of each variable
      for(p in 1:ncol(Var_tp)){
        if( all(Var_tp[,p]>0) ){
          Map[["gamma1_tp"]][,p] = rep( Map[["gamma1_tp"]][1,p], TmbData$n_t )
          Map[["gamma2_tp"]][,p] = rep( Map[["gamma2_tp"]][1,p], TmbData$n_t )
        }
      }
      Map[["gamma1_tp"]] = factor(Map[["gamma1_tp"]])
      Map[["gamma2_tp"]] = factor(Map[["gamma2_tp"]])
    }
    if( "gamma1_ctp" %in% names(TmbParams) ){
      Var_tp = apply( TmbData[["X_xtp"]], MARGIN=2:3, FUN=var )
      Map[["gamma1_ctp"]] = Map[["gamma2_ctp"]] = array( 1:(TmbData$n_c*TmbData$n_t*TmbData$n_p), dim=c(TmbData$n_c,TmbData$n_t,TmbData$n_p) )
      for(t in 1:nrow(Var_tp)){
      for(p in 1:ncol(Var_tp)){
        if( Var_tp[t,p]==0 || sum(DynCovConfig)==0 ){
          Map[["gamma1_ctp"]][,t,p] = NA
          Map[["gamma2_ctp"]][,t,p] = NA
        }
      }}
      # By default, assume constant coefficient for all years of each variable
      for(p in 1:ncol(Var_tp)){
        if( all(Var_tp[,p]>0) ){
          Map[["gamma1_ctp"]][,,p] = rep( Map[["gamma1_ctp"]][1,1,p], TmbData$n_t )
          Map[["gamma2_ctp"]][,,p] = rep( Map[["gamma2_ctp"]][1,1,p], TmbData$n_t )
        }
      }
      Map[["gamma1_ctp"]] = factor(Map[["gamma1_ctp"]])
      Map[["gamma2_ctp"]] = factor(Map[["gamma2_ctp"]])
    }
  }

  # Return
  return(Map)
}

