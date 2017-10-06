#' @export
Make_Map <-
function( TmbData, TmbParams, CovConfig=TRUE, DynCovConfig=TRUE, Q_Config=TRUE, RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0)){
  # Local functions
  fixval_fn <- function( fixvalTF ){
    vec = rep(0,length(fixvalTF))
    vec[which(fixvalTF)] = NA
    vec[which(!is.na(vec))] = 1:sum(!is.na(vec))
    vec = factor( vec ) 
    return( vec )
  }
  seq_pos <- function( length.out, from=1 ) seq(from=from, to=length.out, length.out=max(length.out,0))

  # Create tagged-list in TMB format for fixing parameters
  Map = list()
  # Configurations of spatial and spatiotemporal error
  if(TmbData[["FieldConfig"]]['Omega1'] == -1){
    if("Omegainput1_sc" %in% names(TmbParams)) Map[["Omegainput1_sc"]] = factor( array(NA,dim=dim(TmbParams[["Omegainput1_sc"]])) )
    if("Omegainput1_sf" %in% names(TmbParams)) Map[["Omegainput1_sf"]] = factor( array(NA,dim=dim(TmbParams[["Omegainput1_sf"]])) )
    if("L_omega1_z" %in% names(TmbParams)) Map[["L_omega1_z"]] = factor( rep(NA,length(TmbParams[["L_omega1_z"]])) )
  }
  if(TmbData[["FieldConfig"]]['Epsilon1'] == -1){
    if("Epsiloninput1_sct" %in% names(TmbParams)) Map[["Epsiloninput1_sct"]] = factor( array(NA,dim=dim(TmbParams[["Epsiloninput1_sct"]])) )
    if("Epsiloninput1_sft" %in% names(TmbParams)) Map[["Epsiloninput1_sft"]] = factor( array(NA,dim=dim(TmbParams[["Epsiloninput1_sft"]])) )
    if("L_epsilon1_z" %in% names(TmbParams)) Map[["L_epsilon1_z"]] = factor( rep(NA,length(TmbParams[["L_epsilon1_z"]])) )
  }
  if(TmbData[["FieldConfig"]]['Omega1'] == -1 & TmbData[["FieldConfig"]]['Epsilon1'] == -1){
    Map[["logkappa1"]] = factor(NA)
    if("rho_c1" %in% names(TmbParams)) Map[["rho_c1"]] = factor(NA)
  }
  if(TmbData[["FieldConfig"]]['Omega2'] == -1){
    if("Omegainput2_sc" %in% names(TmbParams)) Map[["Omegainput2_sc"]] = factor( array(NA,dim=dim(TmbParams[["Omegainput2_sc"]])) )
    if("Omegainput2_sf" %in% names(TmbParams)) Map[["Omegainput2_sf"]] = factor( array(NA,dim=dim(TmbParams[["Omegainput2_sf"]])) )
    if("L_omega2_z" %in% names(TmbParams)) Map[["L_omega2_z"]] = factor( rep(NA,length(TmbParams[["L_omega2_z"]])) )
  }
  if(TmbData[["FieldConfig"]]['Epsilon2'] == -1){
    if("Epsiloninput2_sct" %in% names(TmbParams)) Map[["Epsiloninput2_sct"]] = factor( array(NA,dim=dim(TmbParams[["Epsiloninput2_sct"]])) )
    if("Epsiloninput2_sft" %in% names(TmbParams)) Map[["Epsiloninput2_sft"]] = factor( array(NA,dim=dim(TmbParams[["Epsiloninput2_sft"]])) )
    if("L_epsilon2_z" %in% names(TmbParams)) Map[["L_epsilon2_z"]] = factor( rep(NA,length(TmbParams[["L_epsilon2_z"]])) )
  }
  if(TmbData[["FieldConfig"]]['Omega2'] == -1 & TmbData[["FieldConfig"]]['Epsilon2'] == -1){
    Map[["logkappa2"]] = factor(NA)
    if("rho_c2" %in% names(TmbParams)) Map[["rho_c2"]] = factor(NA)
  }

  # Measurement error models
  if(TmbData[["ObsModel"]][1]%in%c(0,1,2)){
    if(ncol(TmbParams[["logSigmaM"]])==2) Map[["logSigmaM"]] = factor( cbind(seq(1,TmbData$n_c),NA) )
    if(ncol(TmbParams[["logSigmaM"]])==3) Map[["logSigmaM"]] = factor( cbind(seq(1,TmbData$n_c),NA,NA) )
  }
  if(TmbData[["ObsModel"]][1]%in%c(5)){
    if(ncol(TmbParams[["logSigmaM"]])==2) Map[["logSigmaM"]] = factor( cbind(seq(1,TmbData$n_c),TmbData$n_c+seq(1,TmbData$n_c)) )
    if(ncol(TmbParams[["logSigmaM"]])==3) Map[["logSigmaM"]] = factor( cbind(seq(1,TmbData$n_c),TmbData$n_c+seq(1,TmbData$n_c),NA) )
    if(TmbData[["ObsModel"]][2]!=0) stop("ObsModel[1]=5 should use ObsModel[2]=0")
  }
  if(TmbData[["ObsModel"]][1]%in%c(6,7)){
    if(ncol(TmbParams[["logSigmaM"]])==2) Map[["logSigmaM"]] = factor( matrix(NA,nrow=TmbData$n_c,ncol=2) )
    if(ncol(TmbParams[["logSigmaM"]])==3) Map[["logSigmaM"]] = factor( matrix(NA,nrow=TmbData$n_c,ncol=3) )
    if(TmbData[["ObsModel"]][2]!=0) stop("ObsModel[1]=6 or 7 should use ObsModel[2]=0")
  }
  if(TmbData[["ObsModel"]][1]%in%c(8,10)){
    if(ncol(TmbParams[["logSigmaM"]])==2) Map[["logSigmaM"]] = factor( cbind(seq(1,TmbData$n_c),NA) )
    if(ncol(TmbParams[["logSigmaM"]])==3) Map[["logSigmaM"]] = factor( cbind(seq(1,TmbData$n_c),NA,NA) )
    if(TmbData[["ObsModel"]][2]!=2) stop("ObsModel[1]=8 and ObsModel[1]=10 should use ObsModel[2]=2")
  }
  if(TmbData[["ObsModel"]][1]%in%c(9)){
    if(ncol(TmbParams[["logSigmaM"]])==2) Map[["logSigmaM"]] = factor( matrix(NA, nrow=TmbData$n_c, ncol=2) )
    if(ncol(TmbParams[["logSigmaM"]])==3) Map[["logSigmaM"]] = factor( matrix(NA, nrow=TmbData$n_c, ncol=3) )
  }
  if(TmbData[["ObsModel"]][1]%in%c(11)){
    if(ncol(TmbParams[["logSigmaM"]])==2) Map[["logSigmaM"]] = factor( cbind(seq(1,TmbData$n_c),NA) )
    if(ncol(TmbParams[["logSigmaM"]])==3) Map[["logSigmaM"]] = factor( cbind(seq(1,TmbData$n_c),NA,NA) )
  }else{
    if( "delta_i" %in% names(TmbParams) ){
      Map[["delta_i"]] = factor( rep(NA,length(TmbParams[["delta_i"]])) )
    }
  }
  if( length(TmbData[["ObsModel"]])==2 && TmbData[["ObsModel"]][2]%in%c(3) ){
    Tmp_ct = tapply(ifelse(TmbData$b_i>0,1,0), INDEX=list(factor(TmbData$c_i,levels=sort(unique(TmbData$c_i))),TmbData$t_i), FUN=mean)
    Map[["beta1_ct"]] = array( 1:prod(dim(Tmp_ct)), dim=dim(Tmp_ct) )
    Map[["beta1_ct"]][which(is.na(Tmp_ct) | Tmp_ct==1)] = NA
    Map[["beta1_ct"]] = factor(Map[["beta1_ct"]])
  }
  # Anisotropy
  if(TmbData[["Options_vec"]]["Aniso"]==0 | all(TmbData[["FieldConfig"]] == -1)) Map[['ln_H_input']] = factor( rep(NA,2) )
  
  # Beta1 -- Fixed
  if( RhoConfig["Beta1"]==0){
    Map[["Beta_mean1"]] = factor( NA )
    Map[["Beta_rho1"]] = factor( NA )
    Map[["logsigmaB1"]] = factor( NA )
  }
  # Beta1 -- White-noise
  if( RhoConfig["Beta1"]==1){
    Map[["Beta_rho1"]] = factor( NA )
  }
  # Beta1 -- Random-walk
  if( RhoConfig["Beta1"]==2){
    Map[["Beta_mean1"]] = factor( NA )
    Map[["Beta_rho1"]] = factor( NA )
  }
  # Beta1 -- Constant over time for each category
  if( RhoConfig["Beta1"]==3){
    Map[["Beta_mean1"]] = factor( NA )
    Map[["Beta_rho1"]] = factor( NA )
    Map[["logsigmaB1"]] = factor( NA )
    Map[["beta1_ct"]] = factor( 1:TmbData$n_c %o% rep(1,TmbData$n_t) )
  }
  # Beta2 -- Fixed
  if( RhoConfig["Beta2"]==0){
    Map[["Beta_mean2"]] = factor( NA )
    Map[["Beta_rho2"]] = factor( NA )
    Map[["logsigmaB2"]] = factor( NA )
  }
  # Beta2 -- White-noise
  if( RhoConfig["Beta2"]==1){
    Map[["Beta_rho2"]] = factor( NA )
  }
  # Beta2 -- Random-walk
  if( RhoConfig["Beta2"]==2){
    Map[["Beta_mean2"]] = factor( NA )
    Map[["Beta_rho2"]] = factor( NA )
  }
  # Beta2 -- Constant over time for each category
  if( RhoConfig["Beta2"]==3){
    Map[["Beta_mean2"]] = factor( NA )
    Map[["Beta_rho2"]] = factor( NA )
    Map[["logsigmaB2"]] = factor( NA )
    Map[["beta2_ct"]] = factor( 1:TmbData$n_c %o% rep(1,TmbData$n_t) )
  }
  # Epsilon1 -- Fixed OR White-noise OR Random walk
  if( RhoConfig["Epsilon1"] %in% c(0,1,2)){
    Map[["Epsilon_rho1"]] = factor( NA )
  }
  # Epsilon2 -- Fixed OR White-noise OR Random walk
  if( RhoConfig["Epsilon2"] %in% c(0,1,2)){
    Map[["Epsilon_rho2"]] = factor( NA )
  }
  # fix betas and/or epsilons for missing years if betas are fixed-effects
  #YearNotInData = !( (1:TmbData$n_t) %in% (unique(TmbData$t_i)+1) )
  Num_ct = tapply( TmbData$b_i, INDEX=list(factor(TmbData$c_i,levels=1:TmbData$n_c-1),factor(TmbData$t_i[,1],levels=1:TmbData$n_t-1)), FUN=function(vec){sum(!is.na(vec))} )
  Num_ct = ifelse( is.na(Num_ct), 0, Num_ct )
  if( sum(Num_ct==0)>0 ){
    # Beta1 -- Fixed
    if( RhoConfig["Beta1"]==0){
      Map[["beta1_ct"]] = fixval_fn( fixvalTF=(Num_ct==0) )
    }
    # Beta1 -- White-noise
    if( RhoConfig["Beta1"]==1){
      # Don't fix because it would affect estimates of variance
      #Map[["beta1_ct"]] = fixval_fn( fixvalTF=rep(YearNotInData,each=TmbData$n_c) )
    }
    # Beta2 -- Fixed
    if( !("beta2_ct" %in% names(Map)) ){
      if( RhoConfig["Beta2"]==0){
        Map[["beta2_ct"]] = fixval_fn( fixvalTF=(Num_ct==0) )
      }
      # Beta2 -- White-noise
      if( RhoConfig["Beta2"]==1){
        # Don't fix because it would affect estimates of variance
        #Map[["beta2_ct"]] = fixval_fn( fixvalTF=rep(YearNotInData,each=TmbData$n_c) )
      }
    }
  }
  # fix AR across bins
  if( TmbData$n_c==1 & ("rho_c1" %in% names(TmbParams)) ){
    Map[["rho_c1"]] = factor(NA)
    Map[["rho_c2"]] = factor(NA)
  }

  # fix variance-ratio for columns of t_iz
  if( "log_sigmaratio1_z" %in% names(TmbParams) ){
    Map[["log_sigmaratio1_z"]] = factor( c(NA, seq_pos(length(TmbParams[["log_sigmaratio1_z"]])-1)) )
  }
  if( "log_sigmaratio2_z" %in% names(TmbParams) ){
    Map[["log_sigmaratio2_z"]] = factor( c(NA, seq_pos(length(TmbParams[["log_sigmaratio2_z"]])-1)) )
  }

  # fix first level of 2nd and higher columns of t_iz
  if( "t_iz"%in%names(TmbData) && ncol(TmbData$t_iz)>=2 ){
    # (Re)start map for intercepts
    if( !("beta1_ct" %in% names(Map)) ){
      Map[["beta1_ct"]] = 1:prod(dim(TmbParams[["beta1_ct"]]))
    }else{ Map[["beta1_ct"]] = as.numeric(Map[["beta1_ct"]]) }
    if( !("beta2_ct" %in% names(Map)) ){
      Map[["beta2_ct"]] = 1:prod(dim(TmbParams[["beta2_ct"]]))
    }else{ Map[["beta2_ct"]] = as.numeric(Map[["beta2_ct"]]) }
    # Add fixed values for lowest value of 2nd and higher columns
    for( zI in 2:ncol(TmbData$t_iz) ){
      Which2Fix = min( TmbData$t_iz[,zI] )
      Map[["beta1_ct"]][Which2Fix+1] = NA
      Map[["beta2_ct"]][Which2Fix+1] = NA
    }
    # Remake as factor
    Map[["beta1_ct"]] = factor(Map[["beta1_ct"]])
    Map[["beta2_ct"]] = factor(Map[["beta2_ct"]])
  }

  # Overdispersion parameters
  if( ("n_f_input"%in%names(TmbData)) && "n_v"%in%names(TmbData) && TmbData[["n_f_input"]]<0 ){
    Map[["L1_z"]] = factor(rep(NA,length(TmbParams[["L1_z"]])))
    Map[["eta1_vf"]] = factor(array(NA,dim=dim(TmbParams[["eta1_vf"]])))
    Map[["L2_z"]] = factor(rep(NA,length(TmbParams[["L1_z"]])))
    Map[["eta2_vf"]] = factor(array(NA,dim=dim(TmbParams[["eta2_vf"]])))
  }
  if( ("OverdispersionConfig"%in%names(TmbData)) && "n_v"%in%names(TmbData) ){
    if( TmbData[["OverdispersionConfig"]][1] == -1 ){
      Map[["L1_z"]] = factor(rep(NA,length(TmbParams[["L1_z"]])))
      Map[["eta1_vf"]] = factor(array(NA,dim=dim(TmbParams[["eta1_vf"]])))
    }
    if( TmbData[["OverdispersionConfig"]][2] == -1 ){
      Map[["L2_z"]] = factor(rep(NA,length(TmbParams[["L1_z"]])))
      Map[["eta2_vf"]] = factor(array(NA,dim=dim(TmbParams[["eta2_vf"]])))
    }
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
    Var_p = apply( TmbData[["X_xtp"]], MARGIN=3, FUN=function(array){var(as.vector(array))})
    Var_tp = apply( TmbData[["X_xtp"]], MARGIN=2:3, FUN=var )
    if( "gamma1_tp" %in% names(TmbParams) ){
      Map[["gamma1_tp"]] = Map[["gamma2_tp"]] = matrix( 1:(TmbData$n_t*TmbData$n_p), nrow=TmbData$n_t, ncol=TmbData$n_p )
      # By default:
      #  1.  turn off coefficient associated with variable having no variance across space and time
      #  2.  assume constant coefficient for all years of each variable and category
      for(p in 1:length(Var_p)){
        if( Var_p[p]==0 || sum(DynCovConfig)==0 ){
          Map[["gamma1_tp"]][,p] = NA
          Map[["gamma2_tp"]][,p] = NA
        }else{
          Map[["gamma1_tp"]][,p] = rep( Map[["gamma1_tp"]][1,p], TmbData$n_t )
          Map[["gamma2_tp"]][,p] = rep( Map[["gamma2_tp"]][1,p], TmbData$n_t )
        }
      }
      Map[["gamma1_tp"]] = factor(Map[["gamma1_tp"]])
      Map[["gamma2_tp"]] = factor(Map[["gamma2_tp"]])
    }
    if( "gamma1_ctp" %in% names(TmbParams) ){
      Map[["gamma1_ctp"]] = Map[["gamma2_ctp"]] = array( 1:(TmbData$n_c*TmbData$n_t*TmbData$n_p), dim=c(TmbData$n_c,TmbData$n_t,TmbData$n_p) )
      # By default:
      #  1.  turn off coefficient associated with variable having no variance across space and time
      #  2.  assume constant coefficient for all years of each variable and category
      for(p in 1:length(Var_p)){
        if( Var_p[p]==0 || sum(DynCovConfig)==0 ){
          Map[["gamma1_ctp"]][,,p] = NA
          Map[["gamma2_ctp"]][,,p] = NA
        }else{
          for(cI in 1:TmbData$n_c){
            Map[["gamma1_ctp"]][cI,,p] = rep( Map[["gamma1_ctp"]][cI,1,p], TmbData$n_t )
            Map[["gamma2_ctp"]][cI,,p] = rep( Map[["gamma2_ctp"]][cI,1,p], TmbData$n_t )
          }
        }
      }
      Map[["gamma1_ctp"]] = factor(Map[["gamma1_ctp"]])
      Map[["gamma2_ctp"]] = factor(Map[["gamma2_ctp"]])
    }
  }

  # Return
  return(Map)
}

