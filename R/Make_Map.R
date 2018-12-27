#' @export
Make_Map <-
function( DataList, TmbParams, CovConfig=TRUE, DynCovConfig=TRUE, Q_Config=TRUE, RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0), Npool=0 ){

  # Local functions
  fixval_fn <- function( fixvalTF ){
    vec = rep(0,length(fixvalTF))
    vec[which(fixvalTF)] = NA
    vec[which(!is.na(vec))] = 1:sum(!is.na(vec))
    vec = factor( vec ) 
    return( vec )
  }
  seq_pos <- function( length.out, from=1 ) seq(from=from, to=length.out, length.out=max(length.out,0))

  # Extract Options and Options_vec (depends upon version)
  if( all(c("Options","Options_vec") %in% names(DataList)) ){
    Options_vec = DataList$Options_vec
    Options = DataList$Options
  }
  if( "Options_list" %in% names(DataList) ){
    Options_vec = DataList$Options_list$Options_vec
    Options = DataList$Options_list$Options
  }

  # Create tagged-list in TMB format for fixing parameters
  Map = list()

  # Anisotropy
  if(Options_vec["Aniso"]==0 | all(DataList[["FieldConfig"]] == -1)) Map[['ln_H_input']] = factor( rep(NA,2) )

  #########################
  # 1. Residual variance ("logSigmaM")
  # 2. Lognormal-Poisson overdispersion ("delta_i")
  #########################

  # Measurement error models
  # NOTE:  Uses DataList$ObsModel_ez, which either exists or is added above
  Map[["logSigmaM"]] = array(NA, dim=dim(TmbParams$logSigmaM))
  if( "delta_i" %in% names(TmbParams)){
    Map[["delta_i"]] = rep(NA, length(TmbParams[["delta_i"]]) )
  }
  for( eI in 1:DataList$n_e ){
    if(DataList$ObsModel_ez[eI,1]%in%c(0,1,2)){
      if(ncol(Map[["logSigmaM"]])==2) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( 1, NA )
      if(ncol(Map[["logSigmaM"]])==3) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( 1, NA, NA )
    }
    if(DataList$ObsModel_ez[eI,1]%in%c(5)){
      if(ncol(Map[["logSigmaM"]])==2) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( 1, 2 )
      if(ncol(Map[["logSigmaM"]])==3) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( 1, 2, NA )
      if( any(DataList$ObsModel_ez[,2]!=0) ) stop("ObsModel[1]=5 should use ObsModel[2]=0")
    }
    if(DataList$ObsModel_ez[eI,1]%in%c(6,7)){
      if(ncol(Map[["logSigmaM"]])==2) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( NA, NA )
      if(ncol(Map[["logSigmaM"]])==3) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( NA, NA, NA )
      if( any(DataList$ObsModel_ez[,2]!=0) ) stop("ObsModel[1]=6 or 7 should use ObsModel[2]=0")
    }
    if(DataList$ObsModel_ez[eI,1]%in%c(8,10)){
      if(ncol(Map[["logSigmaM"]])==2) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( 1, NA )
      if(ncol(Map[["logSigmaM"]])==3) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( 1, NA, NA )
      if( any(DataList$ObsModel_ez[,2]!=2) ) stop("ObsModel[1]=8 and ObsModel[1]=10 should use ObsModel[2]=2")
    }
    if(DataList$ObsModel_ez[eI,1]%in%c(9)){
      if(ncol(Map[["logSigmaM"]])==2) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( NA, NA )
      if(ncol(Map[["logSigmaM"]])==3) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( NA, NA, NA )
    }
    if(DataList$ObsModel_ez[eI,1]%in%c(11)){
      if(ncol(Map[["logSigmaM"]])==2) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( 1, NA )
      if(ncol(Map[["logSigmaM"]])==3) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( 1, NA, NA )
      Map[["delta_i"]][ which((DataList$e_i+1)==eI) ] = max(c(0,Map[["delta_i"]]),na.rm=TRUE) + 1:length(which((DataList$e_i+1)==eI))
    }
    if(DataList$ObsModel_ez[eI,1]%in%c(12,13)){
      if(ncol(Map[["logSigmaM"]])==2) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( NA, NA )
      if(ncol(Map[["logSigmaM"]])==3) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( NA, NA, NA )
      if( any(DataList$ObsModel_ez[,2]!=1) ) stop("ObsModel[1]=12 or 13 should use ObsModel[2]=1")
    }
    if(DataList$ObsModel_ez[eI,1]%in%c(14)){
      if(ncol(Map[["logSigmaM"]])==2) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( 1, NA )
      if(ncol(Map[["logSigmaM"]])==3) Map[["logSigmaM"]][eI,] = max(c(0,Map[["logSigmaM"]]),na.rm=TRUE) + c( 1, NA, NA )
      Map[["delta_i"]][ which((DataList$e_i+1)==eI) ] = max(c(0,Map[["delta_i"]]),na.rm=TRUE) + 1:length(which((DataList$e_i+1)==eI))
      if( any(DataList$ObsModel_ez[,2]!=1) ) stop("ObsModel[1]=14 should use ObsModel[2]=1")
    }
  }
  Map[["logSigmaM"]] = factor(Map[["logSigmaM"]])
  if( "delta_i" %in% names(TmbParams)){
    Map[["delta_i"]] = factor(Map[["delta_i"]])
  }

  #########################
  # Variance for spatial and spatio-temporal
  #########################

  # Configurations of spatial and spatiotemporal error
  if(DataList[["FieldConfig"]]['Omega1'] == -1){
    if("Omegainput1_sc" %in% names(TmbParams)) Map[["Omegainput1_sc"]] = factor( array(NA,dim=dim(TmbParams[["Omegainput1_sc"]])) )
    if("Omegainput1_sf" %in% names(TmbParams)) Map[["Omegainput1_sf"]] = factor( array(NA,dim=dim(TmbParams[["Omegainput1_sf"]])) )
    if("L_omega1_z" %in% names(TmbParams)) Map[["L_omega1_z"]] = factor( rep(NA,length(TmbParams[["L_omega1_z"]])) )
  }
  if(DataList[["FieldConfig"]]['Epsilon1'] == -1){
    if("Epsiloninput1_sct" %in% names(TmbParams)) Map[["Epsiloninput1_sct"]] = factor( array(NA,dim=dim(TmbParams[["Epsiloninput1_sct"]])) )
    if("Epsiloninput1_sft" %in% names(TmbParams)) Map[["Epsiloninput1_sft"]] = factor( array(NA,dim=dim(TmbParams[["Epsiloninput1_sft"]])) )
    if("L_epsilon1_z" %in% names(TmbParams)) Map[["L_epsilon1_z"]] = factor( rep(NA,length(TmbParams[["L_epsilon1_z"]])) )
  }
  if(DataList[["FieldConfig"]]['Omega1'] == -1 & DataList[["FieldConfig"]]['Epsilon1'] == -1){
    Map[["logkappa1"]] = factor(NA)
    if("rho_c1" %in% names(TmbParams)) Map[["rho_c1"]] = factor(NA)
  }
  if(DataList[["FieldConfig"]]['Omega2'] == -1){
    if("Omegainput2_sc" %in% names(TmbParams)) Map[["Omegainput2_sc"]] = factor( array(NA,dim=dim(TmbParams[["Omegainput2_sc"]])) )
    if("Omegainput2_sf" %in% names(TmbParams)) Map[["Omegainput2_sf"]] = factor( array(NA,dim=dim(TmbParams[["Omegainput2_sf"]])) )
    if("L_omega2_z" %in% names(TmbParams)) Map[["L_omega2_z"]] = factor( rep(NA,length(TmbParams[["L_omega2_z"]])) )
  }
  if(DataList[["FieldConfig"]]['Epsilon2'] == -1){
    if("Epsiloninput2_sct" %in% names(TmbParams)) Map[["Epsiloninput2_sct"]] = factor( array(NA,dim=dim(TmbParams[["Epsiloninput2_sct"]])) )
    if("Epsiloninput2_sft" %in% names(TmbParams)) Map[["Epsiloninput2_sft"]] = factor( array(NA,dim=dim(TmbParams[["Epsiloninput2_sft"]])) )
    if("L_epsilon2_z" %in% names(TmbParams)) Map[["L_epsilon2_z"]] = factor( rep(NA,length(TmbParams[["L_epsilon2_z"]])) )
  }
  if(DataList[["FieldConfig"]]['Omega2'] == -1 & DataList[["FieldConfig"]]['Epsilon2'] == -1){
    Map[["logkappa2"]] = factor(NA)
    if("rho_c2" %in% names(TmbParams)) Map[["rho_c2"]] = factor(NA)
  }

  # Epsilon1 -- Fixed OR White-noise OR Random walk
  if( RhoConfig["Epsilon1"] %in% c(0,1,2) ){
    if( "Epsilon_rho1" %in% names(TmbParams) ) Map[["Epsilon_rho1"]] = factor( NA )
    if( "Epsilon_rho1_f" %in% names(TmbParams) ) Map[["Epsilon_rho1_f"]] = factor( rep(NA,length(TmbParams$Epsilon_rho1_f)) )
  }
  if( RhoConfig["Epsilon1"] %in% c(4) ){
    if( "Epsilon_rho1_f" %in% names(TmbParams) ) Map[["Epsilon_rho1_f"]] = factor( rep(1,length(TmbParams$Epsilon_rho1_f)) )
  }
  # Epsilon2 -- Fixed OR White-noise OR Random walk OR mirroring Epsilon_rho1_f
  if( RhoConfig["Epsilon2"] %in% c(0,1,2,6) ){
    if( "Epsilon_rho2" %in% names(TmbParams) ) Map[["Epsilon_rho2"]] = factor( NA )
    if( "Epsilon_rho2_f" %in% names(TmbParams) ) Map[["Epsilon_rho2_f"]] = factor( rep(NA,length(TmbParams$Epsilon_rho2_f)) )
  }
  if( RhoConfig["Epsilon2"] %in% c(4) ){
    if( "Epsilon_rho2_f" %in% names(TmbParams) ) Map[["Epsilon_rho2_f"]] = factor( rep(1,length(TmbParams$Epsilon_rho2_f)) )
  }

  # fix AR across bins
  if( DataList$n_c==1 & ("rho_c1" %in% names(TmbParams)) ){
    Map[["rho_c1"]] = factor(NA)
    Map[["rho_c2"]] = factor(NA)
  }

  #########################
  # Variance for overdispersion
  #########################

  # Overdispersion parameters
  if( ("n_f_input"%in%names(DataList)) && "n_v"%in%names(DataList) && DataList[["n_f_input"]]<0 ){
    Map[["L1_z"]] = factor(rep(NA,length(TmbParams[["L1_z"]])))
    Map[["eta1_vf"]] = factor(array(NA,dim=dim(TmbParams[["eta1_vf"]])))
    Map[["L2_z"]] = factor(rep(NA,length(TmbParams[["L1_z"]])))
    Map[["eta2_vf"]] = factor(array(NA,dim=dim(TmbParams[["eta2_vf"]])))
  }
  if( ("OverdispersionConfig"%in%names(DataList)) && "n_v"%in%names(DataList) ){
    if( DataList[["OverdispersionConfig"]][1] == -1 ){
      Map[["L1_z"]] = factor(rep(NA,length(TmbParams[["L1_z"]])))
      Map[["eta1_vf"]] = factor(array(NA,dim=dim(TmbParams[["eta1_vf"]])))
    }
    if( DataList[["OverdispersionConfig"]][2] == -1 ){
      Map[["L2_z"]] = factor(rep(NA,length(TmbParams[["L2_z"]])))
      Map[["eta2_vf"]] = factor(array(NA,dim=dim(TmbParams[["eta2_vf"]])))
    }
  }


  #########################
  # Npool options
  # Overwrites SigmaM, L_omega, and L_epsilon, so must come after them
  #########################

  # Make all category-specific variances (SigmaM, Omega, Epsilon) constant for models with EncNum_a < Npool
  if( Npool>0 ){
    if( !all(DataList$FieldConfig %in% c(-2)) ){
      stop("Npool should only be specified when using 'IID' variation for `FieldConfig`")
    }
  }
  EncNum_ct = array(0, dim=c(DataList$n_c,DataList$n_t))
  for( tz in 1:ncol(DataList$t_iz) ){
    Temp = tapply( DataList$b_i, INDEX=list(factor(DataList$c_iz[,1],levels=1:DataList$n_c-1),factor(DataList$t_iz[,tz],levels=1:DataList$n_t-1)), FUN=function(vec){sum(vec>0,na.rm=TRUE)} )
    Temp = ifelse( is.na(Temp), 0, Temp )
    EncNum_ct = EncNum_ct + Temp
  }
  EncNum_c = rowSums( EncNum_ct )
  if( any(EncNum_c < Npool) ){
    pool = function(poolTF){
      Return = 1:length(poolTF)
      Return = ifelse( poolTF==TRUE, length(poolTF)+1, Return )
      return(Return)
    }
    # Change SigmaM / L_omega1_z / L_omega2_z / L_epsilon1_z / L_epsilon2_z
    Map[["logSigmaM"]] = array( as.numeric(Map$logSigmaM), dim=dim(TmbParams$logSigmaM) )
    Map[["logSigmaM"]][ which(EncNum_c < Npool), ] = rep(1,sum(EncNum_c<Npool)) %o% Map[["logSigmaM"]][ which(EncNum_c < Npool)[1], ]
    Map[["logSigmaM"]] = factor( Map[["logSigmaM"]] )
    # Change Omegas
    Map[["L_omega1_z"]] = factor(pool(EncNum_c<Npool))
    Map[["L_omega2_z"]] = factor(pool(EncNum_c<Npool))
    # Change Epsilons
    Map[["L_epsilon1_z"]] = factor(pool(EncNum_c<Npool))
    Map[["L_epsilon2_z"]] = factor(pool(EncNum_c<Npool))
  }

  #########################
  # Covariates
  #########################

  # Static covariates
  Var_j = apply( DataList[["X_xj"]], MARGIN=2, FUN=var )
  Map[["gamma1_j"]] = Map[["gamma2_j"]] = 1:ncol(DataList$X_xj)
  for(j in 1:length(Var_j)){
    if( Var_j[j]==0 || sum(CovConfig)==0 ){
      Map[["gamma1_j"]][j] = NA
      Map[["gamma2_j"]][j] = NA
    }
  }
  Map[["gamma1_j"]] = factor(Map[["gamma1_j"]])
  Map[["gamma2_j"]] = factor(Map[["gamma1_j"]])

  # Catchability variables
  Var_k = apply( DataList[["Q_ik"]], MARGIN=2, FUN=var )
  Map[["lambda1_k"]] = Map[["lambda2_k"]] = 1:ncol(DataList$Q_ik)
  for(k in 1:length(Var_k)){
    if( Var_k[k]==0 || sum(Q_Config)==0 ){
      Map[["lambda1_k"]][k] = NA
      Map[["lambda2_k"]][k] = NA
    }
  }
  Map[["lambda1_k"]] = factor(Map[["lambda1_k"]])
  Map[["lambda2_k"]] = factor(Map[["lambda2_k"]])

  # Dynamic covariates
  if( "X_xtp" %in% names(DataList) ){
    Var_p = apply( DataList[["X_xtp"]], MARGIN=3, FUN=function(array){var(as.vector(array))})
    Var_tp = apply( DataList[["X_xtp"]], MARGIN=2:3, FUN=var )
    if( "gamma1_tp" %in% names(TmbParams) ){
      Map[["gamma1_tp"]] = Map[["gamma2_tp"]] = matrix( 1:(DataList$n_t*DataList$n_p), nrow=DataList$n_t, ncol=DataList$n_p )
      # By default:
      #  1.  turn off coefficient associated with variable having no variance across space and time
      #  2.  assume constant coefficient for all years of each variable and category
      for(p in 1:length(Var_p)){
        if( Var_p[p]==0 || sum(DynCovConfig)==0 ){
          Map[["gamma1_tp"]][,p] = NA
          Map[["gamma2_tp"]][,p] = NA
        }else{
          Map[["gamma1_tp"]][,p] = rep( Map[["gamma1_tp"]][1,p], DataList$n_t )
          Map[["gamma2_tp"]][,p] = rep( Map[["gamma2_tp"]][1,p], DataList$n_t )
        }
      }
      Map[["gamma1_tp"]] = factor(Map[["gamma1_tp"]])
      Map[["gamma2_tp"]] = factor(Map[["gamma2_tp"]])
    }
    if( "gamma1_ctp" %in% names(TmbParams) ){
      Map[["gamma1_ctp"]] = Map[["gamma2_ctp"]] = array( 1:(DataList$n_c*DataList$n_t*DataList$n_p), dim=c(DataList$n_c,DataList$n_t,DataList$n_p) )
      # By default:
      #  1.  turn off coefficient associated with variable having no variance across space and time
      #  2.  assume constant coefficient for all years of each variable and category
      for(p in 1:length(Var_p)){
        if( Var_p[p]==0 || sum(DynCovConfig)==0 ){
          Map[["gamma1_ctp"]][,,p] = NA
          Map[["gamma2_ctp"]][,,p] = NA
        }else{
          for(cI in 1:DataList$n_c){
            Map[["gamma1_ctp"]][cI,,p] = rep( Map[["gamma1_ctp"]][cI,1,p], DataList$n_t )
            Map[["gamma2_ctp"]][cI,,p] = rep( Map[["gamma2_ctp"]][cI,1,p], DataList$n_t )
          }
        }
      }
      Map[["gamma1_ctp"]] = factor(Map[["gamma1_ctp"]])
      Map[["gamma2_ctp"]] = factor(Map[["gamma2_ctp"]])
    }
  }

  #########################
  # Seasonal models
  #########################

  # fix variance-ratio for columns of t_iz
  if( "log_sigmaratio1_z" %in% names(TmbParams) ){
    Map[["log_sigmaratio1_z"]] = factor( c(NA, seq_pos(length(TmbParams[["log_sigmaratio1_z"]])-1)) )
  }
  if( "log_sigmaratio2_z" %in% names(TmbParams) ){
    Map[["log_sigmaratio2_z"]] = factor( c(NA, seq_pos(length(TmbParams[["log_sigmaratio2_z"]])-1)) )
  }

  #########################
  # Interactions
  #########################

  if( "VamConfig"%in%names(DataList) & all(c("Chi_fr","Psi_fr")%in%names(TmbParams)) ){
    # Turn off interactions
    if( DataList$VamConfig[1]==0 ){
      Map[["Chi_fr"]] = factor( rep(NA,prod(dim(TmbParams$Chi_fr))) )
      Map[["Psi_fr"]] = factor( rep(NA,prod(dim(TmbParams$Psi_fr))) )
    }
    # Reduce degrees of freedom for interactions
    if( DataList$VamConfig[1] %in% c(1,3) ){
      Map[["Psi_fr"]] = array( 1:prod(dim(TmbParams$Psi_fr)), dim=dim(TmbParams$Psi_fr) )
      Map[["Psi_fr"]][1:ncol(Map[["Psi_fr"]]),] = NA
      Map[["Psi_fr"]] = factor(Map[["Psi_fr"]])
    }
    # Reduce degrees of freedom for interactions
    if( DataList$VamConfig[1]==2 ){
      Map[["Psi_fr"]] = array( 1:prod(dim(TmbParams$Psi_fr)), dim=dim(TmbParams$Psi_fr) )
      Map[["Psi_fr"]][1:ncol(Map[["Psi_fr"]]),] = NA
      Map[["Psi_fr"]] = factor(Map[["Psi_fr"]])
      Map[["Psi_fr"]] = array( 1:prod(dim(TmbParams$Psi_fr)), dim=dim(TmbParams$Psi_fr) )
      Map[["Psi_fr"]][1:ncol(Map[["Psi_fr"]]),] = NA
      Map[["Psi_fr"]][cbind(1:ncol(Map[["Psi_fr"]]),1:ncol(Map[["Psi_fr"]]))] = max(c(0,Map[["Psi_fr"]]),na.rm=TRUE) + 1:ncol(Map[["Psi_fr"]])
      Map[["Psi_fr"]] = factor(Map[["Psi_fr"]])
    }
  }

  #########################
  # 1. Intercepts
  # 2. Hyper-parameters for intercepts
  #########################

  #####
  # Step 1: fix betas and/or epsilons for missing years if betas are fixed-effects
  #####

  Num_ct = array(0, dim=c(DataList$n_c,DataList$n_t))
  for( tz in 1:ncol(DataList$t_iz) ){
    Temp = tapply( DataList$b_i, INDEX=list(factor(DataList$c_iz[,1],levels=1:DataList$n_c-1),factor(DataList$t_iz[,tz],levels=1:DataList$n_t-1)), FUN=function(vec){sum(!is.na(vec))} )
    Temp = ifelse( is.na(Temp), 0, Temp )
    Num_ct = Num_ct + Temp
  }
  if( sum(Num_ct==0)>0 ){
    # Beta1 -- Fixed
    if( RhoConfig["Beta1"]==0){
      Map[["beta1_ct"]] = fixval_fn( fixvalTF=(Num_ct==0) )
    }else{
      # Don't fix because it would affect estimates of variance
    }
    # Beta2 -- Fixed
    if( RhoConfig["Beta2"]==0){
      Map[["beta2_ct"]] = fixval_fn( fixvalTF=(Num_ct==0) )
    }else{
      # Don't fix because it would affect estimates of variance
    }
  }

  #####
  # Step 2: User settings for 100% encounter rates
  # overwrite previous, but also re-checks for missing data
  #####

  # Change beta1_ct if 100% encounters (not designed to work with seasonal models)
  if( any(DataList$ObsModel_ez[,2] %in% c(3)) ){
    if( ncol(DataList$t_iz)==1 ){
      Tmp_ct = tapply(ifelse(DataList$b_i>0,1,0), INDEX=list(factor(DataList$c_iz[,1],levels=sort(unique(DataList$c_iz[,1]))),factor(DataList$t_iz[,1],levels=1:DataList$n_t-1)), FUN=mean)
      Map[["beta1_ct"]] = array( 1:prod(dim(Tmp_ct)), dim=dim(Tmp_ct) )
      Map[["beta1_ct"]][which(is.na(Tmp_ct) | Tmp_ct==1)] = NA
      Map[["beta1_ct"]] = factor(Map[["beta1_ct"]])
    }else{
      stop("`ObsModel[,2]==3` is not implemented to work with seasonal models")
    }
  }

  # Change beta1_ct and beta2_ct if 0% or 100% encounters (not designed to work with seasonal models)
  if( any(DataList$ObsModel_ez[,2] %in% c(4)) ){
    if( ncol(DataList$t_iz)==1 ){
      Tmp_ct = tapply(ifelse(DataList$b_i>0,1,0), INDEX=list(factor(DataList$c_iz[,1],levels=sort(unique(DataList$c_iz[,1]))),factor(DataList$t_iz[,1],levels=1:DataList$n_t-1)), FUN=mean)
      Map[["beta1_ct"]] = array( 1:prod(dim(Tmp_ct)), dim=dim(Tmp_ct) )
      Map[["beta1_ct"]][which(is.na(Tmp_ct) | Tmp_ct==1 | Tmp_ct==0)] = NA
      Map[["beta1_ct"]] = factor(Map[["beta1_ct"]])
      Map[["beta2_ct"]] = array( 1:prod(dim(Tmp_ct)), dim=dim(Tmp_ct) )
      Map[["beta2_ct"]][which(is.na(Tmp_ct) | Tmp_ct==0)] = NA
      Map[["beta2_ct"]] = factor(Map[["beta2_ct"]])
    }else{
      stop("`ObsModel[,2]==3` is not implemented to work with seasonal models")
    }
  }

  #####
  # Step 3: Structure for hyper-parameters
  # overwrites previous structure on intercepts only if temporal structure is specified (in which case its unnecessary)
  #####

  # Hyperparameters for intercepts for <= V5.3.0
  if( "Beta_mean1" %in% names(TmbParams) ){
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
      Map[["beta1_ct"]] = factor( 1:DataList$n_c %o% rep(1,DataList$n_t) )
    }
    # Beta2 -- Fixed (0) or Beta_rho2 mirroring Beta_rho1 (6)
    if( RhoConfig["Beta2"] %in% c(0,6) ){
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
      Map[["beta2_ct"]] = factor( 1:DataList$n_c %o% rep(1,DataList$n_t) )
    }
    # Warnings
    if( DataList$n_c >= 2 ){
      warnings( "This version of VAST has the same hyperparameters for the intercepts of all categories.  Please use CPP version >=5.4.0 for different hyperparameters for each category." )
    }
  }
  # Hyperparameters for intercepts for >= V5.4.0
  if( "Beta_mean1_c" %in% names(TmbParams) ){
    if( RhoConfig["Beta1"]==0){
      Map[["Beta_mean1_c"]] = factor( rep(NA,DataList$n_c) )
      Map[["Beta_rho1_c"]] = factor( rep(NA,DataList$n_c) )
      Map[["logsigmaB1_c"]] = factor( rep(NA,DataList$n_c) )
    }
    # Beta1 -- White-noise
    if( RhoConfig["Beta1"]==1){
      Map[["Beta_rho1_c"]] = factor( rep(NA,DataList$n_c) )
    }
    # Beta1 -- Random-walk
    if( RhoConfig["Beta1"]==2){
      Map[["Beta_mean1_c"]] = factor( rep(NA,DataList$n_c) )
      Map[["Beta_rho1_c"]] = factor( rep(NA,DataList$n_c) )
    }
    # Beta1 -- Constant over time for each category
    if( RhoConfig["Beta1"]==3){
      Map[["Beta_mean1_c"]] = factor( rep(NA,DataList$n_c) )
      Map[["Beta_rho1_c"]] = factor( rep(NA,DataList$n_c) )
      Map[["logsigmaB1_c"]] = factor( rep(NA,DataList$n_c) )
      Map[["beta1_ct"]] = factor( 1:DataList$n_c %o% rep(1,DataList$n_t) )
    }
    # Beta2 -- Fixed (0) or Beta_rho2 mirroring Beta_rho1 (6)
    if( RhoConfig["Beta2"] %in% c(0,6) ){
      Map[["Beta_mean2_c"]] = factor( rep(NA,DataList$n_c) )
      Map[["Beta_rho2_c"]] = factor( rep(NA,DataList$n_c) )
      Map[["logsigmaB2_c"]] = factor( rep(NA,DataList$n_c) )
    }
    # Beta2 -- White-noise
    if( RhoConfig["Beta2"]==1){
      Map[["Beta_rho2_c"]] = factor( rep(NA,DataList$n_c) )
    }
    # Beta2 -- Random-walk
    if( RhoConfig["Beta2"]==2){
      Map[["Beta_mean2_c"]] = factor( rep(NA,DataList$n_c) )
      Map[["Beta_rho2_c"]] = factor( rep(NA,DataList$n_c) )
    }
    # Beta2 -- Constant over time for each category
    if( RhoConfig["Beta2"]==3){
      Map[["Beta_mean2_c"]] = factor( rep(NA,DataList$n_c) )
      Map[["Beta_rho2_c"]] = factor( rep(NA,DataList$n_c) )
      Map[["logsigmaB2_c"]] = factor( rep(NA,DataList$n_c) )
      Map[["beta2_ct"]] = factor( 1:DataList$n_c %o% rep(1,DataList$n_t) )
    }
    # Warnings
    if( DataList$n_c >= 2 ){
      warnings( "This version of VAST has different hyperparameters for each category. Default behavior for CPP version <=5.3.0 was to have the same hyperparameters for the intercepts of all categories." )
    }
  }

  #####
  # Step 4: Structure for seasonal models
  # NOT WELL TESTED!
  #####

  # fix first level of 2nd and higher columns of t_iz
  if( "t_iz"%in%names(DataList) && ncol(DataList$t_iz)>=2 ){
    # (Re)start map for intercepts
    if( !("beta1_ct" %in% names(Map)) ){
      Map[["beta1_ct"]] = 1:prod(dim(TmbParams[["beta1_ct"]]))
    }else{
      Map[["beta1_ct"]] = as.numeric(Map[["beta1_ct"]])
    }
    if( !("beta2_ct" %in% names(Map)) ){
      Map[["beta2_ct"]] = 1:prod(dim(TmbParams[["beta2_ct"]]))
    }else{
      Map[["beta2_ct"]] = as.numeric(Map[["beta2_ct"]])
    }
    # Add fixed values for lowest value of 2nd and higher columns
    for( zI in 2:ncol(DataList$t_iz) ){
      Which2Fix = min( DataList$t_iz[,zI] )
      Which2Fix = matrix( 1:(DataList$n_c*DataList$n_t), ncol=DataList$n_t, nrow=DataList$n_c )[,Which2Fix+1]
      Map[["beta1_ct"]][Which2Fix] = NA
      Map[["beta2_ct"]][Which2Fix] = NA
    }
    # Remake as factor
    Map[["beta1_ct"]] = factor(Map[["beta1_ct"]])
    Map[["beta2_ct"]] = factor(Map[["beta2_ct"]])
  }

  # Return
  return(Map)
}

