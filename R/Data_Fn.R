Data_Fn <-
function( Version, FieldConfig, ObsModel, b_i, a_i, c_i, s_i, t_i, a_xl, v_i=rep(0,length(b_i)), n_f=-1, X_xj=NULL, X_xtp=NULL, Q_ik=NULL, MeshList, Aniso=1, R2_interpretation=0, RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0), Options=c('SD_site_density'=0,'SD_site_logdensity'=0,'Calculate_Range'=0,'Calculate_evenness'=0,'Calculate_effective_area'=0), CheckForErrors=TRUE ){

  # Determine dimensions
  n_t = max(t_i) - min(t_i) + 1
  n_c = max(c_i) + 1
  n_v = length(unique(v_i))   # If n_v=1, then turn off overdispersion later
  n_i = length(b_i)
  n_x = nrow(a_xl)
  n_l = ncol(a_xl)

  # Covariates and defaults
  if( is.null(X_xj) ) X_xj = matrix(0, nrow=n_x, ncol=1)
  if( is.null(X_xtp) ) X_xtp = array(0, dim=c(n_x,n_t,1))
  if( is.null(Q_ik) ) Q_ik = matrix(0, nrow=n_i, ncol=1)
  n_j = ncol(X_xj)
  n_p = dim(X_xtp)[3]
  n_k = ncol(Q_ik)

  # Sort out overdispersion from inputs
    # n_f ranges from 0 to n_c, if n_f=0, then use AR1 structure
  n_f_input = n_f

  # by default, add nothing as Z_xl
  if( Options['Calculate_Range']==FALSE ){
    Z_xm = matrix(0, nrow=nrow(a_xl), ncol=ncol(a_xl) ) # Size so that it works for Version 3g-3j
  }
  if(Options['Calculate_Range']==TRUE ){
    Z_xm = MeshList$loc_x
    message( "Calculating range shift for stratum #1:",colnames(a_xl[1]))
  }

  # Check for bad data entry
  if( CheckForErrors==TRUE ){
    if( !is.matrix(a_xl) | !is.matrix(X_xj) | !is.matrix(Q_ik) ) stop("a_xl, X_xj, and Q_ik should be matrices")
    if( max(s_i)-1 > MeshList$mesh$n | min(s_i)<0 ) stop("s_i exceeds bounds in MeshList")
    if( any(a_i<=0) ) stop("a_i must be greater than zero for all observations, and at least one value of a_i is not")
  }

  # Check for bad data entry
  if( CheckForErrors==TRUE ){
    if( any(c(length(b_i),length(a_i),length(c_i),length(s_i),length(t_i),length(v_i))!=n_i) ) stop("b_i, a_i, c_i, s_i, v_i, or t_i doesn't have length n_i")
    if( nrow(a_xl)!=n_x | ncol(a_xl)!=n_l ) stop("a_xl has wrong dimensions")
    if( nrow(X_xj)!=n_x | ncol(X_xj)!=n_j ) stop("X_xj has wrong dimensions")
    if( nrow(Q_ik)!=n_i | ncol(Q_ik)!=n_k ) stop("Q_ik has wrong dimensions")
  }

  # Output tagged list
  Options_vec = c("Aniso"=Aniso, "R2_interpretation"=0, "Rho_betaTF"=ifelse(RhoConfig[["Beta1"]]|RhoConfig[["Beta2"]],1,0), "Alpha"=0, "AreaAbundanceCurveTF"=0, "CMP_xmax"=30, "CMP_breakpoint"=10 )
  if(Version%in%c("comp_index_v1b","comp_index_v1a")){
    Return = list( "n_i"=n_i, "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig, "ObsModel"=ObsModel, "Options"=Options, "b_i"=b_i, "a_i"=a_i, "c_i"=c_i, "s_i"=s_i, "t_i"=t_i-min(t_i), "a_xl"=a_xl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
  }
  if(Version%in%c("comp_index_v1c")){
    Return = list( "n_i"=n_i, "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=n_t, "n_c"=n_c, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_v"=n_v, "n_f_input"=n_f_input, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig, "ObsModel"=ObsModel, "Options"=Options, "b_i"=b_i, "a_i"=a_i, "c_i"=c_i, "s_i"=s_i, "t_i"=t_i-min(t_i), "v_i"=match(v_i,sort(unique(v_i)))-1, "a_xl"=a_xl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
  }
  if( "spde" %in% names(Return) ) Return[['spde']] = inla.spde2.matern(MeshList$mesh)$param.inla[c("M0","M1","M2")]
  if( "spde_aniso" %in% names(Return) ) Return[['spde_aniso']] = list("n_s"=MeshList$spde$n.spde, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$spde$param.inla$M0, "G0_inv"=inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)) )

  # Check for NAs
  if( CheckForErrors==TRUE ){
    if( any(sapply(Return, FUN=function(Array){any(is.na(Array))})==TRUE) ) stop("Please find and eliminate the NA from your inputs") 
  }

  # Return
  return( Return )
}
