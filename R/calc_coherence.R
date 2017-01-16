
#' Calculate stability metrics
#'
#' \code{Coherence} calculates buffering (`phi`) and coherence (`psi`)
#'
#' @inheritParams Plot_Overdispersion
#' @inheritParams Data_Fn
#' @param covhat estimated covariance used for calculating coherence

#' @return Tagged list containing measures of synchrony
#' \describe{
#'   \item{phi_xz}{Synchrony index for each site (x) and each period (row of \code{yearbounds_zz})}
#'   \item{phi_z}{weighted-average of \code{phi_xz} for each period, weighted by average community-abundance at each site in that period}
#'   \item{psi}{Measure of proportion of variance explained by leading eigen-vectors}
#'   \item{L_c}{Cholesky decomposition of \code{covhat}}
#' }

#' @export
calc_coherence = function( Report, Data, covhat=NULL, yearbounds_zz=matrix(c(1,Data$n_t),nrow=1) ){

  ##################
  # Synchrony
  ##################
  n_z = nrow(yearbounds_zz)

  # Category-specific density ("D"), mean and variance
  varD_xcz = array( NA, dim=c(Data$n_x,Data$n_c,n_z) )

  # Community density ("C") summing across categories
  C_xt = array(0, dim=c(Data$n_x,Data$n_t))
  varC_xz = array(NA, dim=c(Data$n_x,n_z))

  # Area-weighted total biomass ("B") for each category
  B_ct = array(0, dim=c(Data$n_c,Data$n_t))
  B_t = rep(0, Data$n_t)

  # Proportion of total biomass ("P") for each station
  P_xz = array(NA, dim=c(Data$n_x,n_z))

  # Synchrony indices
  phi_xz = array(NA, dim=c(Data$n_x,n_z))
  phi_z = rep(0, n_z)

  # Temporary variables
  temp_xt = array(NA, dim=c(Data$n_x,Data$n_t))
  temp_c = rep(NA, Data$n_c)

  # Derived quantities
  for( tI in 1:Data$n_t ){
    for( cI in 1:Data$n_c ){
      for( xI in 1:Data$n_x ){
        C_xt[xI,tI] = C_xt[xI,tI] + Report$D_xct[xI,cI,tI]
        B_ct[cI,tI] = B_ct[cI,tI] + Report$D_xct[xI,cI,tI] * Data$a_xl[xI,1]
      }
      B_t[tI] = B_t[tI] + B_ct[cI,tI]
    }
  }

  # Loop through periods (only using operations available in TMB, i.e., var, mean, row, col, segment)
  for( zI in 1:n_z){
  for( xI in 1:Data$n_x){
    # Variance for each category
    for( cI in 1:Data$n_c ){
      for( tI in 1:Data$n_t ) temp_xt[xI,tI] = Report$D_xct[xI,cI,tI]
      varD_xcz[xI,cI,zI] = var(temp_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]])
    }
    # Variance for combined across categories
    varC_xz[xI,zI] = var(C_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]])
    # Proportion in each category
    P_xz[xI,zI] = Data$a_xl[xI,1] * mean( C_xt[xI,][yearbounds_zz[zI,1]:yearbounds_zz[zI,2]] / B_t[yearbounds_zz[zI,1]:yearbounds_zz[zI,2]] )
    # Synchrony index
    for( cI in 1:Data$n_c ) temp_c[cI] = varD_xcz[xI,cI,zI]
    phi_xz[xI,zI] = varC_xz[xI,zI] / sum(sqrt(temp_c))^2
    phi_z[zI] = phi_z[zI] + (phi_xz[xI,zI] * P_xz[xI,zI])
  }}

  ##################
  # Coherence
  ##################

  # Coherence index
  if( !is.null(covhat) ){
    L_c = eigen(covhat)$values
    psi = 2 * (mean(cumsum(L_c)/sum(L_c))-0.5)
  }else{
    L_c = psi = NULL
  }

  # Return stuff
  Return = list("phi_xz"=phi_xz, "phi_z"=phi_z, "psi"=psi, "P_xz"=P_xz, "L_c"=L_c)
  return( Return )
}
