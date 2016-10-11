
#' Calculate stability metrics
#'
#' \code{Coherence} calculates buffering (`phi`) and coherence (`psi`)
#'
#' @inheritParams Plot_Overdispersion
#' @param covhat estimated covariance used for calculating coherence
#' @param year_set subset of year-indices used when calculating buffering

#' @return Tagged list containing measures of synchrony
#' \describe{
#'   \item{phi_x}{Synchrony index for each site}
#'   \item{phi}{weighted-average of \code{phi_x} weighted by community-abundance at each site}
#'   \item{psi}{Measure of proportion of variance explained by leading eigen-vectors}
#'   \item{L_c}{Cholesky decomposition of \code{covhat}}
#' }

#' @export
Coherence = function( Report, Data, covhat=NULL, year_set=NULL ){
  # Determine defaults
  if( is.null(year_set)) year_set = 1:Data$n_t

  # Derived quantities
  M_xc = apply( Report$D_xct[,,year_set,drop=FALSE], MARGIN=c(1,2), FUN=mean)
  S_xc = apply( Report$D_xct[,,year_set,drop=FALSE], MARGIN=c(1,2), FUN=var)
  Dsum_xt = apply( Report$D_xct[,,year_set,drop=FALSE], MARGIN=c(1,3), FUN=sum)
  Msum_x = apply( Dsum_xt, MARGIN=1, FUN=mean)
  Ssum_x = apply( Dsum_xt, MARGIN=1, FUN=var)
  I_ct = apply( Report$D_xct[,,year_set,drop=FALSE]*(Data$a_xl[,1]%o%rep(1,Data$n_c)%o%rep(1,length(year_set))), MARGIN=c(2,3), FUN=sum)
  Isum_t = apply( I_ct, MARGIN=2, FUN=sum)
  Psum_x = apply( Dsum_xt*(Data$a_xl[,1]%o%rep(1,length(year_set)))/(rep(1,Data$n_x)%o%Isum_t), MARGIN=c(1), FUN=mean)

  # Synchrony index
  phi_x = Ssum_x / apply(sqrt(S_xc),MARGIN=1,FUN=sum)^2
  phi = weighted.mean( phi_x, w=Psum_x )

  # Coherence index
  if( !is.null(covhat) ){
    L_c = eigen(covhat)$values
    psi = 2 * (mean(cumsum(L_c)/sum(L_c))-0.5)
  }else{
    L_c = psi = NULL
  }

  # Return stuff
  Return = list("phi_x"=phi_x, "phi"=phi, "psi"=psi, "L_c"=L_c)
  return( Return )
}
