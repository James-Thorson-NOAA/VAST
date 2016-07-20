
#' @export
Coherence = function( report, tmbdata, year_set=NULL ){
  # Determine defaults
  if( is.null(year_set)) year_set = 1:tmbdata$n_t

  # Derived quantities
  M_xc = apply( report$D_xct[,,year_set,drop=FALSE], MARGIN=c(1,2), FUN=mean)
  S_xc = apply( report$D_xct[,,year_set,drop=FALSE], MARGIN=c(1,2), FUN=var)
  Dsum_xt = apply( report$D_xct[,,year_set,drop=FALSE], MARGIN=c(1,3), FUN=sum)
  Msum_x = apply( Dsum_xt, MARGIN=1, FUN=mean)
  Ssum_x = apply( Dsum_xt, MARGIN=1, FUN=var)
  I_ct = apply( report$D_xct[,,year_set,drop=FALSE]*(tmbdata$a_xl[,1]%o%rep(1,tmbdata$n_c)%o%rep(1,length(year_set))), MARGIN=c(2,3), FUN=sum)
  Isum_t = apply( I_ct, MARGIN=2, FUN=sum)
  Psum_x = apply( Dsum_xt*(tmbdata$a_xl[,1]%o%rep(1,length(year_set)))/(rep(1,tmbdata$n_x)%o%Isum_t), MARGIN=c(1), FUN=mean)

  # Synchrony index
  phi_x = Ssum_x / apply(sqrt(S_xc),MARGIN=1,FUN=sum)^2
  phi = weighted.mean( phi_x, w=Psum_x )

  # Coherence index
  L_c = eigen(CovSum_hat)$values
  psi = 2 * (mean(cumsum(L_c)/sum(L_c))-0.5)

  # Return stuff
  Return = list("phi_x"=phi_x, "phi"=phi, "psi"=psi, "L_c"=L_c)
  return( Return )
}
