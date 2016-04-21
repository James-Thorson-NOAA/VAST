
calc_cov = function( L_z, n_f, n_c ){
  if( n_f==0 ){
    if( length(L_z)!=2 ) stop("Problem with length of L_z given that n_f=0")
    Dist_cc = outer( 1:n_c, 1:n_c, FUN=function(a,b){abs(a-b)})
    Cov_cc = L_z[2] ^ Dist_cc
  }
  if( n_f>0 ){
    if( length(L_z)!=sum(n_c:(n_c-n_f+1)) ) stop("Problem with length of L_z given that n_f>0")
    L_cf = matrix(0, nrow=n_c, ncol=n_f)
    L_cf[lower.tri(L_cf, diag=TRUE)] = L_z
    Cov_cc = L_cf %*% t(L_cf)
  }
  Return = Cov_cc
  return( Return )
}
