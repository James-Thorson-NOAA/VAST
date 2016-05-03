
calc_cov = function( L_z, n_f, n_c ){
  loadings_matrix = function( L_val, n_rows, n_cols ){
    L_rc = matrix(0, nrow=n_rows, ncol=n_cols);
    Count = 1;
    for( rI in 1:n_rows ){
    for( cI in 1:n_cols ){
      if( rI >= cI ){
        L_rc[rI,cI] = L_val[Count]
        Count = Count + 1
      }
    }}
    return( L_rc)
  }

  if( n_f==0 ){
    if( length(L_z)!=2 ) stop("Problem with length of L_z given that n_f=0")
    Dist_cc = outer( 1:n_c, 1:n_c, FUN=function(a,b){abs(a-b)})
    Cov_cc = L_z[2] ^ Dist_cc
  }
  if( n_f>0 ){
    if( length(L_z)!=sum(n_c:(n_c-n_f+1)) ) stop("Problem with length of L_z given that n_f>0")
    L_cf = loadings_matrix( L_z, n_c, n_f )
    Cov_cc = L_cf %*% t(L_cf)
  }
  Return = Cov_cc
  return( Return )
}
