
#' Rotate results
#'
#' \code{rotate_factors} rotates results from a factor model
#'
#' @param L_pj Loadings matrix for `p` categories and `j` factors (calculated from \code{Cov_pp} if it is provided)
#' @param Psi_sjt Array of factors (1st dimension: spatial knots;  2nd dimension: factors;  3rd dimension:  time)
#' @param Cov_pp Covariance calculated from loadings matrix.
#'        if \code{Cov_pp} is provided and \code{L_pj=NULL}, then \code{L_pj} is calculated from \code{Cov_pp}
#' @param Psi_spt Array of projected factors (1st dimension: spatial knots;  2nd dimension: categories;  3rd dimension:  time)
#'        if \code{Psi_spt} is provided and \code{Psi_sjt=NULL}, then \code{Psi_sjt} is calculated from
#'        \code{Psi_spt} and \code{L_pj}
#' @param RotationMethod Method used for rotation when visualing factor decomposition results,
#'        Options: "PCA" (recommended) or "Varimax"
#' @param testcutoff tolerance for numerical rounding when confirming that rotation doesn't effect results

#' @return tagged list of outputs
#' \describe{
#'   \item{L_pj_rot}{Loadings matrix after rotation}
#'   \item{Psi_rot}{Factors after rotation}
#'   \item{Hinv}{Object used for rotation}
#'   \item{L_pj}{Loadings matrix}
#' }

#' @export
rotate_factors <-
function( L_pj = NULL,
          Psi_sjt = NULL,
          Cov_pp = NULL,
          Psi_spt = NULL,
          RotationMethod = "PCA",
          testcutoff = 1e-10,
          quiet = FALSE,
          ... ){

  # Deprecated inputs
  deprecated_inputs = list(...)
  if( "Cov_jj" %in% names(deprecated_inputs) ){
    Cov_pp = deprecated_inputs$Cov_jj
    warning("argument `Cov_jj` in `rotate_factors` has been renamed `Cov_pp` for naming consistency")
    print(Cov_pp)
  }

  # If missing time, add a third dimension
  if( length(dim(Psi_sjt))==2 ){
    Psi_sjt = array( Psi_sjt, dim=c(dim(Psi_sjt),1) )
  }
  if( length(dim(Psi_spt))==2 ){
    Psi_spt = array( Psi_spt, dim=c(dim(Psi_spt),1) )
  }

  # Local functions
  approx_equal = function(m1,m2,denominator=mean(m1+m2),d=1e-10) (2*abs(m1-m2)/denominator) < d
  trunc_machineprec = function(n) ifelse(n<1e-10,0,n)
  #Nyears = nrow(L_pj)

  # Default option: Use L_pj and Psi_sjt
  if( !is.null(L_pj) ){
    #if(quiet==FALSE) message("Using `L_pj` for loadings matrix")
    Nfactors = ncol(L_pj)
    # Drop singular factors
    #Cov_pp = L_pj %*% t(L_pj)
    #Nfactors = sum(abs(eigen(Cov_pp)$values)>1e-6)
    #L_pj = L_pj[,1:Nfactors,drop=FALSE]
    #Psi_sjt = Psi_sjt[,1:Nfactors,,drop=FALSE]
  }else{
    # Backup:  Calculate L_pj from Cov_pp, and Psi_sjt from Psi_spt
    if( !is.null(Cov_pp) ){
      if(quiet==FALSE) message( "Re-calculating `L_pj` from `Cov_pp`")
    }else{
      stop( "Must provide either `L_pj` or `Cov_pp`" )
    }
    if( any(abs(Cov_pp-t(Cov_pp))>1e-6) ) stop("`Cov_pp` does not appear to be symmetric")
    Nfactors = sum(abs(eigen(Cov_pp)$values)>1e-6)
    #Nfactors = nrow(Cov_pp)
    SVD = svd(Cov_pp)    # SVD$u %*% diag(SVD$d) %*% t(SVD$v) AND SVD$u == t(SVD$v) for a diagonal Cov_pp
    L_pj = SVD$u %*% diag(sqrt(SVD$d))[,1:Nfactors,drop=FALSE]

    # Default option: Use Psi_sjt
      # Backup:  Calculate Psi_sjt from Psi_spt
    if( !is.null(Psi_sjt) ){
      stop("Cannot supply `Cov_pp` and `Psi_sjt` to `rotate_factors`")
    }else{
      if( !is.null(Psi_spt) ){
        if(quiet==FALSE) message( "Re-calculating `Psi_sjt` from `Psi_spt` and `L_pj`")
        Psi_sjt = array(NA,dim=c(dim(Psi_spt)[1],Nfactors,dim(Psi_spt)[3]))
        for( tI in 1:dim(Psi_spt)[3] ){
          tmp = Psi_spt[,,tI] %*% SVD$u %*% diag(1/sqrt(SVD$d))
          Psi_sjt[,,tI] = tmp[,1:Nfactors,drop=FALSE]
        }
      }
    }
  }
  Nknots = dim(Psi_sjt)[1]

  # Varimax
  if( RotationMethod=="Varimax" ){
    Hinv = varimax( L_pj, normalize=FALSE )
    H = solve(Hinv$rotmat)
    L_pj_rot = L_pj %*% Hinv$rotmat
    if( !is.null(Psi_sjt) ){
      Psi_rot = array(NA, dim=dim(Psi_sjt))
      # My new factors
      for( n in 1:Nknots ){
        Psi_rot[n,,] = H %*% Psi_sjt[n,,]
      }
    }
  }

  # PCA
  if( RotationMethod=="PCA" ){
    Cov_tmp = L_pj%*%t(L_pj)
    Cov_tmp = 0.5*Cov_tmp + 0.5*t(Cov_tmp) # Avoid numerical issues with complex eigen-decomposition due to numerical underflow
    Eigen = eigen(Cov_tmp)
    Eigen$values_proportion = Eigen$values / sum(Eigen$values)
    Eigen$values_cumulative_proportion = cumsum(Eigen$values) / sum(Eigen$values)
    # Check decomposition
    #all(approx_equal( Eigen$vectors%*%diag(Eigen$values)%*%t(Eigen$vectors), L_pj%*%t(L_pj)))
    # My attempt at new loadings matrix
    L_pj_rot = (Eigen$vectors%*%diag(sqrt(trunc_machineprec(Eigen$values))))[,1:Nfactors,drop=FALSE]
    rownames(L_pj_rot) = rownames(L_pj)
    # My new factors
    H = corpcor::pseudoinverse(L_pj_rot) %*% L_pj
    Hinv = list("rotmat"=corpcor::pseudoinverse(H))
    if( !is.null(Psi_sjt) ){
      Psi_rot = array(NA, dim=dim(Psi_sjt))
      for( n in 1:Nknots ){
        Psi_rot[n,,] = H %*% Psi_sjt[n,,]
      }
    }
  }

  # Flip around
  for( j in 1:dim(L_pj_rot)[2] ){
    Sign = sign(sum(L_pj_rot[,j]))
    # sign(0) = 0, so switch 0 to 1 in that case
    Sign = ifelse(Sign==0, 1, Sign)
    if( !is.null(Psi_sjt) ){
      Psi_rot[,j,] = Psi_rot[,j,] * Sign
    }
    L_pj_rot[,j] = L_pj_rot[,j] * Sign
  }

  # Check for errors
  # Check covariance matrix
    # Should be identical for rotated and unrotated
  if( !is.na(testcutoff) ){
    if( !all(approx_equal(L_pj%*%t(L_pj),L_pj_rot%*%t(L_pj_rot), d=testcutoff, denominator=1)) ){
      Diff = L_pj%*%t(L_pj) - L_pj_rot%*%t(L_pj_rot)
      stop("Covariance matrix is changed by rotation; maximum absolute difference = ", max(abs(Diff)) )
    }
    # Check linear predictor
      # Should give identical predictions as unrotated
    if( !is.null(Psi_sjt) ){
      for(i in 1:dim(Psi_sjt)[[1]]){
      for(j in 1:dim(Psi_sjt)[[3]]){
        MaxDiff = max(L_pj%*%Psi_sjt[i,,j] - L_pj_rot%*%Psi_rot[i,,j])
        if( !all(approx_equal(L_pj%*%Psi_sjt[i,,j],L_pj_rot%*%Psi_rot[i,,j], d=testcutoff, denominator=1)) ){
          stop(paste0("Linear predictor is wrong for site ",i," and time ",j," with difference ",MaxDiff))
        }
      }}
    }
    # Check rotation matrix
      # Should be orthogonal (R %*% transpose = identity matrix) with determinant one
      # Doesn't have det(R) = 1; determinant(Hinv$rotmat)!=1 ||
    Diag = Hinv$rotmat %*% t(Hinv$rotmat)
    diag(Diag) = ifelse( abs(diag(Diag))<testcutoff, 1, diag(Diag) )
    if( !all(approx_equal(Diag,diag(Nfactors), d=testcutoff, denominator=1)) ){
      message("Rotation Hinv times its transpose:")
      print( Diag )
      stop("Error: Rotation matrix is not a rotation")
    }
  }

  # Return stuff
  Return = list( "L_pj_rot"=L_pj_rot, "Hinv"=Hinv, "L_pj"=L_pj )
  if(RotationMethod=="PCA") Return[["Eigen"]] = Eigen
  if( !is.null(Psi_sjt) ){
    Return[["Psi_rot"]] = Psi_rot
  }
  return( Return )
}
