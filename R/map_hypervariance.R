
map_hypervariance = function( report ){
  solveSubset <- function(Q) {
      require(Matrix)
      require(TMB)
      L <- Cholesky(Q, super=TRUE, perm=TRUE)
      invQ <- .Call("tmb_invQ", L, PACKAGE = "TMB")
      iperm <- invPerm(L@perm + 1L)
      invQ[iperm, iperm, drop=TRUE]
  }

  par( mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
  for( i in 1:2 ){
    Sigma = solveSubset( report[[ c("Q1","Q2")[i] ]] )
    RelVar = (diag(Sigma) - min(diag(Sigma))) / max(diag(Sigma) - min(diag(Sigma)))

    Col = colorRampPalette(colors=c("darkblue","blue","white","red"))
    plot( mesh$loc[,1:2], col=Col(50)[ceiling(RelVar*50)] )
  }
}
