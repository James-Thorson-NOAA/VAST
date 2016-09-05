
#' @export
plot_eigen = function( Cov, which2plot=1:min(3,ncol(Cov)), names=1:nrow(Cov), digits=2, las=2, add=FALSE, ... ){
  Eigen = eigen(Cov)
  if(add==FALSE) par( mfrow=par()$mfrow, oma=par()$oma, mar=par()$mar, tck=par()$tck )
  for( cI in which2plot ){
    plot( 1, type="n", xlim=c(0.5,ncol(Cov)+0.5), ylim=c(-1,1.2), xlab="", ylab="", xaxt="n", ... )
    for( rI in 1:nrow(Cov)){
      lines(y=c(0,Eigen$vectors[rI,cI])*sign(Eigen$vectors[which.max(abs(Eigen$vectors[,cI])),cI]), x=rep(rI,2), lwd=5)
    }
    abline( h=0, lty="dotted" )
    legend( "top", bty="n", legend=paste0("Eigenvalue = ",formatC(Eigen$values[cI],format="f",digits=2)) )
    axis(1, at=1:ncol(Cov), labels=names, las=las  )
  }
  return( invisible(Eigen) )
}
