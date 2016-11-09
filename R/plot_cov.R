
#' Plot covariance matrix
#'
#' \code{plot_cov} plots and formats a covariance matrix
#'
#' @param Cov matrix (covariance or correlation) used for plotting
#' @param zlim numeric-vector (length 2) defining bounds for color-scale of covariance
#' @param names Labels for y-axis
#' @param names2 Labels for x-axis (on top of covariance)
#' @param ncolors Number of colors for color-scale
#' @param digits Number of digits for text-labels of covariance
#' @param ... passed to \code{text} for labelling covariances

#' @export
plot_cov = function( Cov, zlim=NULL, names=1:nrow(Cov), names2=names, ncolors=21, digits=2, ... ){
  Col = colorRampPalette(colors=c("blue","white","red"))
  if(is.null(zlim) ) zlim = c(-1,1)*max(abs(Cov))
  image(z=Cov[1:nrow(Cov),nrow(Cov):1], x=seq(0,1,length=nrow(Cov)), y=seq(0,1,length=nrow(Cov)), col=Col(ncolors), xaxt="n", yaxt="n", zlim=zlim, yaxt="n", ylab="", xlab="" )
  if(length(names)>1) axis(side=2, at=seq(0,1,length=nrow(Cov)), labels=rev(names), las=1)
  if(length(names2)>1) axis(side=3, at=seq(0,1,length=nrow(Cov)), labels=names2, las=ifelse(all(nchar(names2)==1),1,2) )
  for(i in 1:nrow(Cov)){
  for(j in 1:nrow(Cov)){
    Label = formatC(Cov[i,j],digits=digits,format="f")
    Label = ThorsonUtilities::convert_to_unicode( Label, pattern="-", replacement="\u2013" )
    text( y=seq(1,0,length=nrow(Cov))[i], x=seq(0,1,length=nrow(Cov))[j], labels=Label, ...)
  }}
  box()
}
