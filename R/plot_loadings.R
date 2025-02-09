
#' Plotting loadings matrix
#'
#' \code{plot_loadings} plots elements of a column of the loadings matrix
#'
#' @param L_pj Loadings matrix for `p` categories and `j` factors
#' @param whichfactor Integer, giving column of the loadings matrix to plot
#' @param addtitle Boolean, whether to add a title
#' @param LabelPosition Character, where to plot labels for rows of \code{L_pj} (Options: "Right","Above","Below")
#' @param Buffer, how much to pad top and bottom of each panel
#' @param Labels, Labels for rows of L_pj

#' @return tagged list of outputs
#' \describe{
#'   \item{L_pj_rot}{Loadings matrix after rotation}
#'   \item{Psi_rot}{Factors after rotation}
#'   \item{Hinv}{Object used for rotation}
#' }

#' @export
plot_loadings <-
function( L_pj,
          Lsd_pj = NULL,
          At = 1:nrow(L_pj),
          whichfactor = 1,
          addtitle = TRUE,
          LabelPosition = "Right",
          Buffer = c(0,0.1),
          Labels = rownames(L_pj),
          Cex = 1.2,
          legend_text = "Proportion of explained variance",
          interval_width = 1,
          ... ){

  # Modify inputs
  if( is.null(Lsd_pj) ){
    Lsd_pj = array(0, dim=dim(L_pj))
  }

  # Check for bugs
  if(length(At) != nrow(L_pj)) stop("Check argument `At` in `plot_loadings(.)`")

  # Plotting window
  Ylim = range(c(L_pj+Lsd_pj, L_pj-Lsd_pj))
  Ylim[1] = ifelse( Ylim[1]>0, 0, Ylim[1] )
  Ylim[2] = ifelse( Ylim[2]<0, 0, Ylim[2] )
  Ylim = Ylim + diff(Ylim)*Buffer
  Xlim = c(-0.5,0.5) + range(At)

  # Start plot
  plot(1, type="n", xlim=Xlim, ylim=Ylim, xlab="", ylab="", xaxt="n", xaxs="i", ... )
  if(LabelPosition=="Xaxis") axis( side=1, at=1:nrow(L_pj), labels=TRUE, las=3)
  if(addtitle==TRUE) mtext( text=paste("Factor",whichfactor), side=3, line=0.1, adj=0)
  abline(h=0)

  if( all(Lsd_pj==0) ){
    # Loop through categories and plot each
    for(p in 1:nrow(L_pj)){
      lines(y=c(0,L_pj[p,whichfactor]), x=rep(At[p],2), lwd=5)
      if(LabelPosition=="Right") text(x=(At[p]+0.2), y=0+(0.02*max(L_pj[p,whichfactor])), labels=Labels[p], srt=90, pos=4, cex=Cex)
      if(LabelPosition=="Above") text(x=At[p], y=max(L_pj), labels=Labels[p], srt=90, pos=3, cex=Cex)
      if(LabelPosition=="Below") text(x=At[p], y=min(L_pj), labels=Labels[p], srt=90, pos=1, cex=Cex)
      if(LabelPosition=="Side") axis(1)
    }
    if(!is.null(legend_text)){
      legend( "top", legend=paste0(legend_text," ",round(100*sum(L_pj[,whichfactor]^2)/sum(L_pj^2),1),"%"), bty="n")
    }
  }else{
    plot_timeseries( x=At, y=L_pj[,whichfactor], y_sd=Lsd_pj[,whichfactor] ,
      interval_width=interval_width)
    axis(1)
  }
}
