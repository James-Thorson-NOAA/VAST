
#' Plot DHARMa residuals
#'
#' @inheritParams DHARMa::plotResiduals
#'
#' @param dharmaRes Output from \code{\link{plot}} using slot \code{$dharmaRes}
#'
#' @author Zach Oyafuso
#' @export
plot_dharma <-
function( dharmaRes,
          form = NULL,
          rank = TRUE,
          ... ){

  par( mfrow=c(1,2), mgp=c(2,0.5,0), tck=-0.02, mar=c(4,3,2,0), yaxs="r", xaxs="r" )
  #gap::qqunif(dharmaRes$scaledResiduals,
  #            pch = 2,
  #            bty = "n",
  #            logscale = F,
  #            col = "black",
  #            cex = 0.6,
  #            cex.main = 1,
  #            ann = F,
  #            cex.axis = 1.5)
  ECDF = stats::ecdf( dharmaRes$scaledResiduals )
  plot( ECDF,
        xlim=c(0,1),
        ylim=c(0,1),
        xlab = "",
        ylab = "",
        main = "" )
  abline( a=0, b=1, lty="solid" )
  mtext(side = 1, line = 2, text = "Observed")
  mtext(side = 2, line = 1.5, text = "Expected")
  mtext(side = 3, line = 0.5, text = "Quantile-quantile plot")
  box()

  DHARMa::plotResiduals(dharmaRes,
                        rank = rank,
                        ann = FALSE,
                        xlim = c(0, 1),
                        cex = 1.5,
                        form = form )
  mtext(side = 1, line = 2.5, text = "Rank-Transformed\nModel Predictions")
  mtext(side = 2, line = 1.5, text = "Standardized Residual")
  #box(which = "figure")
}
