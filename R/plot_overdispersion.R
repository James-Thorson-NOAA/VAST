
#' Plot estimated overdispersion
#'
#' \code{plot_overdispersion} plots sample and analytic estimates of overdispersion
#'
#' @param filename1 filename (including absolute path) for overdispersion for component #1
#' @param filename2 filename (including absolute path) for overdispersion for component #2
#' @param Data tagged list of input data
#' @param ParHat Tagged list of fitted data, e.g., from \code{ParHat <- obj$env$parList()}
#' @param Report output report, e.g., from \code{Report <- obj$report()}
#' @param SD standard deviation report, e.g., from \code{SD <- sdreport( obj )}
#' @param Map tagged list of map inputs
#' @param ControlList1 plotting size for component #1
#' @param ControlIist2 plotting size for component #2 (default: ControlList2 <- ControlList1)

#' @return Tagged list containing inputs to function Build_TMB_Fn()
#' \describe{
#'   \item{Cov1_cc}{Covariance for component #1 (typically probability of occurence)}
#'   \item{Cov2_cc}{Covariance for component #2 (typically positive catch rates)}
#' }

#' @export
plot_overdispersion <-
function( filename1,
          filename2,
          Data,
          ParHat,
          Report,
          SD = NULL,
          Map = NULL,
          ControlList1 = list("Width"=8, "Height"=4, "Res"=200, "Units"='in'),
          ControlList2 = ControlList1 ){

  # Loop through components
  for(i in 1:2){
    if( Data[["OverdispersionConfig"]][i]>=0 ){

      if(i==1){
        Cov_cc = calc_cov( n_f=Data[["OverdispersionConfig"]][i], n_c=Data[["n_c"]], L_z=ParHat[["L1_z"]])
        cov_eta_vc = cov( Report[["eta1_vc"]] )
      }
      if(i==2){
        Cov_cc = calc_cov( n_f=Data[["OverdispersionConfig"]][i], n_c=Data[["n_c"]], L_z=ParHat[["L2_z"]])
        cov_eta_vc = cov( Report[["eta2_vc"]] )
      }

      # Plot covariance
      png( file=paste0(filename1,"_",i,".png"), width=ControlList1$Width, height=ControlList1$Height, res=ControlList1$Res, units=ControlList1$Units )
        par( mfrow=c(2,1), mar=c(0,2,3,0), mgp=c(2,0.5,0), tck=-0.02, oma=c(0,0,0,0))
        for(j in 1:2){
          Cov_cc = list(Cov_cc, cov_eta_vc)[[j]]
          VAST:::plot_cov( Cov=Cov_cc )
          title( main=c("Encounter prob (pop.)","Positive catch rate (pop.)","Encounter prob (sample)","Positive catch rate (sample)")[j*2-2+i], line=2 )
        }
      dev.off()

      # Plot distributions
      png( file=paste0(filename2,"_",i,".png"), width=ControlList2$Width, height=ControlList2$Height, res=ControlList2$Res, units=ControlList2$Units )
        Mat = list( Report[["eta1_vc"]], Report[["eta2_vc"]] )[[i]]
        par( mfrow=rep(Data$n_c,2), mar=c(0,0,0,0), mgp=c(2,0.5,0), tck=-0.02, oma=c(0,0,0,0))
        for(j1 in 1:Data$n_c){
        for(j2 in 1:Data$n_c){
          if(j1==j2) hist( Mat[,j1], xlab="", ylab="", xaxt="n", yaxt="n" )
          if(j1!=j2) plot( x=Mat[,j1], y=Mat[,j2], xlab="", ylab="", xaxt="n", yaxt="n" )
        }}
      dev.off()
    }else{
      message( "No overdispersion for ", c("presence/absence","positive catch rates")[i], " component so not generating output...")
    }
  }
}
