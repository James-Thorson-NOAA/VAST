
#' Plot estimated overdispersion
#'
#' \code{Plot_Overdispersion} plots sample and analytic estimates of overdispersion
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
Plot_Overdispersion = function( filename1, filename2, Data, ParHat, Report, SD=NULL, Map=NULL, ControlList1=list("Width"=8, "Height"=4, "Res"=200, "Units"='in'), ControlList2=ControlList1 ){
  if( require(ThorsonUtilities)==FALSE ) devtools::install_github("james-thorson/utilities")

  Derived_Quants = NULL
  if( Data[["n_f_input"]]<0 ){
    message("No overdispersion in model")
  }
  if( Data[["n_f_input"]]>=0 ){
    Cov1_cc = VAST:::calc_cov( n_f=Data[["n_f_input"]], n_c=Data[["n_c"]], L_z=ParHat[["L1_z"]])
    Cov2_cc = VAST:::calc_cov( n_f=Data[["n_f_input"]], n_c=Data[["n_c"]], L_z=ParHat[["L2_z"]])
    Derived_Quants[["cov_eta1_vc"]] = cov( Report[["eta1_vc"]] )
    Derived_Quants[["cov_eta2_vc"]] = cov( Report[["eta2_vc"]] )
  }

  if( Data[["n_f_input"]]>=0 ){
    # Plot covariance
    png( file=paste0(filename1,".png"), width=ControlList1$Width, height=ControlList1$Height, res=ControlList1$Res, units=ControlList1$Units )
      par( mfrow=c(2,2), mar=c(0,2,3,0), mgp=c(2,0.5,0), tck=-0.02, oma=c(0,0,0,0))
      for(i in 1:4){
        Cov_cc = list(Cov1_cc, Cov2_cc, Derived_Quants[["cov_eta1_vc"]], Derived_Quants[["cov_eta2_vc"]])[[i]]
        VAST:::plot_cov( Cov=Cov_cc )
        title( main=c("Encounter prob (pop.)","Positive catch rate (pop.)","Encounter prob (sample)","Positive catch rate (sample)")[i], line=2 )
      }
    dev.off()

    # Plot distributions
    for(i in 1:2){
      png( file=paste0(filename2,"_",i,".png"), width=ControlList2$Width, height=ControlList2$Height, res=ControlList2$Res, units=ControlList2$Units )
        Mat = list( Report[["eta1_vc"]], Report[["eta2_vc"]] )[[i]]
        par( mfrow=rep(Data$n_c,2), mar=c(0,0,0,0), mgp=c(2,0.5,0), tck=-0.02, oma=c(0,0,0,0))
        for(i in 1:2){
          for(j1 in 1:Data$n_c){
          for(j2 in 1:Data$n_c){
            if(j1==j2) hist( Mat[,j1], xlab="", ylab="", xaxt="n", yaxt="n" )
            if(j1!=j2) plot( x=Mat[,j1], y=Mat[,j2], xlab="", ylab="", xaxt="n", yaxt="n" )
          }}
        }
      dev.off()
    }

    # Return
    Return = list("Cov1_cc"=Cov1_cc, "Cov2_cc"=Cov2_cc)
    return( invisible(Return) )
  }
}
