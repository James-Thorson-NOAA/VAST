
Plot_Overdispersion = function( filename, Data, ParHat, Report, ControlList=list("Width"=8, "Height"=4, "Res"=200, "Units"='in') ){
  if( Data[["n_f_input"]]<0 ){
    message("No overdispersion in model")
  }
  if( Data[["n_f_input"]]==0 ){
    Dist_cc = outer( 1:Data[["n_c"]], 1:Data[["n_c"]], FUN=function(a,b){abs(a-b)})
    Cov1_cc = ParHat[["L1_z"]] ^ Dist_cc
    Cov2_cc = ParHat[["L2_z"]] ^ Dist_cc
  }
  if( Data[["n_f_input"]]>0 ){
    Cov1_cc = Report[["L1_cf"]] %*% t(Report[["L1_cf"]])
    Cov2_cc = Report[["L2_cf"]] %*% t(Report[["L2_cf"]])
  }

  if( Data[["n_f_input"]]>=0 ){
    # Plot
    png( file=filename, width=ControlList$Width, height=ControlList$Height, res=ControlList$Res, units=ControlList$Units )
      par( mfrow=c(1,2), mar=c(0,2,3,0), mgp=c(2,0.5,0), tck=-0.02, oma=c(0,0,0,0))
      for(i in 1:2){
        Cov_cc = list(Cov1_cc, Cov2_cc)[[i]]
        plot_cov( Cov=Cov_cc )
        title( main=c("Encounter prob","Positive catch rate")[i], line=2 )
      }
    dev.off()

    # Return
    Return = list("Cov1_cc"=Cov1_cc, "Cov2_cc"=Cov2_cc)
    return( invisible(Return) )
  }
}
