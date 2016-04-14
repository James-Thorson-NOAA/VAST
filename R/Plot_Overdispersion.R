
Plot_Overdispersion = function( filename1, filename2, Data, ParHat, Report, SD=NULL, Map=NULL, ControlList1=list("Width"=8, "Height"=4, "Res"=200, "Units"='in'), ControlList2=list("Width"=8, "Height"=8, "Res"=200, "Units"='in') ){
  if( !("ThorsonUtilities" %in% installed.packages()[,1]) ){
    devtools::install_github("james-thorson/utilities")
  }

  Derived_Quants = NULL
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
    # Save Ls
    Derived_Quants[["L1_cf"]] = Report[["L1_cf"]]
    if(!is.null(SD) & !is.null(Map)) Derived_Quants[["SE_L1_cf"]] = ThorsonUtilities::SE_hat_fn( SD=SD, Map=Map, Dim=dim(Report[["L1_ct"]]), parname="L1_cf")
    Derived_Quants[["var_L1_cf"]] = colSums( Derived_Quants[["L1_cf"]]^2 )
    Derived_Quants[["prop_var_L1_cf"]] = cumsum(sort(Derived_Quants[["var_L1_cf"]],decreasing=TRUE)) / sum(Derived_Quants[["var_L1_cf"]])
    Derived_Quants[["L2_cf"]] = Report[["L1_cf"]]
    if(!is.null(SD) & !is.null(Map)) Derived_Quants[["SE_L2_cf"]] = ThorsonUtilities::SE_hat_fn( SD=SD, Dim=dim(Report[["L1_ct"]]), Map=Map, parname="L2_cf")
    Derived_Quants[["var_L2_cf"]] = colSums( Derived_Quants[["L2_cf"]]^2 )
    Derived_Quants[["prop_var_L2_cf"]] = cumsum(sort(Derived_Quants[["var_L2_cf"]],decreasing=TRUE)) / sum(Derived_Quants[["var_L2_cf"]])
    # Save overdispersion
    if( !("eta1_vc" %in% names(Report)) ){
      Report[["eta1_vc"]] = ParHat[["eta1_vf"]] %*% t(Report[["L1_cf"]])
      Report[["eta2_vc"]] = ParHat[["eta2_vf"]] %*% t(Report[["L2_cf"]])
    }
    Derived_Quants[["cov_eta1_vc"]] = cov( Report[["eta1_vc"]] )
    Derived_Quants[["cov_eta2_vc"]] = cov( Report[["eta2_vc"]] )
  }

  if( Data[["n_f_input"]]>=0 ){
    # Plot covariance
    png( file=paste0(filename1,".png"), width=ControlList1$Width, height=ControlList1$Height, res=ControlList1$Res, units=ControlList1$Units )
      par( mfrow=c(2,2), mar=c(0,2,3,0), mgp=c(2,0.5,0), tck=-0.02, oma=c(0,0,0,0))
      for(i in 1:4){
        Cov_cc = list(Cov1_cc, Cov2_cc, Derived_Quants[["cov_eta1_vc"]], Derived_Quants[["cov_eta2_vc"]])[[i]]
        plot_cov( Cov=Cov_cc )
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
    Return = list("Cov1_cc"=Cov1_cc, "Cov2_cc"=Cov2_cc, "Derived_Quants"=Derived_Quants)
    return( invisible(Return) )
  }
}
