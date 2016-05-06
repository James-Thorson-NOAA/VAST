

Plot_Covariance = function( figname, report, tmbdata, parhat, names_set=1:tmbData$n_c, plotTF=c("Omega1"=TRUE,"Epsilon1"=TRUE,"Omega2"=TRUE,"Epsilon2"=TRUE), plot_cor=TRUE, mgp=c(2,0.5,0), tck=-0.02, oma=c(0,5,2,0)){
  # Work out dimensions
  Dim = c(2,2)
  if( sum(plotTF)==1 ) Dim = c(1,1)
  if( all(plotTF==c(1,1,0,0)) | all(plotTF==c(0,0,1,1)) ) Dim=c(1,2)
  if( all(plotTF==c(1,0,1,0)) | all(plotTF==c(0,1,0,1)) ) Dim=c(2,1)

  # Conversion function
  if(plot_cor==TRUE){
    convert = function(Cov) ifelse(is.na(cov2cor(Cov)),0,cov2cor(Cov))
  }else{
    convert = function(Cov) ifelse(is.na(Cov),0,Cov)
  }

  # Plot analytic
  save_fig( file=paste0(figname,"--Analytic"), width=Dim[2]*4+1, height=Dim[1]*4 )
    par(mfrow=Dim, mar=c(0,1,1,0), mgp=mgp, tck=tck, oma=oma)
    for(i in which(plotTF) ){
      Cov_cc = calc_cov( L_z=parhat[c('L_omega1_z','L_epsilon1_z','L_omega2_z','L_epsilon2_z')][[i]], n_f=tmbdata$FieldConfig[i], n_c=tmbdata$n_c )
      plot_cov( Cov=convert(Cov_cc), names=list(names_set,NA)[[ifelse(i==1|i==3|Dim[2]==1,1,2)]], names2=list(1:nrow(Cov_cc),NA)[[ifelse(i==1|i==2,1,2)]], digits=1, font=2 )
      if(i==1 | Dim[1]==1) mtext(side=3, text="Spatial", line=1.5, font=2)
      if(i==2 | Dim[1]==1) mtext(side=3, text="Spatio-temporal", line=1.5, font=2)
      if(i==2 | (Dim[2]==1&i==1)) mtext(side=4, text="Component #1", line=0.5, font=2)
      if(i==4 | (Dim[2]==1&i==3)) mtext(side=4, text="Component #2", line=0.5, font=2)
    }
  dev.off()

  # Plot covariances
  save_fig( file=paste0(figname,"--Sample"), width=Dim[2]*4+1, height=Dim[1]*4 )
    par(mfrow=Dim, mar=c(0,1,1,0), mgp=mgp, tck=tck, oma=oma)
    for(i in which(plotTF) ){
      if(i==1) Cov_cc = cov(report$Omega1_sc)
      if(i==2) Cov_cc = cov(apply(report$Epsilon1_sct,MARGIN=2,FUN=as.vector))
      if(i==3) Cov_cc = cov(report$Omega2_sc)
      if(i==4) Cov_cc = cov(apply(report$Epsilon2_sct,MARGIN=2,FUN=as.vector))
      plot_cov( Cov=convert(Cov_cc), names=list(names_set,NA)[[ifelse(i==1|i==3|Dim[2]==1,1,2)]], names2=list(1:nrow(Cov_cc),NA)[[ifelse(i==1|i==2,1,2)]], digits=1, font=2 )
      if(i==1 | Dim[1]==1) mtext(side=3, text="Spatial", line=1.5, font=2)
      if(i==2 | Dim[1]==1) mtext(side=3, text="Spatio-temporal", line=1.5, font=2)
      if(i==2 | (Dim[2]==1&i==1)) mtext(side=4, text="Component #1", line=0.5, font=2)
      if(i==4 | (Dim[2]==1&i==3)) mtext(side=4, text="Component #2", line=0.5, font=2)
    }
  dev.off()
}

