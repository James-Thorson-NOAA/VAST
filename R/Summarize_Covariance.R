
#' @export
Summarize_Covariance = function( report, tmbdata, parhat, sd_report=NULL, species_order=1:tmbdata$n_c, names_set=1:tmbdata$n_c, figname=NULL, plotTF=c("Omega1"=TRUE,"Epsilon1"=TRUE,"Omega2"=TRUE,"Epsilon2"=TRUE), plot_cor=TRUE, mgp=c(2,0.5,0), tck=-0.02, oma=c(0,5,2,0)){

  # Object to return
  Return = list()

  # Extract standard errors
  if( !is.null(sd_report) ){
    # Object to build
    sd_summary = summary(sd_report)

    # Extract
    for(i in 1:4){
      Par_name = c("omega1", "epsilon1", "omega2", "epsilon2")[i]
      Slot_name = paste0("lowercov_uppercor_",Par_name)
      if( Slot_name %in% rownames(sd_summary) ){
        # Extract covariances
        Cor = Cov = Mat = ThorsonUtilities::Extract_SE( SD=Sdreport, parname=Slot_name, columns=1:2, Dim=c(tmbdata$n_c,tmbdata$n_c) )
        dimnames(Cor) = dimnames(Cov) = list( names_set, names_set, c("Estimate","Std.Error") )
        # Cor
        Cor[,,1][lower.tri(Cor[,,1])] = t(Mat[,,1])[lower.tri(Mat[,,1])]
        diag(Cor[,,1]) = 1
        Cor[,,2][lower.tri(Cor[,,2])] = t(Mat[,,2])[lower.tri(Mat[,,2])]
        diag(Cor[,,2]) = NA
        # Cov
        Cov[,,1][upper.tri(Cov[,,1])] = t(Mat[,,1])[upper.tri(Mat[,,1])]
        Cov[,,2][upper.tri(Cov[,,2])] = t(Mat[,,2])[upper.tri(Mat[,,2])]
        # Add to return
        List = list( Cor, Cov )
        names(List) = paste0( c("Cor_","Cov_"),Par_name)
        Return = c( Return, List )
      }
    }
  }else{
    Return = vector('list',length=8)
    names(Return) = paste0( rep(c("Cor_","Cov_"),4), rep(c("omega1", "epsilon1", "omega2", "epsilon2"),each=2) )
  }

  # Plot covariances
  if( !is.null(figname) ){
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
    png( file=paste0(figname,"--Analytic.png"), width=Dim[2]*4+1, height=Dim[1]*4, units="in", res=200 )
      par(mfrow=Dim, mar=c(0,1,1,0), mgp=mgp, tck=tck, oma=oma)
      for(i in 1:4 ){
        Cov_cc = NULL
        if( i %in% which(plotTF) ){
          Cov_cc = VAST:::calc_cov( L_z=parhat[c('L_omega1_z','L_epsilon1_z','L_omega2_z','L_epsilon2_z')][[i]], n_f=tmbdata$FieldConfig[i], n_c=tmbdata$n_c )
          VAST:::plot_cov( Cov=convert(Cov_cc)[species_order,species_order], names=list(names_set[species_order],NA)[[ifelse(i==1|i==3|Dim[2]==1,1,2)]], names2=list(1:nrow(Cov_cc),NA)[[ifelse(i==1|i==2,1,2)]], digits=1, font=2 )
          if(i==1 | Dim[1]==1) mtext(side=3, text="Spatial", line=1.5, font=2)
          if(i==2 | Dim[1]==1) mtext(side=3, text="Spatio-temporal", line=1.5, font=2)
          if(i==2 | (Dim[2]==1&i==1)) mtext(side=4, text=ifelse(length(tmbdata$ObsModel)==1||tmbdata$ObsModel[2]==0,"Encounter probability","Component #1"), line=0.5, font=2)
          if(i==4 | (Dim[2]==1&i==3)) mtext(side=4, text=ifelse(length(tmbdata$ObsModel)==1||tmbdata$ObsModel[2]==0,"Positive catch rate","Component #2"), line=0.5, font=2)
        }
        if( length(Return[[paste0( "Cov_", c("omega1", "epsilon1", "omega2", "epsilon2")[i])]])==0 ){
          Return[[paste0( "Cov_", c("omega1", "epsilon1", "omega2", "epsilon2")[i])]] = Cov_cc
          if( !is.null(Cov_cc)) Return[[paste0( "Cor_", c("omega1", "epsilon1", "omega2", "epsilon2")[i])]] = cov2cor(Cov_cc)
        }
      }
    dev.off()

    # Plot sample
    png( file=paste0(figname,"--Sample.png"), width=Dim[2]*4+1, height=Dim[1]*4, units="in", res=200 )
      par(mfrow=Dim, mar=c(0,1,1,0), mgp=mgp, tck=tck, oma=oma)
      for(i in which(plotTF) ){
        if(i==1) Cov_cc = cov(report$Omega1_sc)
        if(i==2) Cov_cc = cov(apply(report$Epsilon1_sct,MARGIN=2,FUN=as.vector))
        if(i==3) Cov_cc = cov(report$Omega2_sc)
        if(i==4) Cov_cc = cov(apply(report$Epsilon2_sct,MARGIN=2,FUN=as.vector))
        VAST:::plot_cov( Cov=convert(Cov_cc)[species_order,species_order], names=list(names_set[species_order],NA)[[ifelse(i==1|i==3|Dim[2]==1,1,2)]], names2=list(1:nrow(Cov_cc),NA)[[ifelse(i==1|i==2,1,2)]], digits=1, font=2 )
        if(i==1 | Dim[1]==1) mtext(side=3, text="Spatial", line=1.5, font=2)
        if(i==2 | Dim[1]==1) mtext(side=3, text="Spatio-temporal", line=1.5, font=2)
        if(i==2 | (Dim[2]==1&i==1)) mtext(side=4, text=ifelse(length(tmbdata$ObsModel)==1||tmbdata$ObsModel[2]==0,"Encounter probability","Component #1"), line=0.5, font=2)
        if(i==4 | (Dim[2]==1&i==3)) mtext(side=4, text=ifelse(length(tmbdata$ObsModel)==1||tmbdata$ObsModel[2]==0,"Positive catch rate","Component #2"), line=0.5, font=2)
      }
    dev.off()
  }

  # Return
  return( invisible(Return) )
}

