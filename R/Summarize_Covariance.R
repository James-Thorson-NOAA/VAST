
#' Explore spatio-temporal covariance
#'
#' \code{Summarize_Covariance} plots and returns spatio-temporal covariance among categories
#'
#' @inheritParams Plot_Overdispersion
#' @param category_order integer-vector to re-order categories while plotting
#' @param category_names character-vector listing name for each category
#' @param plotdir directory (absolute path) for plots
#' @param figname character for first-part of figure names
#' @param plotTF integer-vector (length 4) specifying which covariances to plot (recommended: use \code{plotTF=FieldConfig})
#' @param plot_cor Boolean, whether to plot correlation or covariance
#' @param mgp passed to \code{par}
#' @param tck passed to \code{par}
#' @param oma passed to \code{par}

#' @export
Summarize_Covariance = function( Report, Data, ParHat, SD=NULL, category_order=1:Data$n_c, category_names=1:Data$n_c,
  plotdir=paste0(getwd(),"/"), figname="Cov", plotTF=c("Omega1"=TRUE,"Epsilon1"=TRUE,"Omega2"=TRUE,"Epsilon2"=TRUE), plot_cor=TRUE,
  mgp=c(2,0.5,0), tck=-0.02, oma=c(0,5,2,0)){

  # Object to return
  Return = list()

  # Extract standard errors
  if( !is.null(SD) ){
    # Object to build
    sd_summary = summary(SD)

    # Extract
    for(i in 1:4){
      Par_name = c("omega1", "epsilon1", "omega2", "epsilon2")[i]
      Slot_name = paste0("lowercov_uppercor_",Par_name)
      if( Slot_name %in% rownames(sd_summary) ){
        # Extract covariances
        Cor = Cov = Mat = ThorsonUtilities::Extract_SE( SD=SD, parname=Slot_name, columns=1:2, Dim=c(Data$n_c,Data$n_c) )
        dimnames(Cor) = dimnames(Cov) = list( category_names, category_names, c("Estimate","Std.Error") )
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
    if( sum(ifelse(plotTF>0,1,0))==1 ) Dim = c(1,1)
    if( all(ifelse(plotTF>0,1,0)==c(1,1,0,0)) | all(ifelse(plotTF>0,1,0)==c(0,0,1,1)) ) Dim=c(1,2)
    if( all(ifelse(plotTF>0,1,0)==c(1,0,1,0)) | all(ifelse(plotTF>0,1,0)==c(0,1,0,1)) ) Dim=c(2,1)

    # Conversion function
    if(plot_cor==TRUE){
      convert = function(Cov) ifelse(is.na(cov2cor(Cov)),0,cov2cor(Cov))
    }else{
      convert = function(Cov) ifelse(is.na(Cov),0,Cov)
    }

    # Plot analytic
    png( file=paste0(plotdir,figname,"--Analytic.png"), width=Dim[2]*4+1, height=Dim[1]*4, units="in", res=200 )
      par(mfrow=Dim, mar=c(0,1,1,0), mgp=mgp, tck=tck, oma=oma)
      for(i in 1:4 ){
        Cov_cc = NULL
        if( i %in% which(plotTF>0) ){
          Cov_cc = VAST:::calc_cov( L_z=ParHat[c('L_omega1_z','L_epsilon1_z','L_omega2_z','L_epsilon2_z')][[i]], n_f=Data$FieldConfig[i], n_c=Data$n_c )
          VAST:::plot_cov( Cov=convert(Cov_cc)[category_order,category_order], names=list(category_names[category_order],NA)[[ifelse(i==1|i==3|Dim[2]==1,1,2)]], names2=list(1:nrow(Cov_cc),NA)[[ifelse(i==1|i==2,1,2)]], digits=1, font=2 )
          if(i==1 | Dim[1]==1) mtext(side=3, text="Spatial", line=1.5, font=2)
          if(i==2 | Dim[1]==1) mtext(side=3, text="Spatio-temporal", line=1.5, font=2)
          if(i==2 | (Dim[2]==1&i==1)) mtext(side=4, text=ifelse(length(Data$ObsModel)==1||Data$ObsModel[2]==0,"Encounter probability","Component #1"), line=0.5, font=2)
          if(i==4 | (Dim[2]==1&i==3)) mtext(side=4, text=ifelse(length(Data$ObsModel)==1||Data$ObsModel[2]==0,"Positive catch rate","Component #2"), line=0.5, font=2)
        }
        if( length(Return[[paste0( "Cov_", c("omega1", "epsilon1", "omega2", "epsilon2")[i])]])==0 ){
          Return[[paste0( "Cov_", c("omega1", "epsilon1", "omega2", "epsilon2")[i])]] = Cov_cc
          if( !is.null(Cov_cc)) Return[[paste0( "Cor_", c("omega1", "epsilon1", "omega2", "epsilon2")[i])]] = cov2cor(Cov_cc)
        }
      }
    dev.off()

    # Plot sample
    png( file=paste0(plotdir,figname,"--Sample.png"), width=Dim[2]*4+1, height=Dim[1]*4, units="in", res=200 )
      par(mfrow=Dim, mar=c(0,1,1,0), mgp=mgp, tck=tck, oma=oma)
      for(i in which(plotTF>0) ){
        if(i==1) Cov_cc = cov(Report$Omega1_sc)
        if(i==2) Cov_cc = cov(apply(Report$Epsilon1_sct,MARGIN=2,FUN=as.vector))
        if(i==3) Cov_cc = cov(Report$Omega2_sc)
        if(i==4) Cov_cc = cov(apply(Report$Epsilon2_sct,MARGIN=2,FUN=as.vector))
        VAST:::plot_cov( Cov=convert(Cov_cc)[category_order,category_order], names=list(category_names[category_order],NA)[[ifelse(i==1|i==3|Dim[2]==1,1,2)]], names2=list(1:nrow(Cov_cc),NA)[[ifelse(i==1|i==2,1,2)]], digits=1, font=2 )
        if(i==1 | Dim[1]==1) mtext(side=3, text="Spatial", line=1.5, font=2)
        if(i==2 | Dim[1]==1) mtext(side=3, text="Spatio-temporal", line=1.5, font=2)
        if(i==2 | (Dim[2]==1&i==1)) mtext(side=4, text=ifelse(length(Data$ObsModel)==1||Data$ObsModel[2]==0,"Encounter probability","Component #1"), line=0.5, font=2)
        if(i==4 | (Dim[2]==1&i==3)) mtext(side=4, text=ifelse(length(Data$ObsModel)==1||Data$ObsModel[2]==0,"Positive catch rate","Component #2"), line=0.5, font=2)
      }
    dev.off()
  }

  # Return
  return( invisible(Return) )
}

