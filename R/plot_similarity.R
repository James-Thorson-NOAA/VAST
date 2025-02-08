
#' Plot similarity resulting from estimated covariance
#'
#' \code{plot_similarity} plots metrics of similarity derived from estimated covariance matrices
#'
#' @param fit output from \code{fit_model}
#' @param similarity_metric approach used to visualize similarity among years/categories
#'        resulting from estimated loadings matrices.  Available options include
#'        \code{"hclust", "Correlation", "Dissimilarity", "Covariance"}
#'
#' @export
plot_similarity <-
function( fit,
          year_labels = fit$year_labels,
          category_names = fit$category_names,
          similarity_metric = c("hclust", "Correlation", "Dissimilarity", "Covariance")[1],
          working_dir = getwd(),
          file_name = similarity_metric,
          panel_size = 3,
          Res = 200,
          Format = "png" ){

  # Explore seriation::dissplot

  # Update labels
  fit$Report = amend_output( fit = fit,
                             year_labels = year_labels,
                             category_names = category_names )

  if( fit$data_list$n_c >= 2 ){
    #
    if(Format=="png"){
      png(file=file.path(working_dir,paste0(file_name,".png")),
           width=2 * panel_size,
           height=4 * panel_size, res=Res, units='in')
      on.exit( dev.off() )
    }
    if(Format=="jpg"){
      jpeg(file=file.path(working_dir,paste0(file_name,".jpg")),
           width=2 * panel_size,
           height=4 * panel_size, res=Res, units='in')
      on.exit( dev.off() )
    }
    if(Format%in%c("tif","tiff")){
      tiff(file=file.path(working_dir,paste0(file_name,".tif")),
           width = 2 * panel_size,
           height = 4 * panel_size, res=Res, units='in')
      on.exit( dev.off() )
    }
    par( mfcol=c(4,2), mgp=c(2,0.5,0), mar=c(0,4,4,0), oma=c(0,2,2,0) )
    for(Col in 1:2){
    for(Row in 1:4){
    #for(i in 1:8){

      # Variable names
      i = Row + (Col-1)*4
      Par_name = c("Omega1", "Epsilon1", "Beta1", "EpsilonTime1", "Omega2", "Epsilon2", "Beta2", "EpsilonTime2")[i]
      Lpar_name = c("L_omega1_z", "L_epsilon1_z", "L_beta1_z", "Ltime_epsilon1_z", "L_omega2_z", "L_epsilon2_z", "L_beta2_z", "Ltime_epsilon2_z")[i]

      # Backwards compatible loading of variables and names
      if(Par_name == "Omega1"){ Var_name = "Omegainput1_sf"; Var2_name = "Omegainput1_gf"; L_name = "L_omega1_cf" }
      if(Par_name == "Epsilon1"){ Var_name = "Epsiloninput1_sft"; Var2_name = "Epsiloninput1_gft"; L_name = "L_epsilon1_cf" }
      if(Par_name == "Beta1"){ Var_name = "beta1_ft"; Var2_name = "missing"; L_name = "L_beta1_cf" }
      if(Par_name == "EpsilonTime1"){ Var_name = "Epsiloninput1_sff"; Var2_name = "Epsiloninput1_gff"; L_name = "Ltime_epsilon1_tf" }
      if(Par_name == "Omega2"){ Var_name = "Omegainput2_sf"; Var2_name = "Omegainput2_gf"; L_name = "L_omega2_cf" }
      if(Par_name == "Epsilon2"){ Var_name = "Epsiloninput2_sft"; Var2_name = "Epsiloninput2_gft"; L_name = "L_epsilon2_cf" }
      if(Par_name == "Beta2"){ Var_name = "beta2_ft"; Var2_name = "missing"; L_name = "L_beta2_cf" }
      if(Par_name == "EpsilonTime2"){ Var_name = "Epsiloninput2_sff"; Var2_name = "Epsiloninput2_gff"; L_name = "Ltime_epsilon2_tf" }

      Cov = fit$Report[[L_name]] %*% t(fit$Report[[L_name]])
      Dist = dist(fit$Report[[L_name]], diag=TRUE, upper=TRUE)     #
      # equivalent to: sqrt(outer( diag(Cov), diag(Cov), "+" ) - 2*Cov)
      if( (ncol(fit$Report[[L_name]])==0) || all(Cov==diag(ncol(Cov))) ){
        diag(Cov) = 0
        Dist[] = 0
        Cor = array(1, dim=dim(Cov))
      }else{
        Cor = cov2cor(Cov)
      }

      if( nrow(Cov) <= 2 ){
        plot.new()
        #legend( "center", bty="n", legend = "Skipped: covariance is diagonal")
      }else{
        if( tolower(similarity_metric) %in% c("cor","correlation") ){
          corrplot::corrplot.mixed( Cor, tl.pos="lt" )
        } else
        if( tolower(similarity_metric) %in% c("cov","covariance") ){
          corrplot::corrplot( Cov, is.corr=FALSE, tl.pos="lt", cl.lim = range(Cov) )
        } else
        if( tolower(similarity_metric) == "hclust" ){
          # Throws error with two groups
          # X - diag(diag(X)) throws error when Cov is 1-by-1 matrix with value 0
          offdiag = Cov - diag(diag(Cov))
          if( all( offdiag==0 ) ){
            plot.new()
          }else{
            Hclust = hclust( Dist )
            plot(Hclust, main="", ylab="")
          }
        } else
        if( tolower(similarity_metric) == "dissimilarity" ){
          Order = seriation::seriate( Dist )
          Dist2 = seriation::permute( Dist, Order )
          Dist2 = as.matrix(Dist2)
          #Dist2 = Dist2 / max(abs(Dist))
          #gclus::plotcolors(dmat.color(Dist2, viridisLite::viridis(4)), rlabels=rownames(as.matrix(Dist2)) )
          corrplot::corrplot( Dist2, is.corr=FALSE, tl.pos="lt", cl.lim = range(Dist2) )
        } else { stop("Check `similarity_metric`") }
        #if( tolower(similarity_metric) == "correlation" ) corrplot::corrplot( cov2cor(Cov), method="pie", type="lower" )
        #if( tolower(similarity_metric) == "dissimilarity" ) plot( as.matrix(Dist) )
        #if( tolower(similarity_metric) == "dissimilarity" ) pimage( Dist, order=seriate(Dist), axes="x", newpage=FALSE, pop=FALSE )
      }
      if(Col==1) mtext( side=2, line=1, text=c("Omega","Epsilon","Beta","EpsilonTime")[Row] )
      if(Row==1) mtext( side=3, line=1, text=c("Component 1","Component 2")[Col] )
    }}
  }
}
