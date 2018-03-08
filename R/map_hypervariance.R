
#' Plot variance of GMRF knots
#'
#' \code{map_hypervariance} Plot variance of GMRF knots
#'
#' @param report Report from TmbObj (e.g., TmbList[["Obj"]])
#' @param Spatial_List Output from SpatialDeltaGLMM::Spatial_Information_Fn()
#' @param method Choose whether to plot anisotropic or isotropic covariance matrix
#' @examples
#' map_hypervariance(report = Save$Report, 
#'                   Spatial_List = Spatial_List, 
#'                   method = "anisotropic")

#' @export
map_hypervariance <- function(report, Spatial_List, method){
  solveSubset <- function(Q) {
    require(Matrix)
    require(TMB)
    L <- Cholesky(Q, super=TRUE, perm=TRUE)
    invQ <- .Call("tmb_invQ", L, PACKAGE = "TMB")
    iperm <- invPerm(L@perm + 1L)
    invQ[iperm, iperm, drop=TRUE]
  }
  
  if (method == "anisotropic") {
    mesh <- Spatial_List$MeshList$anisotropic_mesh  
  } else if (method == "isotropic") {
    mesh <- Spatial_List$MeshList$isotropic_mesh
  }
  
  
  # Q1 is the covariance matrix for the GMRF of the encounter model
  # Q2 is the covariance matrix for the GMRF of the catch-rate model
  # Mesh is the knots associated with the GMRF
  
  diag_sigma <- matrix(data = NA, nrow = length(diag(report$Q1)), ncol = 2)
  for( i in 1:2 ){
    Sigma = solveSubset( report[[ c("Q1","Q2")[i] ]] )
    diag_sigma[,i] <- log(diag(Sigma))
  }
  
  boundary_id <- mesh$segm$bnd$idx
  df2plot0 <-
    data.frame(E_km = mesh$loc[,1],
               N_km = mesh$loc[,2],
               sigma_q1 = diag_sigma[,1],
               sigma_q2 = diag_sigma[,2],
               label = rep(c("Encounter GMRF", "Catch-rate GMRF"), 
                           each = length(RelVar[,1])))
  df2plot <-
    tidyr::gather(df2plot0, variable, value, -E_km, -N_km, -label)
  
  p <-
    ggplot2::ggplot(df2plot, ggplot2::aes(x = E_km, y = N_km)) +
    ggplot2::geom_point(ggplot2::aes(color = value)) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::facet_wrap(~label) +
    ggplot2::ggtitle("Log-variance of GMRF knots")
  
  print(p)
}
