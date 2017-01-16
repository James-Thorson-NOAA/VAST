
#' Calculate synchrony among species or locations
#'
#' \code{calc_synchrony} calculates synchrony
#'
#' @inheritParams Plot_Overdispersion
#' @inheritParams Data_Fn

#' @return Tagged list containing measures of synchrony
#' \describe{
#'   \item{phi_xz}{Synchrony index for each site (x) and each period (row of \code{yearbounds_zz})}
#'   \item{phi_z}{weighted-average of \code{phi_xz} for each period, weighted by average community-abundance at each site in that period}
#' }

#' @export
calc_synchrony = function( Report, Data, yearbounds_zz=matrix(c(1,Data$n_t),nrow=1) ){

  # Index lengths
  n_z = nrow(yearbounds_zz)
  n_t = Data$n_t
  n_c = Data$n_c
  n_x = Data$n_x
  n_y = nrow(Data$t_yz)

  # Extract elements
  pow = function(a,b) a^b
  D_xcy = Report$D_xcy
  a_xl = Data$a_xl

  # Density ("D") or area-expanded total biomass ("B") for each category (use B when summing across sites)
  D_xy = matrix( 0, nrow=n_x, ncol=n_y );
  B_cy = matrix( 0, nrow=n_c, ncol=n_y );
  B_y = rep( 0, n_y );
  # Mean
  meanD_xcz = array( 0, dim=c(n_x,n_c,n_z) );
  meanD_xz = matrix( 0, nrow=n_x, ncol=n_z );
  meanB_cz = matrix( 0, nrow=n_c, ncol=n_z );
  meanB_z = rep( 0, n_z );
  # Sample variance in category-specific density ("D") and biomass ("B")
  varD_xcz = array( 0, dim=c(n_x,n_c,n_z) );
  varD_xz = matrix( 0, nrow=n_x, ncol=n_z );
  varB_cz = matrix( 0, nrow=n_c, ncol=n_z );
  varB_z = rep( 0, n_z );
  varB_xbar_z = rep( 0, n_z );
  varB_cbar_z = rep( 0, n_z );
  maxsdD_xz = matrix( 0, nrow=n_x, ncol=n_z );
  maxsdB_cz = matrix( 0, nrow=n_c, ncol=n_z );
  maxsdB_z = rep( 0, n_z );
  # Proportion of total biomass ("P") for each location or each category
  propB_xz = matrix( 0, nrow=n_x, ncol=n_z );
  propB_cz = matrix( 0, nrow=n_c, ncol=n_z );
  # Synchrony indices
  phi_xz = matrix( 0, nrow=n_x, ncol=n_z );
  phi_cz = matrix( 0, nrow=n_c, ncol=n_z );
  phi_xbar_z = rep( 0, n_z );
  phi_cbar_z = rep( 0, n_z );
  phi_z = rep( 0, n_z );

  # Calculate total biomass for different categories
  for( yI in 1:n_y ){
    for( cI in 1:n_c ){
      for( xI in 1:n_x ){
        D_xy[xI,yI] = D_xy[xI,yI] + D_xcy[xI,cI,yI];
        B_cy[cI,yI] = B_cy[cI,yI] + D_xcy[xI,cI,yI] * a_xl[xI,1];
        B_y[yI] = B_y[yI] + D_xcy[xI,cI,yI] * a_xl[xI,1];
      }
    }
  }
  # Loop through periods (only using operations available in TMB, i.e., var, mean, row, col, segment)
  for( zI in 1:n_z ){
    for( xI in 1:n_x ){
      # Variance for biomass in each category, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
      for( cI in 1:n_c ){
        temp_mean = 0;
        for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) meanD_xcz[xI,cI,zI] = meanD_xcz[xI,cI,zI] + D_xcy[xI,cI,yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
        for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ){
          varD_xcz[xI,cI,zI] = varD_xcz[xI,cI,zI] + pow(D_xcy[xI,cI,yI]-meanD_xcz[xI,cI,zI],2) / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]);
        }
      }
      # Variance for combined biomass across categories, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) meanD_xz[xI,zI] = meanD_xz[xI,zI] + D_xy[xI,yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
        varD_xz[xI,zI] = varD_xz[xI,zI] + pow(D_xy[xI,yI]-meanD_xz[xI,zI],2) / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]);
      }
    }
    for( cI in 1:n_c ){
      # Variance for combined biomass across sites, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) meanB_cz[cI,zI] = meanB_cz[cI,zI] + B_cy[cI,yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
        varB_cz[cI,zI] = varB_cz[cI,zI] + pow(B_cy[cI,yI]-meanB_cz[cI,zI],2) / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]);
      }
    }
    # Variance for combined biomass across sites and categories, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
    for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) meanB_z[zI] = meanB_z[zI] + B_y[yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
    for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
      varB_z[zI] = varB_z[zI] + pow(B_y[yI]-meanB_z[zI],2) / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]);
    }
    # Proportion in each site
    for( xI in 1:n_x ){
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
        propB_xz[xI,zI] = propB_xz[xI,zI] + a_xl[xI,1] * D_xy[xI,yI] / B_y[yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
      }
    }
    # Proportion in each category
    for( cI in 1:n_c ){
      for( yI in yearbounds_zz[zI,1]:yearbounds_zz[zI,2] ) {
        propB_cz[cI,zI] = propB_cz[cI,zI] + B_cy[cI,yI] / B_y[yI] / (yearbounds_zz[zI,2]-yearbounds_zz[zI,1]+1);
      }
    }
    # Species-buffering index (calculate in Density so that areas with zero area are OK)
    for( xI in 1:n_x ){
      for( cI in 1:n_c ){
        maxsdD_xz[xI,zI] = maxsdD_xz[xI,zI] + pow(varD_xcz[xI,cI,zI], 0.5);
      }
      phi_xz[xI,zI] = varD_xz[xI,zI] / pow( maxsdD_xz[xI,zI], 2);
      varB_xbar_z[zI] = varB_xbar_z[zI] + pow(a_xl[xI,1],2) * varD_xz[xI,zI] * propB_xz[xI,zI];
      phi_xbar_z[zI] = phi_xbar_z[zI] + phi_xz[xI,zI] * propB_xz[xI,zI];
    }
    # Spatial-buffering index
    for( cI in 1:n_c ){
      for( xI in 1:n_x ){
        maxsdB_cz[cI,zI] = maxsdB_cz[cI,zI] + a_xl[xI,1] * pow(varD_xcz[xI,cI,zI], 0.5);
      }
      phi_cz[cI,zI] = varB_cz[cI,zI] / pow( maxsdB_cz[cI,zI], 2);
      varB_cbar_z[zI] = varB_cbar_z[zI] + varB_cz[cI,zI] * propB_cz[cI,zI];
      phi_cbar_z[zI] = phi_cbar_z[zI] + phi_cz[cI,zI] * propB_cz[cI,zI];
    }
    # Spatial and species-buffering index
    for( cI in 1:n_c ){
      for( xI in 1:n_x ){
        maxsdB_z[zI] = maxsdB_z[zI] + a_xl[xI,1] * pow(varD_xcz[xI,cI,zI], 0.5);
      }
    }
    phi_z[zI] = varB_z[zI] / pow( maxsdB_z[zI], 2);
  }

  # Return stuff
  Return = list( "phi_z"=phi_z, "phi_xz"=phi_xz, "phi_cz"=phi_cz, "phi_xbar_z"=phi_xz, "phi_cbar_z"=phi_cz, "B_y"=B_y, "B_cy"=B_cy, "D_xy"=D_xy, "D_xcy"=D_xcy)
  return( Return )
}
