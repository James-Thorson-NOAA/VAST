

#' Format density covariate matrix
#'
#' \code{make_covariates} uses a formula interface to generate covariates
#'
#' This function generates 3D arrays \code{Cov_gtp}, \code{Cov_ip}, and  \code{Cov_stp} as
#' required by \code{VAST::make_data} to incorporate density covariates.
#' The user must supply a data frame \code{covariate_data} of covariate values, with columns named \code{Lat}, \code{Lon}, and \code{Year},
#' as well as values for all covariates as additional named columns.
#' This data frame is then used as a "look-up table", and is matched against variables listed in \code{formula}.
#'
#' Covariates then affect the linear predictor via coefficients that are estimated separately for every category.
#' Therefore, in a multivariate model the formula specifies the structure of covariates that is applied separately
#' to each category.
#'
#' Specifically,
#' \enumerate{
#' \item for every observation \code{i} at location \code{Lat_i[i]} and \code{Lon_i[i]} in year \code{t_i[t]}, the nearest
#'       Lat-Lon observation in that year is identified in \code{covariate_data}, and covariate
#'       values in that row of \code{covariate_data} are assigned to observation \code{i}.
#'       Covariate values are then used to predict densities at each sample location.
#' \item Similarly, for every extrapolation-grid cell \code{g} at location \code{spatial_list$latlon_g[g,]} in each year,
#'       the nearest row of \code{covariate_data} in that year
#'       is used to assign covariate values.
#'       Covariate values are then used to predict densities at each extrapolation-grid cell.
#' \item Finally, for every mesh-location \code{s} at location \code{spatial_list$latlon_s[s,]} in each year,
#'       the nearest row of \code{covariate_data} in that year
#'       is used to assign covariate values.
#'       Covariate values could be used in the future to predict habitat-specific variance-covariance
#'       for spatial variables (although this feature has not been built yet).
#' }
#' \code{make_covariates} then formats these covariate values appropriately and returns them.
#'
#' If all covariates as "static" (not changing among years),
#' then set \code{Year = NA} to cause values to be duplicated internally for all values of Year.
#' If using a mix of static and dynamic covariates,
#' then duplicate rows for static covariates for every value of Year
#'
#' @inheritParams stats::model.matrix
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. Similar specification to \code{\link[stats]{lm}}
#' @param covariate_data data frame of covariate values with columns \code{Lat}, \code{Lon}, and \code{Year}, and other columns matching names in \code{formula}; \code{Year=NA} can be used for covariates that do not change among years (e.g., depth)
#'
#' @return Tagged list of useful output
#' \describe{
#'   \item{\code{Cov_gtp}}{3-dimensional array for use in \code{VAST::make_data}}
#'   \item{\code{Cov_ip}}{2-dimensional array for use in \code{VAST::make_data}}
#'   \item{\code{Cov_stp}}{3-dimensional array for use in \code{VAST::make_data}}
#' }

#' @export
make_covariates <-
function( formula,
          covariate_data,
          Year_i,
          spatial_list,
          contrasts.arg = NULL ){

  # Check for bad entries
  if( !is.data.frame(covariate_data) ) stop("Please ensure that `covariate_data` is a data frame")
  if( !all(c("Lat","Lon","Year") %in% names(covariate_data)) ){
    stop( "`data` in `make_covariates(.)` must include columns `Lat`, `Lon`, and `Year`" )
  }

  # set of years needed
  Year_Set = min(Year_i):max(Year_i)

  # Make `covariate_df` by expanding for rows with Year=NA
  covariate_df = covariate_data[ which(!is.na(covariate_data[,'Year'])), ]
  for( tI in seq_along(Year_Set) ){
    newrows = covariate_data[ which(is.na(covariate_data[,'Year'])), ]
    newrows[,"Year"] = rep( Year_Set[tI], nrow(newrows) )
    covariate_df = rbind( covariate_df, newrows )
  }

  # Make model.matrix
    # To ensure identifiability given betas (intercepts), add intercept to formula
    # and then remove that term from model.matrix. This will not fix identifiability
    # issues arising when both conditions are met:
    # factor(Year) has an interaction with another factor, and
    # betas vary among years (are not constant)
  Model_matrix = model.matrix( update.formula(formula, ~.+1), data=covariate_df, contrasts.arg=contrasts.arg )
  Columns_to_keep = which( attr(Model_matrix,"assign") != 0 )
  coefficient_names = attr(Model_matrix,"dimnames")[[2]][Columns_to_keep]
  X = Model_matrix[,Columns_to_keep,drop=FALSE]
  dimnames(X) = list(NULL, coefficient_names)

  # transform data inputs
  # Used to define covariates at samples
  sample_i = data.frame( "Year"=Year_i, "Lat"=spatial_list$latlon_i[,'Lat'], "Lon"=spatial_list$latlon_i[,'Lon'] )

  # extract latitude and longitude for extrapolation grid
  # Used to define covariates in extrapolation grid
  latlon_g = spatial_list$latlon_g

  #
  latlon_s = spatial_list$latlon_s

  # Create data frame of necessary size
  X_gtp = array( NA, dim=c(nrow(latlon_g),length(Year_Set),ncol(X)), dimnames=list(NULL,Year_Set,colnames(X)) )
  X_ip = array( NA, dim=c(nrow(sample_i),ncol(X)), dimnames=list(NULL,colnames(X)) )
  X_stp = array( NA, dim=c(nrow(latlon_s),length(Year_Set),ncol(X)), dimnames=list(NULL,Year_Set,colnames(X)) )

  # Loop through data and extrapolation-grid
  for( tI in seq_along(Year_Set) ){

    # Subset to same year
    tmp_covariate_df = covariate_df[ which(Year_Set[tI]==covariate_df[,'Year']), , drop=FALSE]
    tmp_X = X[ which(Year_Set[tI]==covariate_df[,'Year']), , drop=FALSE]
    if( nrow(tmp_covariate_df)==0 ){
      stop("Year ", Year_Set[tI], " not found in `covariate_data` please specify covariate values for all years" )
    }

    # Fill in values in X_ip
    Which = which(Year_Set[tI]==sample_i[,'Year'])
    # Do nearest neighbors to define covariates for observations, skipping years without observations
    if( length(Which) > 0 ){
      NN = RANN::nn2( data=tmp_covariate_df[,c("Lat","Lon")], query=sample_i[Which,c("Lat","Lon")], k=1 )
      # Fill in values
      X_ip[Which, ] = tmp_X[ NN$nn.idx[,1], , drop=FALSE]
    }

    # Do nearest neighbors to define covariates for extrapolation grid, including years without observations
    NN = RANN::nn2( data=tmp_covariate_df[,c("Lat","Lon")], query=latlon_g[,c("Lat","Lon")], k=1 )
    # Add rows
    X_gtp[,tI,] = tmp_X[ NN$nn.idx[,1], , drop=FALSE ]

    # Do nearest neighbors to define covariates for extrapolation grid, including years without observations
    NN = RANN::nn2( data=tmp_covariate_df[,c("Lat","Lon")], query=latlon_s[,c("Lat","Lon")], k=1 )
    # Add rows
    X_stp[,tI,] = tmp_X[ NN$nn.idx[,1], , drop=FALSE ]
  }

  # Check for obvious problems
  if( any(is.na(X_ip)) ) stop("Problem with `X_ip` in `make_covariates(.)")
  if( any(is.na(X_gtp)) ) stop("Problem with `X_gtp` in `make_covariates(.)")
  if( any(is.na(X_stp)) ) stop("Problem with `X_stp` in `make_covariates(.)")

  # warnings
  if( any(apply(X_gtp, MARGIN=2:3, FUN=sd)>10 | apply(X_ip, MARGIN=2, FUN=sd)>10) ){
    warning("The package author recommends that you rescale covariates in `covariate_data` to have mean 0 and standard deviation 1.0")
  }

  # return stuff
  Return = list( "X_gtp" = X_gtp,
           "X_ip" = X_ip,
           "X_stp" = X_stp,
           "coefficient_names" = coefficient_names )
  return( Return )
}

