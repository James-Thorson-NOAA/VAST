#' Check fit for VAST model
#'
#' \code{check_fit} checks bounds and throws an informative message if any look bad
#'
#' If `check_fit` identifies an issue in estimated parameters, then the model structure should typically be changed.
#' Recommended model changes differ somewhat for univariate and multivariate models as explained below.
#'
#' For univariate models:
#' \itemize{
#' \item If \code{ln_H_input} are approaching extreme values (i.e., > 5 or < -5), then turn consider turning off anisotropy, \code{make_data(..., Aniso=FALSE)}
#' \item If \code{L_beta1_z} is approaching zero (i.e., +/- 0.001), then turn off random-effects for temporal variation in the first intercept \code{RhoConfig["Beta1"]=3}
#' \item If \code{L_beta2_z} is approaching zero (i.e., +/- 0.001), then turn off random-effects for temporal variation in the first intercept \code{RhoConfig["Beta2"]=3}
#' \item If \code{L_omega1_z} is approaching zero (i.e., +/- 0.001), then turn off spatial effects for the 1st linear predictor `FieldConfig["Omega1"]=0`
#' \item If \code{L_omega2_z} is approaching zero (i.e., +/- 0.001), then turn off spatial effects for the 1st second predictor `FieldConfig["Omega2"]=0`
#' \item If \code{L_epsilon1_z} is approaching zero (i.e., +/- 0.001), then turn off spatio-temporal effects for the 1st linear predictor `FieldConfig["Epsilon1"]=0`
#' \item If \code{L_epsilon2_z} is approaching zero (i.e., +/- 0.001), then turn off spatio-temporal effects for the 1st second predictor `FieldConfig["Epsilon2"]=0`
#' \item If \code{Beta_rho1_f} is approaching one (i.e., > 0.999), then turn consider reducing to a random-walk structure for the intercept of the 1st linear predictor `RhoConfig["Beta1"]=2`
#' \item If \code{Beta_rho2_f} is approaching one (i.e., > 0.999), then turn consider reducing to a random-walk structure for the intercept of the 2st linear predictor `RhoConfig["Beta2"]=2`
#' \item If \code{Epsilon_rho1_f} is approaching one (i.e., > 0.999), then turn consider reducing to a random-walk structure for spatio-temporal variation of the 1st linear predictor `RhoConfig["Epsilon1"]=2`
#' \item If \code{Epsilon_rho2_f} is approaching one (i.e., > 0.999), then turn consider reducing to a random-walk structure for spatio-temporal variation of the 2nd linear predictor `RhoConfig["Epsilon2"]=2`
#' }
#'
#' For multivariate models, these same principles apply, but there are more options to simplify model structure.
#' For example, if any `L_beta1_z` is approaching zero (i.e., +/- 0.001), then consider using `fit_model(...,Map=[custom-map])` to turn off individual parameters;
#' or if using a factor model then reduce the number of factors by decreasing \code{FieldConfig["Beta1"]}
#'
#' @param parameter_estimates output from \code{\link{fit_tmb}}
#' @param check_gradients Boolean stating whether to check bounds as well as other issues
#' @param quiet Boolean stating whether to print warnings to terminal
#' @return Did an automated check find an obvious problem code (TRUE is bad; FALSE is good)
#'
#' @export
#' @md
# Using https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html for guidance on markdown-enabled documentation
check_fit <-
function( parameter_estimates,
          check_gradients=FALSE,
          quiet=FALSE ){

  # Initialize code for good model
  problem_found = FALSE

  # Check for informative bounds
  On_bounds = ifelse( parameter_estimates$par<(parameter_estimates$diagnostics[,'Lower']+0.0001) | parameter_estimates$par>(parameter_estimates$diagnostics[,'Upper']-0.0001), TRUE, FALSE )
  if( any(On_bounds) ){
    problem_found = TRUE
    if(quiet==FALSE){
      message("\nCheck bounds for the following parameters:")
      print( parameter_estimates$diagnostics[which(On_bounds),] )
    }
  }

  # Check for stuff at zero
  At_zero = ifelse( abs(parameter_estimates$par) < 0.0001, TRUE, FALSE )
  if( any(At_zero) ){
    problem_found = TRUE
    if(quiet==FALSE){
      message("\nThe following parameters appear to be approaching zero:")
      print( parameter_estimates$diagnostics[which(At_zero),] )
      if( length(grep( "L_", names(parameter_estimates$par[which(At_zero)]))) ){
        message("Please turn off factor-model variance parameters `L_` that are approaching zero and re-run the model")
      }
    }
  }

  # Check bad gradients
  Bad_gradient = ifelse( abs(parameter_estimates$diagnostics[,'final_gradient']) > 0.001, TRUE, FALSE )
  if( check_gradients==TRUE && any(Bad_gradient) ){
    problem_found = TRUE
    if(quiet==FALSE){
      message("\nThe following parameters has a bad final gradient:")
      print( parameter_estimates$diagnostics[which(Bad_gradient),] )
    }
  }

  # Return
  return( invisible(problem_found) )
}

