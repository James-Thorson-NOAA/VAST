#' Check fit for VAST model
#'
#' \code{check_fit} checks bounds and throws an informative message if any look bad
#'
#' @param parameter_estimates output from \code{TMBhelper::fit_tmb}
#' @param check_bounds Boolean stating whether to check bounds as well as other issues
#' @param quiet Boolean stating whether to print warnings to terminal
#' @return Did an automated check find an obvious problem code (TRUE is bad; FALSE is good)

#' @export
check_fit = function( parameter_estimates, check_bounds=FALSE, quiet=FALSE ){

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
  if( check_bounds==TRUE && any(Bad_gradient) ){
    problem_found = TRUE
    if(quiet==FALSE){
      message("\nThe following parameters has a bad final gradient:")
      print( parameter_estimates$diagnostics[which(Bad_gradient),] )
    }
  }

  # Return
  return( invisible(problem_found) )
}

