#' Check fit for VAST model
#'
#' \code{check_fit} checks bounds and throws an informative message if any look bad
#'
#' @param parameter_estimates output from \code{TMBhelper::fit_tmb}
#' @param quiet Boolean stating whether to print warnings to terminal
#' @return Return Error code (TRUE is good; FALSE is bad)

#' @export
check_fit = function( parameter_estimates, quiet=FALSE ){

  # Initialize code for good model
  Code = TRUE

  # Check for informative bounds
  On_bounds = ifelse( parameter_estimates$par<(parameter_estimates$diagnostics[,'Lower']+0.0001) | parameter_estimates$par>(parameter_estimates$diagnostics[,'Upper']-0.0001), TRUE, FALSE )
  if( any(On_bounds) ){
    Code = FALSE
    if(quiet==FALSE){
      message("\nCheck bounds for the following parameters:")
      print( parameter_estimates$diagnostics[which(On_bounds),] )
    }
  }

  # Check for stuff at zero
  At_zero = ifelse( abs(parameter_estimates$par) < 0.0001, TRUE, FALSE )
  if( any(At_zero) ){
    Code = FALSE
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
  if( any(Bad_gradient) ){
    Code = FALSE
    if(quiet==FALSE){
      message("\nThe following parameters has a bad final gradient:")
      print( parameter_estimates$diagnostics[which(Bad_gradient),] )
    }
  }

  # Return
  return( invisible(Code) )
}

