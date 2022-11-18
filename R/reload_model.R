#' Reload a VAST model
#'
#' \code{reload_model} allows a user to save a fitted model, reload it in a new
#'      R terminal, and then relink the DLLs so that it functions as expected.
#'
#' @inheritParams make_model
#' @param x Output from \code{\link{fit_model}}, potentially with DLLs not linked
#' @param check_gradient Whether to check the gradients of the reloaded model
#'
#' @return Output from \code{\link{fit_model}} with DLLs relinked
#'
#' @examples
#' \dontrun{
#' # Run model
#' fit = fit_model( ... )
#' saveRDS( object=fit, file="path_and_name.rds" )
#'
#' # Reload and relink
#' fit_new = readRDS( file="path_and_name.rds" )
#' fit_new = reload_model( x = fit_new )
#' }
#'
#' @export
reload_model <-
function( x,
          check_gradient = TRUE,
          CompileDir = system.file("executables",package = "VAST"),
          Version = x$settings$Version,
          framework = x$input_args$model_args_input$framework,
          Obj = x$tmb_list$Obj ){

  # Load old one
  if( is.null(framework) ){
    Version_framework = Version
  }else{
    Version_framework = paste0( Version, "_", framework )
  }
  origwd = getwd()
  on.exit( setwd(origwd), add=TRUE )
  setwd(CompileDir)
  dyn.load( TMB::dynlib(Version_framework) ) # random=Random,

  # Retape
  Obj$retape()

  # Ensure that last.par and last.par.best are right
  Obj$fn(x$parameter_estimates$par)

  # Check gradient
  if( check_gradient==TRUE ){
    Gr = Obj$gr(x$parameter_estimates$par)
    if(max(abs(Gr))>1){
      warning("Maximum absolute gradient of ", signif(max(abs(Gr)),3), ": does not seem converged")
    }else if(max(abs(Gr))>0.01){
      warning("Maximum absolute gradient of ", signif(max(abs(Gr)),3), ": might not be converged")
    }else{
      message("Maximum absolute gradient of ", signif(max(abs(Gr)),3), ": No evidence of non-convergence")
    }
  }

  return(x)
}

