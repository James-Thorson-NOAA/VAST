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
          CompileDir = system.file("executables",package = "VAST") ){

  # Load old one
  origwd = getwd()
  on.exit( setwd(origwd), add=TRUE )
  setwd(CompileDir)
  Version = x$settings$Version
  TMB::compile( paste0(Version,".cpp"), flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign" )
  dyn.load( TMB::dynlib(Version) ) # random=Random,

  # Retape
  x$tmb_list$Obj$retape()

  # Ensure that last.par and last.par.best are right
  x$tmb_list$Obj$fn(x$par$par)

  # Check gradient
  if( check_gradient==TRUE ){
    Gr = x$tmb_list$Obj$gr(x$par$par)
    if(max(abs(Gr))>1){
      stop("Maximum absolute gradient of ", max(abs(Gr)), ": does not seem converged")
    }else if(max(abs(Gr))>0.01){
      warning("Maximum absolute gradient of ", max(abs(Gr)), ": might not be converged")
    }else{
      print("Maximum absolute gradient of ", max(abs(Gr)), ": might be converged")
    }
  }

  return(x)
}

