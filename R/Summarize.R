
#' Summarize and explain results
#'
#' \code{Summarize} is intended to help document results, by providing explanation for outputs from \code{obj$report()}
#'
#' @param obj the TMB object after it has been previously optimized

#' @return Tagged list of useful output
#' \describe{
#'   \item{parhat}{a list of fixed-effect estimates, with same format as starting values}
#'   \item{D_xct}{an array (dimensions: knots; categories, e.g., species/age;  year) of predicted density (doesn't include bias-correction)}
#' }
#'
#' @details
#' Please post on the GitHub issue tracker if you'd like for another specific output to be added and documented here

#' @export
Summarize = function( obj ){

  # extract elements
  report = obj$report()
  parhat = obj$env$parList()

  # Slots worth describing
  D_xct = Report$D_xct

  # Bundle and return
  Return = list( "parhat"=parhat, "D_xct"=D_xct )
  return( Return )
}

