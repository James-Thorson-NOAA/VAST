
#' @title
#' Remove units from object
#'
#' @description
#' \code{strip_units} facilitates logical comparisons and other generic operators by stripping units from an object but otherwise leaving it unchanged
#'
#' @export
strip_units = function(x){
  if("units" %in% class(x)) units(x) = NULL
  return(x)
}
