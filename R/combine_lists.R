
#' Combine lists
#'
#' \code{combine_lists} combines two lists with precident given to one over the other
#'
#' @param default default values for combined list
#' @param input preferred values for combined list
#'
#' @return combined list.
#'
#' @export
combine_lists = function( default, input, args_to_use=NULL, use_partial_matching=FALSE ){
  output = default
  for( i in seq_along(input) ){
    if( names(input)[i] %in% names(default) ){
      output[[names(input)[i]]] = input[[i]]
    }else{
      output = c( output, input[i] )
    }
  }
  if( !is.null(args_to_use) ){
    # Exact matching
    if( use_partial_matching==FALSE ){
      output = output[ intersect(names(output),args_to_use) ]
    }
    # Partial matching
    if( use_partial_matching==TRUE ){
      stop( "`use_partial_matching=TRUE` is not implemented" )
    }
  }
  return( output )
}

