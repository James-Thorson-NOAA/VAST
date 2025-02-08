#' Convert string to semantic version number
#'
#' @param char
#'
#' @return A semantic version numvber
#'
#'
#' @export
convert_version_name <- function( char, split="_" ) {

  # Split into pieces
  convert = function(part) paste0(na.omit(as.numeric(strsplit(part,"")[[1]])),collapse="")

  # Reduce to numbers within each piece
  Split = strsplit( char, split=split )[[1]]
  Codes = suppressWarnings(sapply( Split, FUN=convert))
  Codes = paste0( na.omit(ifelse(Codes=="",NA,Codes)), collapse="." )

  # Convert to semantic number and return
  if( any(Codes != "") ){
    semantic_version = numeric_version(Codes)
  }else{
    semantic_version = numeric_version("0.0.0")
  }
  return(semantic_version)
}
