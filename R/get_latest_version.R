#' Determine latest version of VAST
#'
#' @param version The default is \code{NULL}, which will cause
#'   the function to look for the latest version of the
#'   \code{.cpp} file to compile.  If a version is supplied, then
#'   \code{R} will look for that version name within the
#'   appropriate folder the disk.
#' @param path Optional path for running tests locally and with
#'   continuous integration since the default looks in the wrong
#'   place (installed library).
#'
#' @return A full file path as a character value supplying the location of
#' the latest, or supplied, version of the VAST code to compile.
#'
#' @author Kelli Faye Johnson
#'
#' @export
get_latest_version <- function(version = NULL, package = "VAST", path=NULL) {

  # Add file type to avoid compiled renames
  if (!is.null(version)) {
    version_extension = paste0( version, ".cpp" )
  }else{
    version_extension = NULL
  }

  # Determine location of files on machine
  if(!is.null(path)){
    thedir <- path
    if(!dir.exists(path)){
      warning("Manually specified path to executables is not valid, trying system files")
      thedir <- system.file("executables", package = package)
    }
  } else {
    thedir <- system.file("executables", package = package)
  }

  if(thedir=="")
    stop("Could not find an executables folder for ", package,
         "\n Something is likely wrong with the installation. Try manually specifiying\n",
         "the version number, or try reinstalling ", package)

  # Determine list of available files
  if (!is.null(version_extension)) {
    thefile <- dir(thedir, pattern = version_extension, ignore.case = TRUE, full.names = FALSE)
    if( length(thefile)==0 ){
      stop("The file ", version, " was not found in the dir, ", thedir, ".")
    }
  }else{
    thefile <- dir(thedir, full.names = FALSE, pattern = "\\.cpp")
    if( length(thefile)==0 ){
      stop("cpp files were not found in the dir, ", thedir, ".")
    }
    # remove _TMBad and _CppAD
    thefile = gsub( thefile, pattern="_TMBad", replacement="", ignore.case=TRUE )
    thefile = gsub( thefile, pattern="_CppAD", replacement="", ignore.case=TRUE )
    thefile = unique(thefile)
  }

  # Determine which is latest version
  for(i in 1:length(thefile)){
    if(i==1) semantic_version = convert_version_name(thefile[1])
    if(i>=2) semantic_version = c( semantic_version, convert_version_name(thefile[i]) )
  }
  #semantic_version = sapply( thefile, FUN=convert_version_name )
  if( max(semantic_version) == numeric_version("0.0.0") ) stop("Problem with `get_latest_version`")
  thefile = thefile[which(semantic_version==max(semantic_version))]

  # Remove .cpp from end and return
  thefile <- strsplit( thefile, ".", fixed=TRUE )[[1]][1]
  return(thefile)
}
