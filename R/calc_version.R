#' Determine latest version of VAST
#'
#' @param version The default is \code{NULL}, which will cause the function
#' to look for the latest version of the \code{.cpp} file to compile.
#' If a version is supplied, then \code{R} will look for that version name
#' within the appropriate folder the disk.
#'
#' @return A full file path as a character value supplying the location of
#' the latest, or supplied, version of the VAST code to compile.
#'
#' @author Kelli Faye Johnson
#' @export
#'
calc_version <- function(version = NULL) {
  thedir <- system.file("executables", package = "VAST")

  if (!is.null(version)) {
    thefile <- dir(thedir, pattern = version, ignore.case = TRUE,
      full.names = TRUE)
    if (length(thefile) == 0) stop("The file ", version, " was not",
      " found in the dir, ", thedir, ".")
  } else {
    thefile <- dir(thedir, full.names = TRUE, pattern = "\\.cpp")
    if (length(thefile) == 0) stop("cpp files were not found in the ",
      "dir, ", thedir, ".")
    thefile <- tail(thefile, 1)
  }

  return(thefile)
}
