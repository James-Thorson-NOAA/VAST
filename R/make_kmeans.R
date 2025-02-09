
#' Calculate location for knots approximating spatial variation
#'
#' \code{make_kmeans} determines the location for a set of knots for approximating spatial variation
#'
#' @param n_x the number of knots to select
#' @param loc_orig a matrix with two columns where each row gives the 2-dimensional coordinates to be approximated
#' @param nstart the number of times that the k-means algorithm is run while searching for the best solution (default=100)
#' @param randomseed a random number seed
#' @param iter.max the number of iterations used per k-means algorithm (default=1000)
#' @param DirPath a directory where the function looks for a previously-saved output (default is working directory)
#' @param Save_Results a boolean stating whether to save the output (Default=TRUE)
#' @param kmeans_purpose a character representing whether the call is to calculate "extrapolation" or "spatial" information
#' @param backwards_compatible_kmeans a boolean stating how to deal with changes in the kmeans algorithm implemented in R version 3.6.0,
#'        where \code{backwards_compatible_kmeans==TRUE} modifies the default algorithm to maintain backwards compatibility, and
#'        where \code{backwards_compatible_kmeans==FALSE} breaks backwards compatibility between R versions prior to and after R 3.6.0.

#' @return Tagged list containing outputs
#' \describe{
#'   \item{centers}{a matrix with 2 columns and n_x rows}
#'   \item{cluster}{A vector with length \code{nrow(loc_orig)} specifying which row of \code{centers} corresponds to each row of loc_orig}
#' }

#' @export
make_kmeans <-
function( n_x,
          loc_orig,
          nstart = 100,
          randomseed = 1,
          iter.max = 1000,
          DirPath = getwd(),
          Save_Results = TRUE,
          kmeans_purpose = "spatial",
          backwards_compatible_kmeans = FALSE ){

  # get old seed
  oldseed = ceiling(runif(1,min=1,max=1e6))
  # fix new seed
  if( !is.null(randomseed) ) set.seed( round(randomseed) )
  old.options <- options()
  options( "warn" = -1 )
  on.exit( options(old.options) )

  # warnings
  if( paste0("Kmeans-",n_x,".RData") %in% list.files(DirPath) ){
    warning("Found `Kmeans-",n_x,".RData` in directory; this might indicate a desire to import previously constructed k-means under deprecated naming conventions")
  }

  if(kmeans_purpose=="spatial"){
    tmpfile <- paste0("Kmeans_knots-",n_x,".RData")
  } else if(kmeans_purpose=="extrapolation"){
    ## n_x is really max_cells in this case
    tmpfile <- paste0("Kmeans_extrapolation-",n_x,".RData")
  } else {
    stop("Invalid kmeans_purpose for make_kmeans:", kmeans_purpose)
  }

  # Backwards compatibility
  if( backwards_compatible_kmeans==TRUE ){
    if( identical(formalArgs(RNGkind), c("kind","normal.kind","sample.kind")) ){
      RNGkind_orig = RNGkind()
      on.exit( RNGkind(kind=RNGkind_orig[1], normal.kind=RNGkind_orig[2], sample.kind=RNGkind_orig[3]), add=TRUE )
      RNGkind( sample.kind="Rounding" )
    }else if( !identical(formalArgs(RNGkind), c("kind","normal.kind")) ){
      stop("Assumptions about `RNGkind` are not met within `make_kmeans`; please report problem to package developers")
    }
  }

  # Calculate knots for SPDE mesh
  if( length(unique(paste(loc_orig[,1],loc_orig[,2],sep="_")))<=n_x ){
    # If number of knots is less than number of sample locations
    Kmeans = NULL
    Kmeans[["centers"]] = unique( loc_orig )
    Kmeans[["cluster"]] = RANN::nn2( data=Kmeans[["centers"]], query=loc_orig, k=1)$nn.idx[,1]
    message( "n_x greater than or equal to n_unique so no calculation necessary" )
  }else{
    if(tmpfile  %in% list.files(DirPath) ){
      # If previously saved knots are available
      load( file=file.path(DirPath, tmpfile))
      message( "Loaded from ", file.path(DirPath, tmpfile) )
    }else{
      # Multiple runs to find optimal knots
      message("Using ", nstart, " iterations to find optimal ",
              ifelse(kmeans_purpose=="extrapolation", "extrapolation grid", "spatial knot"),
              " placement because no saved file found...")
      Kmeans = list( "tot.withinss"=Inf )
      tries <- 1
      for(i in 1:nstart){
        Tmp = stats::kmeans( x=loc_orig, centers=n_x, iter.max=iter.max, nstart=1, trace=0)
        if(i==1) message("Iter=1: Current=", round(Tmp$tot.withinss,0))
        if(i!=1 & i!=nstart & i %% 10 ==0)
          message( 'Iter=',i,': Current=',round(Kmeans$tot.withinss,0),' Proposed=',round(Tmp$tot.withinss,0) )#,' Time=',round(Time,4)) )
        if( Tmp$tot.withinss < Kmeans$tot.withinss ){
          Kmeans = Tmp
          tries <- i # which iteration was the optimal
        }
      }
      message("Iter=", nstart, ': Final=', round(Kmeans$tot.withinss,0), " after ", tries, " iterations")
      if(Save_Results==TRUE){
        save( Kmeans, file=file.path(DirPath, tmpfile))
        message( "Results saved to ", file.path(DirPath, tmpfile), "\n for subsequent runs by default (delete it to override)")
      }
    }
  }

  # fix to old seed
  if( !is.null(randomseed) ) set.seed( oldseed )

  # Return stuff
  Return = list("centers"=Kmeans[["centers"]], "cluster"=Kmeans[["cluster"]] )
  return( Return )
}
