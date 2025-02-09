
#' @title
#' Sample from predictive distribution of a variable
#'
#' @description
#' \code{sample_variable} samples from the joint distribution of random and fixed effects to approximate the predictive distribution for a variable
#'
#' Using \code{sample_fixed=TRUE} (the default) in \code{\link{sample_variable}} is similar to using \code{type=3} in \code{\link{simulate_data}}, while
#'       using \code{sample_fixed=TRUE} in \code{\link{sample_variable}} is similar to using \code{type=4} in \code{\link{simulate_data}}.
#'       Sampling fixed effects will sometimes cause numerical under- or overflow (i.e., output values of \code{NA}) in cases when
#'       variance parameters are estimated imprecisely.  In these cases, the multivariate normal approximation being used is a poor
#'       representation of the tail probabilities, and results in some samples with implausibly high (or negative) variances,
#'       such that the associated random effects then have implausibly high magnitude.
#'
#' @param Sdreport TMB output from `\code{TMB::sdreport(Obj)}`
#' @param Obj Fitted TMB object from package `VAST`, i.e., output from `\code{fit_model(...)$tmb_list$Obj}`
#' @param variable_name name of variable available in report using \code{Obj$report()} or parameters using \code{Obj$env$parList()}
#' @param n_samples number of samples from the joint predictive distribution for fixed and random effects.  Default is 100, which is slow.
#' @param seed integer used to set random-number seed when sampling variables, as passed to \code{set.seed(.)}
#' @param sample_fixed whether to sample fixed and random effects, \code{sample_fixed=TRUE} as by default, or just sample random effects, \code{sample_fixed=FALSE}
#'
#' @examples
#' \dontrun{
#' # Run model using selected inputs, but also with getJointPrecision=TRUE
#' fit = fit_model( ...,
#'     getJointPrecision=TRUE )
#'
#' # Run sample_variable
#' sample = sample_variable( Sdreport=fit$parameter_estimates$SD,
#'     Obj=fit$tmb_list$Obj, variable_name="D_gct" )
#' }
#'
#' @export
sample_variable <-
function( Sdreport,
          Obj,
          variable_name,
          n_samples = 100,
          sample_fixed = TRUE,
          seed = 123456 ){

  # Informative error messages
  if( !("jointPrecision" %in% names(Sdreport)) ){
    stop("jointPrecision not present in Sdreport; please re-run with `getJointPrecision=TRUE`")
  }
  # Combine Report and ParHat, and check for issues
  Report = Obj$report()
  ParHat = Obj$env$parList()
  Intersect = intersect(names(Report), names(ParHat))
  if( isFALSE(all.equal(Report[Intersect],ParHat[Intersect])) ){
    stop("Duplicate entries in `Obj$report()` and `Obj$env$parList()` are not identical when calling `sample_variable`")
  }
  Output = c( Report, ParHat )
  # Check that variable_name is available
  if( !(variable_name %in% names(Output)) ){
    stop( variable_name, " not found in `Obj$report()` or `Obj$env$parList()`; please choose check your requested variable name from available list: ", paste(names(Output),collapse=", ") )
  }

  #### Local function
  # Sample from GMRF using sparse precision
  rmvnorm_prec <- function(mu, prec, n.sims, seed) {
    set.seed(seed)
    z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- Matrix::Cholesky(prec, super=TRUE)
    z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    return(mu + z)
  }

  # Sample from joint distribution
  if( sample_fixed==TRUE ){
    u_zr = rmvnorm_prec( mu=Obj$env$last.par.best, prec=Sdreport$jointPrecision, n.sims=n_samples, seed=seed)
    # apply( u_zr, MARGIN=2, FUN=function(vec){sum(abs(vec)==Inf)})
    # u_zr[-Obj$env$random,1]
  }else{
    u_zr = Obj$env$last.par.best %o% rep(1, n_samples)
    MC = Obj$env$MC( keep=TRUE, n=n_samples, antithetic=FALSE )
    u_zr[Obj$env$random,] = attr(MC, "samples")
  }

  # Extract variable for each sample
  message( "# Obtaining samples from predictive distribution for variable ", variable_name )
  for( rI in 1:n_samples ){
    if( rI%%max(1,floor(n_samples/10)) == 0 ){
      message( "  Finished sample ", rI, " of ",n_samples )
    }
    Report = Obj$report( par=u_zr[,rI] )
    ParHat = Obj$env$parList( x=u_zr[,rI][Obj$env$lfixed()], par=u_zr[,rI] )
    if( isFALSE(all.equal(Report[Intersect],ParHat[Intersect])) ){
      stop("Duplicate entries in `Obj$report()` and `Obj$env$parList()` are not identical when calling `sample_variable`")
    }
    Output = c( Report, ParHat )
    Var = Output[[variable_name]]
    if(is.vector(Var)) Var = as.array(Var)
    if(rI==1) Var_zr = Var
    if(rI>=2){
      Var_zr = abind::abind( Var_zr, Var, along=length(dim(Var))+1 )
    }
  }

  # Return
  return( Var_zr )
}


