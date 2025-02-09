
#' @title
#' Simulate new dynamics and sampling data
#'
#' @description
#' \code{simulate_data} conducts a parametric bootstrap to simulate new data and potentially simulate new population dynamics and associated variables
#'
#' Simulate new data given various potential procedures to propagate uncertainty about parameters.
#'
#' Using \code{sample_fixed=TRUE} (the default) in \code{\link{sample_variable}} is similar to using \code{type=3} in \code{\link{simulate_data}}, while
#'       using \code{sample_fixed=TRUE} in \code{\link{sample_variable}} is similar to using \code{type=4} in \code{\link{simulate_data}}.
#'       Sampling fixed effects will sometimes cause numerical under- or overflow (i.e., output values of \code{NA}) in cases when
#'       variance parameters are estimated imprecisely.  In these cases, the multivariate normal approximation being used is a poor
#'       representation of the tail probabilities, and results in some samples with implausibly high (or negative) variances,
#'       such that the associated random effects then have implausibly high magnitude.
#'
#' @param fit output from \code{\link{fit_model}}
#' @param type integer stating what type of simulation to use from the following options:
#' \itemize{
#' \item \code{type=1} is a "measurement error" or "conditional" simulator that simulates new data conditional upon estimated fixed and random effects.
#' \item \code{type=2} is an "unconditional" simulator that simulates new random effects conditional upon fixed effects
#' (but not otherwise conditioning upon original data), and new data conditional upon both.
#' \item \code{type=3} simulates new fixed and random effects from the joint precision matrix (i.e., conditioning upon the original data), and new data conditional upon these values.
#' \item \code{type=4} simulates new random effects from the internal Hessian matrix evaluated at the MLE (i.e., conditional on fixed effects estimates and the original data),
#' and new data conditional upon these values.
#' }
#' @param random_seed integer passed to \code{set.seed}, where the default value \code{random_seed=NULL} resets the random-number seed.
#'

#' @return Report object containing new data and population variables including
#' \describe{
#'   \item{b_i}{New simulated data}
#'   \item{D_gcy}{Density for each grid cell g, category c, and year y}
#'   \item{Index_cyl}{Index of abundance for each category c, year y, and stratum l}
#' }

#' @export
simulate_data <-
function( fit,
          type = 1,
          random_seed = NULL ){

  # Sample from GMRF using sparse precision
  rmvnorm_prec <- function(mu, prec, n.sims, random_seed ) {
    set.seed( random_seed )
    z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L = Matrix::Cholesky(prec, super=TRUE)
    z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z = as.matrix(z)
    return(mu + z)
  }

  # Check for loaded VAST
  # Modified from TMB:::getUserDLL
  dlls <- getLoadedDLLs()
  isTMBdll <- function(dll) !is(try(getNativeSymbolInfo("MakeADFunObject",dll), TRUE), "try-error")
  TMBdll <- sapply(dlls, isTMBdll)
  if( sum(TMBdll)==0 ){
    stop("VAST is not linked as a DLL, so `simulate_data` will not work.
    Please re-load model using `reload_model` to use `simulate_data`")
  }else if(sum(TMBdll)>=2){
    warning("VAST is linked to multiple DLLs. Please consider using dyn.unload() to unload
    earlier VAST runs to avoid potentially ambiguous behavior when running `simulate_data`")
  }

  # Extract stuff
  Obj = fit$tmb_list$Obj
  simulate_random_effects_orig = Obj$env$data$Options_list$Options['simulate_random_effects']

  # Revert settings when done
  revert_settings = function(simulate_random_effects){Obj$env$data$Options_list$Options['simulate_random_effects'] = simulate_random_effects}
  on.exit( revert_settings(simulate_random_effects_orig) )

  # Simulate conditional upon fixed and random effect estimates
  if( type==1 ){
    # Change and revert settings
    Obj$env$data$Options_list$Options['simulate_random_effects'] = FALSE
    set.seed(random_seed)
    Return = Obj$simulate( complete=TRUE )
  }

  # Simulate new random effects and data
  if( type==2 ){
    Obj$env$data$Options_list$Options['simulate_random_effects'] = TRUE
    set.seed(random_seed)
    Return = Obj$simulate( complete=TRUE )
  }

  # Simulate from predictive distribution of fixed AND random effects, and then new data
  # Could instead explore fit$tmb_list$Obj$env$MC(.) for sampling-importance-resampling approach
  if( type==3 ){
    # Informative error messages
    if( !("jointPrecision" %in% names(fit$parameter_estimates$SD)) ){
      stop("jointPrecision not present in fit$parameter_estimates$SD; please re-run with `getJointPrecision=TRUE`")
    }

    # Sample from joint distribution
    newpar = rmvnorm_prec( mu=Obj$env$last.par.best, prec=fit$parameter_estimates$SD$jointPrecision, n.sims=1, random_seed=random_seed )[,1]

    # Simulate
    Obj$env$data$Options_list$Options['simulate_random_effects'] = FALSE
    Return = Obj$simulate( par=newpar, complete=TRUE )
  }

  # Simulate from predictive distribution of random effects and NOT fixed effects, and then new data
  if( type==4 ){
    set.seed( random_seed )
    warning( "Type-4 residuals are still under development, please use with care and note that they may change at any point.")
    newpar = Obj$env$last.par.best
    #Hess = Obj$env$spHess( par=Obj$env$last.par.best, random=TRUE )
    #newrandom = rmvnorm_prec( mu=rep(0,nrow(Hess)), prec=Hess, n.sims=1, random_seed=random_seed )[,1]
    #newpar[Obj$env$random] = newrandom
    MC = Obj$env$MC( keep=TRUE, n=1, antithetic=FALSE )
    newpar[Obj$env$random] = attr(MC, "samples")

    # Simulate
    Obj$env$data$Options_list$Options['simulate_random_effects'] = FALSE
    Return = Obj$simulate( par=newpar, complete=TRUE )
  }

  # Return
  return( Return )
}


