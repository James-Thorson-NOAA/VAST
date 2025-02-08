
#' Optimize a TMB model
#'
#' \code{fit_tmb} runs a TMB model and generates standard diagnostics
#'
#' @inheritParams stats::nlminb
#' @inheritParams TMB::sdreport
#' @param obj The compiled TMB object
#' @param startpar Starting values for fixed effects (default NULL uses \code{obj$par})
#' @param control A list of control parameters. For details see \code{\link[stats]{nlminb}}
#' @param getsd Boolean indicating whether to run standard error calculation; see \code{\link[TMB]{sdreport}} for details
#' @param bias.correct Boolean indicating whether to do epsilon bias-correction;
#'        see \code{\link[TMB]{sdreport}} and \code{fit_tmb} for details
#' @param bias.correct.control tagged list of options for epsilon bias-correction,
#'        where \code{vars_to_correct} is a character-vector of ADREPORT variables that should be bias-corrected
#' @param savedir directory to save results (if \code{savedir=NULL}, then results aren't saved)
#' @param loopnum number of times to re-start optimization (where \code{loopnum=3}
#'        sometimes achieves a lower final gradient than \code{loopnum=1})
#' @param newtonsteps Integer specifying the number of extra newton steps to take
#'        after optimization (alternative to \code{loopnum}).
#'        Each newtonstep requires calculating the Hessian matrix and is therefore slow.
#'        But for well-behaved models, each Newton step will typically
#'        decrease the maximum gradient of the loglikelihood with respect to each fixed effect,
#'        and therefore this option can be used to achieve an arbitrarily low final gradient
#'        given sufficient time for well-behaved models.  However, this option will also
#'        perform strangely or have unexpected consequences for poorly-behaved models, e.g.,
#'        when fixed effects are at upper or lower bounds.
#' @param n sample sizes (if \code{n!=Inf} then \code{n} is used to calculate BIC and AICc)
#' @param getHessian return Hessian for usage in later code
#' @param quiet Boolean whether to print additional messages results to terminal
#' @param start_time_elapsed how much time has elapsed prior to calling fit_tmb,
#'        for use, e.g., when calling \code{fit_tmb} multiple times in sequence,
#'        where \code{start_time_elapsed = opt_previous$time_for_run}
#' @param ... list of settings to pass to \code{\link[TMB]{sdreport}}
#'
#' @return the standard output from \code{\link[stats]{nlminb}}, except with additional diagnostics and timing info,
#'         and a new slot containing the output from \code{\link[TMB]{sdreport}}
#'
#' @examples
#' fit_tmb( Obj ) # where Obj is a compiled TMB object
#'
#' @references For more details see \url{https://doi.org/10.1016/j.fishres.2015.11.016}
#' @export
fit_tmb <-
function( obj,
          fn = obj$fn,
          gr = obj$gr,
          startpar = NULL,
          lower = -Inf,
          upper = Inf,
          getsd = TRUE,
          control = list(eval.max = 1e4, iter.max = 1e4, trace = 1),
          bias.correct = FALSE,
          bias.correct.control = list(sd = FALSE, split = NULL, nsplit = NULL, vars_to_correct = NULL),
          savedir = NULL,
          loopnum = 2,
          newtonsteps = 0,
          n = Inf,
          getReportCovariance = FALSE,
          getJointPrecision = FALSE,
          getHessian = FALSE,
          quiet = FALSE,
          start_time_elapsed = as.difftime("0:0:0"),
          ... ){

  # Defaults
  if(is.null(startpar)) startpar = obj$par

  # Check for issues
  if( bias.correct==TRUE & is.null(obj$env$random) ){
    message( "No random effects detected in TMB model, so overriding user input to `fit_tmb` to instead specify `bias.correct=FALSE`")
    bias.correct = FALSE
  }
  if( getReportCovariance==FALSE ){
    message( "Note that `getReportCovariance=FALSE` causes an error in `TMB::sdreport` when no ADREPORTed variables are present")
  }

  # Check for issues
  List = list(...)
  #if( !is.null(List$jointJointPrecision) ){
  #  if( getJointPrecision==FALSE & getReportCovariance==TRUE ){
  #    stop("Some versions of TMB appear to throw an error when `getJointPrecision=TRUE` and `getReportCovariance=FALSE`")
  #  }
  #}

  # Local function -- combine two lists
  combine_lists = function( default, input ){
    output = default
    for( i in seq_along(input) ){
      if( names(input)[i] %in% names(default) ){
        output[[names(input)[i]]] = input[[i]]
      }else{
        output = c( output, input[i] )
      }
    }
    return( output )
  }

  # Replace defaults for `BS.control` with provided values (if any)
  BS.control = list(sd=FALSE, split=NULL, nsplit=NULL, vars_to_correct=NULL)
  BS.control = combine_lists( default=BS.control, input=bias.correct.control )

  # Replace defaults for `control` with provided values if any
  nlminb.control = list(eval.max=1e4, iter.max=1e4, trace=0)
  nlminb.control = combine_lists( default=nlminb.control, input=control )

  # Run first time
  start_time = Sys.time()
  parameter_estimates = nlminb( start=startpar, objective=fn, gradient=gr, control=nlminb.control, lower=lower, upper=upper )

  # Re-run to further decrease final gradient
  for( i in seq(2,loopnum,length=max(0,loopnum-1)) ){
    Temp = parameter_estimates[c('iterations','evaluations')]
    parameter_estimates = nlminb( start=parameter_estimates$par, objective=fn, gradient=gr, control=nlminb.control, lower=lower, upper=upper )
    parameter_estimates[['iterations']] = parameter_estimates[['iterations']] + Temp[['iterations']]
    parameter_estimates[['evaluations']] = parameter_estimates[['evaluations']] + Temp[['evaluations']]
  }

  ## Run some Newton steps
  for(i in seq_len(newtonsteps)) {
    g = as.numeric( gr(parameter_estimates$par) )
    h = optimHess(parameter_estimates$par, fn=fn, gr=gr)
    parameter_estimates$par = parameter_estimates$par - solve(h, g)
    parameter_estimates$objective = fn(parameter_estimates$par)
  }

  # Exclude difficult-to-interpret messages
  parameter_estimates = parameter_estimates[c('par','objective','iterations','evaluations')]

  # Add diagnostics
  parameter_estimates[["time_for_MLE"]] = Sys.time() - start_time
  parameter_estimates[["max_gradient"]] = max(abs(gr(parameter_estimates$par)))
  parameter_estimates[["Convergence_check"]] = ifelse( parameter_estimates[["max_gradient"]]<0.0001, "There is no evidence that the model is not converged", "The model is likely not converged" )
  parameter_estimates[["number_of_coefficients"]] = c("Total"=length(unlist(obj$env$parameters)), "Fixed"=length(startpar), "Random"=length(unlist(obj$env$parameters))-length(startpar) )
  parameter_estimates[["AIC"]] = TMBAIC( opt=parameter_estimates )
  if( n!=Inf ){
    parameter_estimates[["AICc"]] = TMBAIC( opt=parameter_estimates, n=n )
    parameter_estimates[["BIC"]] = TMBAIC( opt=parameter_estimates, p=log(n) )
  }
  parameter_estimates[["diagnostics"]] = data.frame( "Param"=names(startpar), "starting_value"=startpar, "Lower"=lower, "MLE"=parameter_estimates$par, "Upper"=upper, "final_gradient"=as.vector(gr(parameter_estimates$par)) )

  # Get standard deviations
  if(getsd==TRUE){
    # Compute hessian
    sd_time = Sys.time()
    h = optimHess(parameter_estimates$par, fn=fn, gr=gr)
    # Check for problems
    if( is.character(try(chol(h),silent=TRUE)) ){
      warning("Hessian is not positive definite, so standard errors are not available")
      if( !is.null(savedir) ){
        capture.output( parameter_estimates, file=file.path(savedir,"parameter_estimates.txt"))
      }
      return( list("opt"=parameter_estimates, "h"=h) )
    }
    # Compute standard errors
    if( bias.correct==FALSE | is.null(BS.control[["vars_to_correct"]]) ){
      if( !is.null(BS.control[["nsplit"]]) ) {
        if( BS.control[["nsplit"]] == 1 ) BS.control[["nsplit"]] = NULL
      }
      parameter_estimates[["SD"]] = TMB::sdreport( obj=obj, par.fixed=parameter_estimates$par, hessian.fixed=h, bias.correct=bias.correct,
        bias.correct.control=BS.control[c("sd","split","nsplit")], getReportCovariance=getReportCovariance,
        getJointPrecision=getJointPrecision, ... )
    }else{
      if( "ADreportIndex" %in% names(obj$env) ){
        Which = as.vector(unlist( obj$env$ADreportIndex()[ BS.control[["vars_to_correct"]] ] ))
      }else{
        # Run first time to get indices
        parameter_estimates[["SD"]] = TMB::sdreport( obj=obj, par.fixed=parameter_estimates$par, hessian.fixed=h, bias.correct=FALSE,
           getReportCovariance=FALSE, getJointPrecision=FALSE, ... )
        # Determine indices
        Which = which( rownames(summary(parameter_estimates[["SD"]],"report")) %in% BS.control[["vars_to_correct"]] )
      }
      # Split up indices
      if( !is.null(BS.control[["nsplit"]]) && BS.control[["nsplit"]]>1 ){
        Which = split( Which, cut(seq_along(Which), BS.control[["nsplit"]]) )
      }
      Which = Which[sapply(Which,FUN=length)>0]
      if(length(Which)==0) Which = NULL
      # Repeat SD with indexing
      message( paste0("Bias correcting ", length(Which), " derived quantities") )
      parameter_estimates[["SD"]] = TMB::sdreport( obj=obj, par.fixed=parameter_estimates$par, hessian.fixed=h, bias.correct=TRUE,
        bias.correct.control=list(sd=BS.control[["sd"]], split=Which, nsplit=NULL), getReportCovariance=getReportCovariance,
        getJointPrecision=getJointPrecision, ... )
    }
    # Update
    parameter_estimates[["Convergence_check"]] = ifelse( parameter_estimates$SD$pdHess==TRUE, parameter_estimates[["Convergence_check"]], "The model is definitely not converged" )
    parameter_estimates[["time_for_sdreport"]] = Sys.time() - sd_time
    # Add hessian to parameter_estimates
    if( getHessian==TRUE ){
      parameter_estimates[["hessian"]] = h
    }
  }
  parameter_estimates[["time_for_run"]] = Sys.time() - start_time + start_time_elapsed

  # Save results
  if( !is.null(savedir) ){
    save( parameter_estimates, file=file.path(savedir,"parameter_estimates.RData"))
    capture.output( parameter_estimates, file=file.path(savedir,"parameter_estimates.txt"))
  }

  # Print warning to screen
  if( quiet==FALSE & parameter_estimates[["Convergence_check"]] != "There is no evidence that the model is not converged" ){
    message( "#########################" )
    message( parameter_estimates[["Convergence_check"]] )
    message( "#########################" )
  }

  # Return stuff
  return( parameter_estimates )
}

