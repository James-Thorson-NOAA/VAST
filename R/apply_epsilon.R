
#' Custom epsilon-correct method
#'
#' \code{apply_epsilon} uses updates to TMB to implement a faster calculation for epsilon-correction
#'
#' @inheritParams TMB::MakeADFun
#' @param fit output from \code{\link[FishStatsUtils]{fit_model}}, specifically using
#'        slotes \code{tmb_list}, \code{input_args}, \code{parameter_estimates$SD}
#' @param ADREPORT_name string indicating name of ADREPORT'ed variable
#' @param eps_name string indicating name of PARAMETER used internally by TMB
#'        for calculating desired gradient
#'
#' @return Standard output from \code{\link[TMB]{sdreport}}, but with slot
#'         \code{x$unbiased} added if needed, and adding or replacing values for
#'         \code{x$unbiased$value} corresponding to \code{ADREPORT_name}
#'
#' @export
apply_epsilon <-
function( fit,
          ADREPORT_name = "Index_ctl",
          eps_name = "eps_Index_ctl",
          inner.control = list(sparse=TRUE, lowrank=TRUE, trace=FALSE) ){

  # Extract stuff
  Obj = fit$tmb_list$Obj
  if(is.null(fit$Report)) fit$Report = Obj$report()
  if(is.null(fit$ParHat)) fit$ParHat = Obj$env$parList()

  # Error checks
  if( !all(c("tmb_list","input_args","parameter_estimates") %in% names(fit)) ) stop("Check `fit` in `apply_epsilon`: function is only designed to work with output from VAST using `fit_model`")
  if( !(ADREPORT_name %in% names(fit$Report)) ) stop("Check `ADREPORT_name` in `apply_epsilon`")
  if( fit$input_args$model_args_input$framework!="TMBad" ) stop("`apply_epsilon` requires that the CPP be compiled using framework=`TMBad`")
  if( is.null(fit$parameter_estimates$SD) ) stop("Please re-run with `getsd=TRUE`")

  # Simple extractions
  Data = Obj$env$data
  Map = Obj$env$map
  Random = fit$tmb_list$Random

  # Extract and modify parameters
  New_params = fit$ParHat
  New_params[[eps_name]] = array(0, dim=dim(fit$Report[[ADREPORT_name]]) )

  # Change MLE
  fixed = fit$parameter_estimates$par
  new_values = rep( 0, prod(dim(New_params[[eps_name]])) )
  names(new_values) = rep( eps_name, length(new_values))
  fixed = c( fixed, new_values )

  # detect sparse + lowrank hessian ... appears to freeze with lowrank=FALSE
  obj = MakeADFun( data = Data,
                    parameters = New_params,
                    map = Map,
                    random = Random,
                    intern = TRUE,
                    inner.control = inner.control )
  obj$env$beSilent()
  gradient = obj$gr(fixed)

  # Expand SD
  SD = fit$parameter_estimates$SD
  if( is.null(SD$unbiased) ){
    SD$unbiased = list( "value"=SD$value, "sd"=NA, "cov"=array(NA,c(1,1)) )
    SD$unbiased$value[] = NA
  }else{
    if( any(!is.na(SD$unbiased$value[ADREPORT_name])) ){
      warning( paste0("it appears that `", ADREPORT_name,"` is already bias-corrected; using `apply_epsilon` seems inefficient and will replace existing values") )
    }
  }

  # Merge gradient into SD
  i = which( names(SD$value) == ADREPORT_name )
  j = which( names(obj$par) == eps_name )
  if( length(i)==length(j) ){
    SD$unbiased$value[i] = gradient[j]
  }else{
    warning("Check `apply_epsilon` for bugs")
  }

  # return SD
  class(SD) = "sdreport"
  return(SD)
}
