
#' Custom epsilon-correct method
#'
#' \code{apply_epsilon} uses updates to TMB to implement a faster calculation for epsilon-correction
#'
#' @inheritParams TMB::MakeADFun
#' @param fit output from \code{\link[FishStatsUtils]{fit_model}}
#' @param ADREPORT_name string indicating name of ADREPORT'ed variable
#' @param eps_name string indicating name of PARAMETER used internally by TMB
#'        for calculating desired gradient
#' @param data_function function to apply to \code{fit$data_list} prior to passing
#'        to \code{\link[TMB]{MakeADFun}}, e.g., for VAST use \code{strip_units}
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
          data_function = identity,
          inner.control = list(sparse=TRUE, lowrank=TRUE, trace=TRUE) ){

  # Checks
  if(class(fit)!="fit_model") stop("Check `fit` in `apply_epsilon`")
  if(!(ADREPORT_name %in% names(fit$Report))) stop("Check `ADREPORT_name` in `apply_epsilon`")

  # Simple extractions
  Data = fit$data_list
  Random = fit$tmb_list$Random
  Map = fit$tmb_list$Map

  # Extract and modify parameters
  New_params = fit$ParHat
  New_params[[eps_name]] = array(0, dim=dim(fit$Report[[ADREPORT_name]]) )

  # Change MLE
  fixed = fit$par$par
  new_values = rep( 0, prod(dim(New_params[[eps_name]])) )
  names(new_values) = eps_name
  fixed = c( fixed, new_values )

  # detect sparse + lowrank hessian ... appears to freeze with lowrank=FALSE
  obj = MakeADFun( data = lapply(Data, FUN=data_function),
                    parameters = New_params,
                    map = Map,
                    random = Random,
                    intern = TRUE,
                    DLL = fit$settings$Version,
                    inner.control = inner.control )
  obj$env$beSilent()
  gradient = obj$gr(fixed)

  # Expand SD
  SD = fit$par$SD
  if( is.null(SD$unbiased) ){
    SD$unbiased = list( "value"=NA, "sd"=NA, "cov"=array(NA,c(1,1)) )
    SD$unbiased$value = SD$value
    SD$unbiased$value[] = NA
  }else{
    if( any(!is.na(SD$unbiased$value[ADREPORT_name])) ){
      warning( paste0("it appears that `", ADREPORT_name,"` is already bias-corrected; using `apply_epsilon` seems inefficient and will replace existing values") )
    }
  }

  # Merge gradient into SD
  i = which( names(SD$unbiased$value) == ADREPORT_name )
  j = which( names(obj$par) == eps_name )
  if( length(i)==length(j) ){
    SD$unbiased$value[i] = gradient[j]
  }

  # return SD
  class(SD) = "sdreport"
  return(SD)
}
