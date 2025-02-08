
#' Calculate marginal AIC for a fitted model
#'
#' \code{TMBAIC} calculates AIC for a given model fit
#'
#' @param opt the output from \code{nlminb} or \code{optim}
#' @param p the penalty on additional fixed effects (default=2, for AIC)
#' @param n the sample size, for use in AICc calculation (default=Inf, for which AICc=AIC)
#'
#' @return AIC, where a parsimonious model has a AIC relative to other candidate models

#' @export
TMBAIC=function(opt, p=2, n=Inf){
  k = length(opt[["par"]])
  if( all(c("par","objective") %in% names(opt)) ) negloglike = opt[["objective"]]
  if( all(c("par","value") %in% names(opt)) ) negloglike = opt[["value"]]
  Return = p*k + 2*negloglike + 2*k*(k+1)/(n-k-1)
  return( Return )
}

