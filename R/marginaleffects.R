#' @method get_coef fit_model
#' @export
get_coef.fit_model = function(model, covariate, ...){
  # Select covariate
  if(covariate=="X1") param = "gamma1_cp"
  if(covariate=="X2") param = "gamma2_cp"
  if(covariate=="Q1") param = "lambda1_k"
  if(covariate=="Q2") param = "lambda2_k"
  # extract params
  out = model$last.par.best[grep(param,names(model$last.par.best))]
  return(out)
}

#' @method get_vcov fit_model
#' @export
get_vcov.fit_model = function(model, covariate, ...){
  if(covariate=="X1") param = "gamma1_cp"
  if(covariate=="X2") param = "gamma2_cp"
  if(covariate=="Q1") param = "lambda1_k"
  if(covariate=="Q2") param = "lambda2_k"
  # Extract covariance
  whichrows = which( names(model$parameter_estimates$par) == param )
  if( is.null(model$parameter_estimates$SD) ){
    out = NULL
  }else{
    out = array(model$parameter_estimates$SD$cov.fixed[whichrows,whichrows],dim=rep(length(whichrows),2))
  }
  return(out)
}

#' @method set_coef fit_model
#' @export
set_coef.fit_model = function(model, newpar, covariate, ...){
  if(covariate=="X1") param = "gamma1_cp"
  if(covariate=="X2") param = "gamma2_cp"
  if(covariate=="Q1") param = "lambda1_k"
  if(covariate=="Q2") param = "lambda2_k"
  # modify and return params
  if( length(newpar) != length(grep(param,names(model$last.par.best))) ){
    stop("Check length of 'newpar'")
  }
  model$last.par.best[grep(param,names(model$last.par.best))] <- newpar
  return(model)
}

#' @method get_predict fit_model
#' @export
get_predict.fit_model = function(model, newdata, covariate, center=FALSE, ...){
  # update formula (following logic in make_covariates)
  if(covariate=="X1"){
    formula = update.formula(model$X1_formula, ~.+1)
    param = "gamma1_cp"
    data = model$covariate_data
  }
  if(covariate=="X2"){
    formula = update.formula(model$X2_formula, ~.+1)
    param = "gamma2_cp"
    data = model$covariate_data
  }
  if(covariate=="Q1"){
    formula = update.formula(model$Q1_formula, ~.+1)
    param = "lambda1_k"
    data = model$catchability_data
  }
  if(covariate=="Q2"){
    formula = update.formula(model$Q2_formula, ~.+1)
    param = "lambda2_k"
    data = model$catchability_data
  }

  # build original model.frame
  frame0 = model.frame( formula=formula, data=data )
  terms0 = terms( frame0 )
  xlevels = .getXlevels( terms0, frame0 )

  # get new design matrix
  terms1 = delete.response( terms0 )
  frame1 = model.frame( terms1, newdata, xlev=xlevels )
  X_ip = model.matrix( terms1, frame1 )

  # Drop intercept (following logic in make_covariates)
  Columns_to_keep = which( attr(X_ip,"assign") != 0 )
  X_ip = X_ip[,Columns_to_keep,drop=FALSE]

  # Multiply and center
  ParHat = model$tmb_list$Obj$env$parList( par=model$last.par.best )
  gamma_cp = ParHat[[param]]
  yhat_ic = X_ip %*% t(gamma_cp)
  if(center==TRUE) yhat_ic = yhat_ic - outer(rep(1,nrow(yhat_ic)),colMeans(yhat_ic))

  # Return
  if( covariate %in% c("X1","X2") ){
    out = expand.grid( rowid=seq_along(yhat_ic[,1]),
                       category=model$category_names )
    out$estimate = as.vector(yhat_ic)
  }
  if( covariate %in% c("Q1","Q2") ){
    out = data.frame( rowid=seq_along(yhat_ic[,1]),
                      estimate=yhat_ic[,1] )
  }
  return(out)
}
