
#' Calculate effects for plotting
#'
#' @title Adapts package \code{effects}
#'
#' @rawNamespace S3method(effects::Effect, fit_model)
#' @export
Effect.fit_model = function (focal.predictors, mod, which_formula="X1", ...) {

  # Error checks
  if( mod$data_list$n_c>1 ){
    stop("`Effect.fit_model` is not currently designed for multivariate models")
  }
  if( !all(c("covariate_data_full","catchability_data_full") %in% ls(.GlobalEnv)) ){
    stop("Please load `covariate_data_full` and `catchability_data_full` into global memory")
  }

  # Basic
  if( which_formula=="X1" ){
    formula_orig = mod$X1_formula
    parname = "gamma1_cp"
    mod$call = mod$effects$call_X1
  }else if( which_formula=="X2" ){
    formula_orig = mod$X2_formula
    parname = "gamma2_cp"
    mod$call = mod$effects$call_X2
  }else{
    stop("Other options are not implemented")
  }

  # Extract parameters
  whichnum = which(names(mod$parameter_estimates$par)==parname)
  mod$parhat = mod$parameter_estimates$par[whichnum]
    mod$covhat = mod$parameter_estimates$SD$cov.fixed[whichnum,whichnum]
  # Fill in values that are mapped off
  if( parname %in% names(mod$tmb_list$Obj$env$map) ){
    mod$parhat = mod$parhat[ mod$tmb_list$Obj$env$map[[parname]] ]
      mod$covhat = mod$covhat[ mod$tmb_list$Obj$env$map[[parname]], mod$tmb_list$Obj$env$map[[parname]] ]
    mod$parhat = ifelse( is.na(mod$parhat), 0, mod$parhat)
      mod$covhat = ifelse( is.na(mod$covhat), 0, mod$covhat)
  }
  # add names
  names(mod$parhat)[] = parname
    rownames(mod$covhat) = colnames(mod$covhat) = names(mod$parhat)

  # Augment stuff
  formula_full = update.formula(formula_orig, linear_predictor~.+0)
  mod$coefficients = mod$parhat
  mod$vcov = mod$covhat
  mod$formula = formula_full
  mod$family = gaussian(link = "identity")

  # Functions for package
  family.fit_model = function(x,...) x$family
  vcov.fit_model = function(x,...) x$vcov

  # Do stuff
  fam = mod$family
  ## dummy functions to make Effect.default work
  dummyfuns = list(variance = function(mu) mu,
                    initialize = expression(mustart = y + 0.1),
                    dev.resids = function(...) poisson()$dev.res(...) )
  for( i in names(dummyfuns) ){
    if( is.null(fam[[i]]) ) fam[[i]] = dummyfuns[[i]]
  }
  ## allow calculation of effects ...
  if (length(formals(fam$variance))>1) {
    warning("overriding variance function for effects: computed variances may be incorrect")
    fam$variance = dummyfuns$variance
  }
  args = list(call = mod$call,
               coefficients = mod$coefficients,
               vcov = mod$vcov,
               family = fam,
               formula = formula_full)
  if(!requireNamespace("effects")) stop("please install the effects package")
  effects::Effect.default(focal.predictors,
    mod,
    ...,
    sources = args)
}
