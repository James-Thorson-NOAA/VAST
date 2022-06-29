#' Project a fitted VAST model forward in time
#'
#' \code{project_model} simulates random effects forward in time, for use to
#'        generate a predictive interval without actually re-fitting the model.
#'        This is useful e.g., to generate end-of-century projections.
#'
#' The function specifically simulates new values for random effects occurring
#'      during forecasted years.  This includes some combination of intercepts
#'      {beta1/beta2} and spatio-temporal terms {epsilon1/epsilon2} depending on which
#'      are treated as random during estimation.  It does *not* generate new values of
#'      covariates or random-effects that are not indexed by time {omega1/omega2}
#'
#' Note that the model may behave poorly when \code{historical_uncertainty="both"}
#'      and the estimation model includes an AR1 process for any component.
#'      Given this combination of features, some samples may have a `rho` value >1
#'      or <1, which will result in exponential growth for any such sampled value.
#'      This behavior could be improved in future code updates by using \code{tmbstan}
#'      instead of the normal approximation to generate parametric uncertainty
#'      during the historical period.
#'
#' Similarly, estimating a RW process for epsilon will result in an exponential increase
#'      in forecasted total abundance over time.  This occurs because the variance across locations
#'      of epsilon increases progressively during the forecast period, such that
#'      the index is again dominated by the forecasted density at a few sites.
#'
#'
#' @param x Output from \code{\link{fit_model}}
#' @param n_proj Number of time-steps to include in projection
#' @param new_covariate_data New covariates to include for future intervals
#' @param historical_uncertainty Whether to incorporate uncertainty about fitted interval
#' \describe{
#'    \item{\code{historical_uncertainty="both"}}{Include uncertainty in fixed and random effects using joint precision matrix}
#'    \item{\code{historical_uncertainty="random"}}{Include uncertainty in random effects using inner Hessian matrix}
#'    \item{\code{historical_uncertainty="none"}}{Condition upon MLE for fixed and Empirical Bayes for random effects}
#' }
#'
#' @return All \code{obj$report()} output for a single simulation of historical period
#'         as well as \code{n_proj} forecast intervals
#'
#' @examples
#' \dontrun{
#' # Run model
#' fit = fit_model( ... )
#'
#' # Add projection
#' project_model( x = fit,
#'                n_proj = 80,
#'                new_covariate_data = NULL,
#'                historical_uncertainty = "both",
#'                seed = NULL )
#' }
#'
#' @export
project_model <-
function( x,
          n_proj,
          new_covariate_data = NULL,
          historical_uncertainty = "both",
          seed = 123456,
          working_dir = paste0(getwd(),"/") ){

  # Unpack
  Obj = x$tmb_list$Obj
  Sdreport = x$parameter_estimates$SD

  # Warnings
  # REVISE: remove historical years from new_covariate_data to avoid new data changing fit in earlier years
  if( is.null(new_covariate_data) ){
    new_covariate_data = x$covariate_data
  }else{
    # Confirm all columns are available
    if( !all(colnames(x$covariate_data) %in% colnames(new_covariate_data)) ){
      stop("Please ensure that all columns of `x$covariate_data` are present in `new_covariate_data`")
    }
    # Eliminate unnecessary columns
    new_covariate_data = new_covariate_data[,match(colnames(x$covariate_data),colnames(new_covariate_data))]
    # Eliminate old-covariates that are also present in new_covariate_data
    NN = RANN::nn2( query=x$covariate_data[,c('Lat','Lon','Year')], data=new_covariate_data[,c('Lat','Lon','Year')], k=1 )
    if( any(NN$nn.dist==0) ){
      x$covariate_data = x$covariate_data[-which(NN$nn.dist==0),,drop=FALSE]
    }
  }

  ##############
  # Step 1: Generate uncertainty in historical period
  ##############

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
  if( historical_uncertainty == "both" ){
    u_z = rmvnorm_prec( mu=Obj$env$last.par.best, prec=Sdreport$jointPrecision, n.sims=1, seed=seed)[,1]
  }else if( historical_uncertainty == "random" ){
    # Retape and call once to get last.par.best to work
    Obj$retape()
    Obj$fn(x$parameter_estimates$par)
    u_z = Obj$env$last.par.best
    # Simulate random effects
    set.seed(seed)
    MC = Obj$env$MC( keep=TRUE, n=1, antithetic=FALSE )
    u_z[Obj$env$random] = attr(MC, "samples")[,1]
    #Hess = Obj$env$spHess(par=u_z, random=TRUE)
    #u_z[-Obj$env$lfixed()] = rmvnorm_prec( mu=u_z[-Obj$env$lfixed()], prec=Hess, n.sims=1, seed=seed)[,1]
  }else if( historical_uncertainty == "none" ){
    u_z = Obj$env$last.par.best
  }else{
    stop("Check `historical_uncertainty` argument")
  }

  # Get ParList
  ParList = Obj$env$parList( par = u_z )

  ##############
  # Step 2: Generate uncertainty in historical period
  ##############

  t_i = c( x$data_frame$t_i, max(x$data_frame$t_i)+rep(1:n_proj,each=2) )
  b_i = c( x$data_list$b_i, as_units(rep(c(0,mean(x$data_frame$b_i)),n_proj), units(x$data_list$b_i)) )
  v_i = c( x$data_frame$v_i, rep(0,2*n_proj) )
  Lon_i = c( x$data_frame$Lon_i, rep(mean(x$data_frame$Lon_i),2*n_proj) )
  Lat_i = c( x$data_frame$Lat_i, rep(mean(x$data_frame$Lat_i),2*n_proj) )
  a_i = c( x$data_list$a_i, as_units(rep(mean(x$data_frame$a_i),2*n_proj), units(fit$data_list$a_i)) )
  PredTF_i = c( x$data_list$PredTF_i, rep(1,2*n_proj) )
  c_iz = rbind( x$data_list$c_iz, x$data_list$c_iz[rep(1:n_proj,each=2),,drop=FALSE] )
  new_catchability_data = rbind( x$catchability_data, x$catchability_data[rep(1:n_proj,each=2),,drop=FALSE] )

  ##############
  # Step 3: Build object with padded bounds
  ##############

  x1 = fit_model( settings = x$settings,
    Lat_i = Lat_i,
    Lon_i = Lon_i,
    t_i = t_i,
    b_i = b_i,
    a_i = a_i,
    v_i = v_i,
    c_iz = c_iz,
    PredTF_i = PredTF_i,
    covariate_data = new_covariate_data,
    X1_formula = x$X1_formula,
    X2_formula = x$X2_formula,
    X1config_cp = x$X1config_cp,
    X2config_cp = x$X2config_cp,
    catchability_data = new_catchability_data,
    Q1config_k = x$Q1config_k,
    Q2config_k = x$Q2config_k,
    Q1_formula = x$Q1_formula,
    Q2_formula = x$Q2_formula,
    build_model = FALSE,
    working_dir = working_dir )

  # Get full size
  #ParList1 = x1$tmb_list$Obj$env$parList()
  ParList1 = x1$tmb_list$Parameters

  ##############
  # Step 4: Merge ParList and ParList1
  ##############

  for( i in seq_along(ParList) ){
    dim = function(x) if(is.vector(x)){return(length(x))}else{return(base::dim(x))}
    dim_match = ( dim(ParList[[i]]) == dim(ParList1[[i]]) )
    if( sum(dim_match==FALSE)==0 ){
      ParList1[[i]] = ParList[[i]]
    }else if( sum(dim_match==FALSE)==1 ){
      dim_list = lapply( dim(ParList[[i]]), FUN=function(x){seq_len(x)} )
      ParList1[[i]][as.matrix(expand.grid(dim_list))] = ParList[[i]][as.matrix(expand.grid(dim_list))]
    }else if( sum(dim_match==FALSE)>=2 ){
      stop("Check matching")
    }
  }

  ##############
  # Step 5: Re-build model
  ##############

  x2 = fit_model( settings = x$settings,
    Lat_i = Lat_i,
    Lon_i = Lon_i,
    t_i = t_i,
    b_i = b_i,
    a_i = a_i,
    v_i = v_i,
    c_iz = c_iz,
    PredTF_i = PredTF_i,
    covariate_data = new_covariate_data,
    X1_formula = x$X1_formula,
    X2_formula = x$X2_formula,
    X1config_cp = x$X1config_cp,
    X2config_cp = x$X2config_cp,
    catchability_data = new_catchability_data,
    Q1config_k = x$Q1config_k,
    Q2config_k = x$Q2config_k,
    Q1_formula = x$Q1_formula,
    Q2_formula = x$Q2_formula,
    run_model = FALSE,
    Parameters = ParList1,
    working_dir = working_dir )

  ##############
  # Step 5: Simulate random effects
  ##############

  # Simulate Epsiloninput / Betainput for projection years
  x2$tmb_list$Obj$env$data$Options_list$simulate_t[] = c( rep(0,x$data_list$n_t), rep(1,n_proj) )

  # Simulate type=1 so Omegas and other random effects are held fixed
  Sim = simulate_data( fit = x2,
                       type = 1,
                       random_seed = NULL )

  # Amend labels
  x2$Report = Sim
  Sim = amend_output(x2)

  return(Sim)
}
