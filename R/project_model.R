#' Project a fitted VAST model forward in time
#'
#' \code{project_model} simulates random effects forward in time, for use to
#'        generate a predictive interval without actually re-fitting the model.
#'        This is useful e.g., to generate end-of-century projections.
#'
#' The function specifically simulates new values for random effects occurring
#'      during forecasted years.  This includes some combination of intercepts
#'      \code{beta1} or \code{beta2} and spatio-temporal terms
#'      \code{epsilon1} or \code{epsilon2} depending on which
#'      are treated as random during estimation.  It does *not* generate new values of
#'      covariates or random-effects that are not indexed by time \code{omega1} or \code{omega2}
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
#' @param n_samples Number of samples to include.  If \code{n_samples=1} then \code{project_model}
#'        just returns the list of REPORTed variables.  If \code{n_samples>1} then \code{project_model}
#'        returns a list of lists, where each element is the list of REPORTed variables.
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
          n_samples = 1,
          new_covariate_data = NULL,
          historical_uncertainty = "both",
          seed = 123456,
          working_dir = paste0(getwd(),"/"),
          what = NULL ){

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

  # Errors
  if( any(x$data_list$RhoConfig %in% c(0,3)) ){
    #stop("`project_model` is currently designed to work only with temporally varying epsilon and beta terms")
  }
  if( any(x$data_list$RhoConfig[c("Beta1","Beta2")] %in% c(0)) ){
    stop("`project_model` is currently designed to work only with temporally varying or constant beta terms")
  }
  if( convert_version_name(x$settings$Version) < convert_version_name("VAST_v14_0_1") ){
    stop("`project_model` requires version >= `VAST_v14_0_1`")
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
    u_zr = rmvnorm_prec( mu=Obj$env$last.par.best, prec=Sdreport$jointPrecision, n.sims=n_samples, seed=seed)
  }else if( historical_uncertainty == "random" ){
    # Retape and call once to get last.par.best to work
    Obj$retape()
    Obj$fn(x$parameter_estimates$par)
    u_zr = Obj$env$last.par.best %o% rep(1, n_samples)
    # Simulate random effects
    set.seed(seed)
    MC = Obj$env$MC( keep=TRUE, n=n_samples, antithetic=FALSE )
    u_zr[Obj$env$random,] = attr(MC, "samples")
    #Hess = Obj$env$spHess(par=u_z, random=TRUE)
    #u_z[-Obj$env$lfixed()] = rmvnorm_prec( mu=u_z[-Obj$env$lfixed()], prec=Hess, n.sims=1, seed=seed)[,1]
  }else if( historical_uncertainty == "none" ){
    u_zr = Obj$env$last.par.best %o% rep(1, n_samples)
  }else{
    stop("Check `historical_uncertainty` argument")
  }

  ##############
  # Step 2: Generate uncertainty in historical period
  ##############

  t_i = c( x$data_frame$t_i, max(x$data_frame$t_i)+rep(seq_len(n_proj),each=2) )
  b_i = c( x$data_list$b_i, as_units(rep(c(0,mean(x$data_frame$b_i)),n_proj), units(x$data_list$b_i)) )
  v_i = c( x$data_frame$v_i, rep(0,2*n_proj) )
  Lon_i = c( x$data_frame$Lon_i, rep(mean(x$data_frame$Lon_i),2*n_proj) )
  Lat_i = c( x$data_frame$Lat_i, rep(mean(x$data_frame$Lat_i),2*n_proj) )
  a_i = c( x$data_list$a_i, as_units(rep(mean(x$data_frame$a_i),2*n_proj), units(x$data_list$a_i)) )
  PredTF_i = c( x$data_list$PredTF_i, rep(1,2*n_proj) )
  c_iz = rbind( x$data_list$c_iz, x$data_list$c_iz[rep(seq_len(n_proj),each=2),,drop=FALSE] )
  new_catchability_data = rbind( x$catchability_data, x$catchability_data[rep(seq_len(n_proj),each=2),,drop=FALSE] )
  proj_t = x$data_list$n_t + seq_len( n_proj )

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

  # Object to keep output
  out = vector("list", length=n_samples)

  # Loop through 1:n_samples
  for( sampleI in seq_len(n_samples) ){
    ##############
    # Step 4: Merge ParList and ParList1
    ##############

    # Get full size
    #ParList1 = x1$tmb_list$Obj$env$parList()
    ParList1 = x1$tmb_list$Parameters

    # Get ParList
    ParList = Obj$env$parList( par = u_zr[,sampleI] )

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

    # Deal with beta1/beta2 = 3
    if( x$data_list$RhoConfig["Beta1"]==3 ){
      tmp = ParList1$beta1_ft
      tmp[,proj_t] = NA
      ParList1$beta1_ft = ifelse( is.na(tmp), rowMeans(tmp,na.rm=TRUE)%o%rep(1,ncol(tmp)), ParList1$beta1_ft )
    }
    if( x$data_list$RhoConfig["Beta2"]==3 ){
      tmp = ParList1$beta2_ft
      tmp[,proj_t] = NA
      ParList1$beta2_ft = ifelse( is.na(tmp), rowMeans(tmp,na.rm=TRUE)%o%rep(1,ncol(tmp)), ParList1$beta2_ft )
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

    # Debugging comparisons
    if( FALSE ){
      # Samples
      range(u_zr[,1] - Obj$env$last.par.best)
      # Compare fixed effects
      v0 - x$tmb_list$Obj$env$last.par.best[-x$tmb_list$Obj$env$random]
      v2 = x2$tmb_list$Obj$env$last.par.best[-x2$tmb_list$Obj$env$random]
      range(v0 - v2)
      # Compare parameters
      v0 = unlist(x$tmb_list$Parameters[-match(x$tmb_list$Random,names(ParList))])
      #v1 = unlist(ParList1[-match(x1$tmb_list$Random,names(ParList))])
      v2 = unlist(x2$tmb_list$Parameters[-match(x2$tmb_list$Random,names(ParList))])
      range(v1-v2)
      # compare last.par.best
      v0 = x$tmb_list$Obj$env$last.par.best[-x$tmb_list$Obj$env$random]
      v2 = x2$tmb_list$Obj$env$last.par.best[-x2$tmb_list$Obj$env$random]
      range( v0 - v2 )
    }

    ##############
    # Step 5: Simulate random effects
    ##############

    # Simulate Epsiloninput / Betainput for projection years
    x2$tmb_list$Obj$env$data$Options_list$simulate_t[] = c( rep(0,x$data_list$n_t), rep(1,n_proj) )

    # Simulate type=1 so Omegas and other random effects are held fixed
    out[[sampleI]] = simulate_data( fit = x2,
                         type = 1,
                         random_seed = NULL )

    # Amend labels
    x2$Report = out[[sampleI]]
    out[[sampleI]] = amend_output(x2)

    # Subset to values specified in "what"
    if( !is.null(what) ){
      out[[sampleI]] = out[[sampleI]][which(names(out[[sampleI]]) %in% what)]
    }
  }

  if( n_samples==1 ){
    out = out[[1]]
  }
  return(out)
}

