Build_TMB_Fn <-
function( TmbData, Version, Q_Config=TRUE, CovConfig=TRUE,
  RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0),
  ConvergeTol=1, Use_REML=FALSE, loc_x=NULL, Parameters="generate", Random="generate", Map="generate",
  DiagnosticDir=NULL, TmbDir=system.file("executables", package="SpatialCompData"), RunDir=getwd() ){
                                            
  # Compile TMB software
  #dyn.unload( paste0(RunDir,"/",dynlib(TMB:::getUserDLL())) ) # random=Random,
  file.copy( from=paste0(TmbDir,"/",Version,".cpp"), to=paste0(RunDir,"/",Version,".cpp"), overwrite=FALSE)
  setwd( RunDir )
  compile( paste0(Version,".cpp") )

  # Local functions
  boundsifpresent_fn = function( par, map, name, lower, upper, bounds ){
    if( name %in% names(par) ){
      bounds[grep(name,names(par)),c('Lower','Upper')] = c(lower,upper)
    }
    return( bounds )
  }
  
  # Parameters
    # DataList=TmbData
  if( length(Parameters)==1 && Parameters=="generate" ) Parameters = Param_Fn( Version=Version, DataList=TmbData, RhoConfig=RhoConfig )

  # Which are random
  if( length(Random)==1 && Random=="generate" ) Random = c("Epsiloninput1_sct", "Omegainput1_sc", "Epsiloninput2_sct", "Omegainput2_sc")
  if( RhoConfig[["Beta1"]]!=0 ) Random = c(Random, "beta1_ct")
  if( RhoConfig[["Beta2"]]!=0 ) Random = c(Random, "beta2_ct")
  if( Use_REML==TRUE ){
    Random = union(Random, c("beta1_ct","gamma1_j","gamma1_tp","lambda1_k","beta2_ct","gamma2_j","gamma2_tp","lambda2_k"))
    Random = Random[which(Random %in% names(Parameters))]
  }

  # Which parameters are turned off
  if( length(Map)==1 && Map=="generate" ) Map = Make_Map( Version=Version, TmbData=TmbData, CovConfig=CovConfig, Q_Config=Q_Config, RhoConfig=RhoConfig, Aniso=TmbData[['Options_vec']]['Aniso'])

  # Build object
  dyn.load( paste0(RunDir,"/",dynlib(Version)) ) # random=Random,
  Obj <- MakeADFun(data=TmbData, parameters=Parameters, hessian=FALSE, map=Map, random=Random, inner.method="newton")
  Obj$control <- list(trace=1, parscale=1, REPORT=1, reltol=1e-12, maxit=100)

  # Diagnostic functions (optional)
  if( !is.null(DiagnosticDir) ){
    Obj$gr_orig = Obj$gr
    Obj$fn_orig = Obj$fn
    Obj$fn = function( vec ){
      capture.output( matrix(vec,ncol=1,dimnames=list(names(Obj$par),NULL)), file=paste0(DiagnosticDir,"fn.txt") )
      write.table( matrix(vec,nrow=1), row.names=FALSE, sep=",", col.names=FALSE, append=TRUE, file=paste0(DiagnosticDir,"trace.csv"))
      return( Obj$fn_orig(vec) )
    }
    Obj$gr = function( vec ){
      capture.output( matrix(vec,ncol=1,dimnames=list(names(Obj$par),NULL)), file=paste0(DiagnosticDir,"gr.txt") )
      return( Obj$gr_orig(vec) )
    }
    write.table( matrix(Obj$par,nrow=1), row.names=FALSE, sep=",", col.names=FALSE, file=paste0(DiagnosticDir,"trace.csv"))
  }
  
  # Declare upper and lower bounds for parameter search
  Bounds = matrix( NA, ncol=2, nrow=length(Obj$par), dimnames=list(names(Obj$par),c("Lower","Upper")) )
  Bounds[,'Lower'] = rep(-50, length(Obj$par))
  Bounds[,'Upper'] = rep( 50, length(Obj$par))
  Bounds[grep("logtau",names(Obj$par)),'Upper'] = 10   # Version < v2i
  Bounds[grep("logeta",names(Obj$par)),'Upper'] = log(1/(1e-2*sqrt(4*pi))) # Version >= v2i: Lower bound on margSD = 1e-4
  Bounds[grep("SigmaM",names(Obj$par)),'Upper'] = 10 # ZINB can crash if it gets > 20
  if( !is.null(loc_x) ){
    Dist = dist(loc_x)
    Bounds[grep("logkappa",names(Obj$par)),'Lower'] = log( sqrt(8)/max(Dist) ) # Range = nu*sqrt(8)/kappa
    Bounds[grep("logkappa",names(Obj$par)),'Upper'] = log( sqrt(8)/min(Dist) ) # Range = nu*sqrt(8)/kappa
  }
  Bounds = boundsifpresent_fn( par=Obj$par, name="gamma1", lower=-20, upper=20, bounds=Bounds)
  Bounds = boundsifpresent_fn( par=Obj$par, name="gamma2", lower=-20, upper=20, bounds=Bounds)
  Bounds = boundsifpresent_fn( par=Obj$par, name="lambda1", lower=-20, upper=20, bounds=Bounds)
  Bounds = boundsifpresent_fn( par=Obj$par, name="lambda2", lower=-20, upper=20, bounds=Bounds)
  Bounds = boundsifpresent_fn( par=Obj$par, name="Beta_rho1", lower=-0.99, upper=0.99, bounds=Bounds)
  Bounds = boundsifpresent_fn( par=Obj$par, name="Beta_rho2", lower=-0.99, upper=0.99, bounds=Bounds)
  Bounds = boundsifpresent_fn( par=Obj$par, name="Epsilon_rho1", lower=-0.99, upper=0.99, bounds=Bounds)
  Bounds = boundsifpresent_fn( par=Obj$par, name="Epsilon_rho2", lower=-0.99, upper=0.99, bounds=Bounds)
  Bounds = boundsifpresent_fn( par=Obj$par, name="rho_c1", lower=-0.99, upper=0.99, bounds=Bounds)
  Bounds = boundsifpresent_fn( par=Obj$par, name="rho_c2", lower=-0.99, upper=0.99, bounds=Bounds)

  # Change convergence tolerance
  Obj$env$inner.control$step.tol <- c(1e-8,1e-12,1e-15)[ConvergeTol] # Default : 1e-8  # Change in parameters limit inner optimization
  Obj$env$inner.control$tol10 <- c(1e-6,1e-8,1e-12)[ConvergeTol]  # Default : 1e-3     # Change in pen.like limit inner optimization
  Obj$env$inner.control$grad.tol <- c(1e-8,1e-12,1e-15)[ConvergeTol] # # Default : 1e-8  # Maximum gradient limit inner optimization

  # Return stuff
  Return = list("Obj"=Obj, "Upper"=Bounds[,'Upper'], "Lower"=Bounds[,'Lower'], "Parameters"=Parameters, "Map"=Map, "Random"=Random)
  return( Return )
}
