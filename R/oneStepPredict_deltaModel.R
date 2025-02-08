
#' Calculate one-step-ahead (OSA) residuals for a mixed-effects delta-model.
#'
#' \code{oneStepPredict_deltaModel} is a wrapper for \code{oneStepPredict} for distributions with a mixture of discrete and continuous distributions
#'
#' It is convenient to compute one-step-ahead residuals for data that arise as a mixture of continuous and discrete distributions.
#' One common example is a delta-model, which arises as a mixture of an encounter probability and a continuous distribution for
#' biomass given an encounter.  In these cases, it is possible to apply \code{\link[TMB]{oneStepPredict}} twice, once for
#' observations falling within the continuous domain, and again for observations in the discrete domain, and then combining the two.
#' This function provides an example of doing so.  It is designed to use the `method="cdf"` feature in \code{\link[TMB]{oneStepPredict}}.
#' This example also shows a proof-of-concept for uniform residuals under a (sufficiently-close-to) correctly specified model.
#'
#' @inheritParams TMB::oneStepPredict
#' @param deltaSupport integer-vector, listing values that have a dirac-delta within an otherwise continuous distribution
#' @param ... list of arguments to pass to \code{\link[TMB]{oneStepPredict}}
#'
#' @return the standard output from \code{\link[TMB]{oneStepPredict}}
#'
#' @examples
#' \dontrun{
#'
#' library(TMB)
#' library(RandomFields)
#' library(INLA) # FROM: http://www.r-inla.org/download
#'
#' ###################
#' # Poisson-link gamma distribution
#' ###################
#'
#' # n = numbers density
#' # w = weight-per-number
#' # cv = CV of gamma
#' dpoislinkgamma = function(x, n, w, cv){
#'   pow = function(a,b) a^b
#'   enc_prob = 1 - exp(-n)
#'   posmean = n * w / enc_prob
#'   if( x==0 ){
#'     dens = 1 - enc_prob
#'   }else{
#'     dens = enc_prob * dgamma(x, shape=pow(cv,-2), scale=posmean*pow(cv,2))
#'   }
#'   if(log==FALSE) return(dens)
#'   if(log==TRUE) return(log(dens))
#' }
#' ppoislinkgamma = function(x, n, w, cv){
#'   pow = function(a,b) a^b
#'   enc_prob = 1 - exp(-n)
#'   posmean = n * w / enc_prob
#'   dist = 1 - enc_prob
#'   if( x>0 ){
#'     posmean = n*w
#'     dist = dist + enc_prob * pgamma(x, shape=pow(cv,-2), scale=posmean*pow(cv,2))
#'   }
#'   return(dist)
#' }
#' rpoislinkgamma = function(n, w, cv){
#'   pow = function(a,b) a^b
#'   enc_prob = 1 - exp(-n)
#'   posmean = n * w / enc_prob
#'   enc = rbinom(n=1, prob=enc_prob, size=1)
#'   x = enc * rgamma(n=1, shape=pow(cv,-2), scale=posmean*pow(cv,2))
#'   return(x)
#' }
#'
#' ###################
#' # Simulate data
#' ###################
#'
#' Dim = c("n_x"=10, "n_y"=10)
#' loc_xy = expand.grid("x"=1:Dim['n_x'], "y"=1:Dim['n_y'])
#' Scale = 2
#' Sigma2 = (0.5) ^ 2
#' beta0 = 1
#' w = 1
#' cv = 0.1
#'
#' # Simulate spatial process
#' RMmodel = RMgauss(var=Sigma2, scale=Scale)
#' epsilon_xy = array(RFsimulate(model=RMmodel, x=loc_xy[,'x'], y=loc_xy[,'y'])@data[,1], dim=Dim)
#'
#' # Simulate samples
#' c_xy = array(NA, dim=dim(epsilon_xy))
#' for(x in 1:nrow(c_xy)){
#' for(y in 1:ncol(c_xy)){
#'   c_xy[x,y] = rpoislinkgamma( n=exp(beta0 + epsilon_xy[x,y]), w=w, cv=cv )
#' }}
#'
#' #' ###################
#' #' # SPDE-based
#' ###################
#'
#' # create mesh
#' mesh = inla.mesh.create( loc_xy, plot.delay=NULL, refine=FALSE)
#' # Create matrices in INLA
#' spde <- inla.spde2.matern(mesh, alpha=2)
#'
#' # COmpile
#' compile( "deltaModel.cpp" )
#' dyn.load( dynlib("deltaModel") )
#'
#' # Build object
#' Data = list("c_i"=as.vector(c_xy), "j_i"=mesh$idx$loc-1, "M0"=spde$param.inla$M0, "M1"=spde$param.inla$M1, "M2"=spde$param.inla$M2 )
#' Params = list( "beta0"=0, "ln_tau"=0, "ln_kappa"=0, "ln_w"=0, "ln_cv"=0, "epsilon_j"=rep(0,nrow(spde$param.inla$M0)) )
#' Map = list( "ln_tau"=factor(NA), "ln_kappa"=factor(NA), "epsilon_j"=factor(rep(NA,length(Params$epsilon_j))) )
#' Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_j", map=Map )
#'
#' # Optimize
#' Opt = fit_tmb( obj=Obj, newtonsteps=0, getsd=FALSE )
#' report = Obj$report()
#'
#' # Run
#' osa = oneStepPredict_deltaModel( obj=Obj, observation.name="c_i", method="cdf",
#'   data.term.indicator="keep", deltaSupport=0, trace=TRUE, seed=1 ) #discreteSupport = seq(0,max(Data$c_i),by=1) )
#' qqnorm(osa$residual); abline(0,1)
#'
#' # should be uniform from 0 to mean(c_xy==0) when mapping off random effects
#' qresid = NULL
#' for(i in 1:1000){
#'   osa = oneStepPredict_deltaModel( obj=Obj, observation.name="c_i", method="cdf",
#'     data.term.indicator="keep", deltaSupport=0, trace=FALSE, seed=i ) #discreteSupport = seq(0,max(Data$c_i),by=1) )
#'   qresid = c( qresid, pnorm(osa[which(Obj$env$data[["c_i"]]==0),'residual']) )
#' }
#' hist(qresid)
#' abline( v=mean(c_xy==0), lwd=3, lty="dotted" )
#'
#' }
#'
#' @export
oneStepPredict_deltaModel = function( obj, observation.name, deltaSupport=0, ... ){
  osa1 = oneStepPredict( obj=obj, discrete=TRUE, observation.name=observation.name, ... )
  osa = osa2 = oneStepPredict( obj=obj, discrete=FALSE, observation.name=observation.name, ... )
  osa[ which(obj$env$data[[observation.name]]%in%deltaSupport), na.omit(match(c("U","residual"),colnames(osa)))] =
    osa1[ which(obj$env$data[[observation.name]]%in%deltaSupport), na.omit(match(c("U","residual"),colnames(osa)))]
  return(osa)
}



