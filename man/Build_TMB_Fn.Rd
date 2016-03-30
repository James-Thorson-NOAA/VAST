\name{Build_TMB_Fn}
\alias{Build_TMB_Fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Build TMB object and set convergence tolerances and upper/lower limits
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Build_TMB_Fn(TmbData, TmbFile, Version, VesselConfig, CovConfig, Aniso, ConvergeTol = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TmbData}{
%%     ~~Describe \code{TmbData} here~~
}
  \item{TmbFile}{
%%     ~~Describe \code{TmbFile} here~~
}
  \item{Version}{
%%     ~~Describe \code{Version} here~~
}
  \item{VesselConfig}{
%%     ~~Describe \code{VesselConfig} here~~
}
  \item{CovConfig}{
%%     ~~Describe \code{CovConfig} here~~
}
  \item{Aniso}{
%%     ~~Describe \code{Aniso} here~~
}
  \item{ConvergeTol}{
%%     ~~Describe \code{ConvergeTol} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (TmbData, TmbFile, Version, VesselConfig, CovConfig, 
    Aniso, ConvergeTol = 2) 
{
    Parameters = list(ln_H_input = c(0, 0), beta1_t = qlogis(tapply(ifelse(TmbData$b_i > 
        0, 1, 0), INDEX = TmbData$t_i, FUN = mean)), gamma1_j = rep(0, 
        TmbData$n_j), lambda1_k = rep(0, TmbData$n_k), logetaE1 = 0, 
        logetaO1 = 0, logkappa1 = 0, logsigmaV1 = log(1), logsigmaVT1 = log(1), 
        nu1_v = rep(0, TmbData$n_v), nu1_vt = matrix(0, nrow = TmbData$n_v, 
            ncol = TmbData$n_t), Omegainput1_s = rep(0, TmbData$n_s), 
        Epsiloninput1_st = matrix(0, nrow = TmbData$n_s, ncol = TmbData$n_t), 
        beta2_t = log(tapply(ifelse(TmbData$b_i > 0, TmbData$b_i/TmbData$a_i, 
            NA), INDEX = TmbData$t_i, FUN = mean, na.rm = TRUE)), 
        gamma2_j = rep(0, TmbData$n_j), lambda2_k = rep(0, TmbData$n_k), 
        logetaE2 = 0, logetaO2 = 0, logkappa2 = 0, logsigmaV2 = log(1), 
        logsigmaVT2 = log(1), logSigmaM = c(log(5), qlogis(0.8), 
            log(2), log(5)), nu2_v = rep(0, TmbData$n_v), nu2_vt = matrix(0, 
            nrow = TmbData$n_v, ncol = TmbData$n_t), Omegainput2_s = rep(0, 
            TmbData$n_s), Epsiloninput2_st = matrix(0, nrow = TmbData$n_s, 
            ncol = TmbData$n_t))
    Random = c("Epsiloninput1_st", "Omegainput1_s", "Epsiloninput2_st", 
        "Omegainput2_s", "nu1_v", "nu2_v", "nu1_vt", "nu2_vt")
    Map = Make_Map(TmbData = TmbData, VesselConfig = VesselConfig, 
        CovConfig = CovConfig, Aniso = Aniso)
    dyn.load(paste0(TmbFile, "/", dynlib(Version)))
    if (any(FieldConfig != 0) | any(VesselConfig != 0)) {
        Obj <- MakeADFun(data = TmbData, parameters = Parameters, 
            random = Random, hessian = FALSE, map = Map, inner.method = "newton")
    }
    else {
        Obj <- MakeADFun(data = TmbData, parameters = Parameters, 
            hessian = FALSE, map = Map)
    }
    Obj$control <- list(trace = 1, parscale = 1, REPORT = 1, 
        reltol = 1e-12, maxit = 100)
    Lower = rep(-50, length(Obj$par))
    Lower[grep("logsigmaV", names(Obj$par))] = log(0.01)
    Upper = rep(50, length(Obj$par))
    Upper[grep("logtau", names(Obj$par))] = 10
    Upper[grep("logeta", names(Obj$par))] = log(1/(0.01 * sqrt(4 * 
        pi)))
    Upper[grep("SigmaM", names(Obj$par))] = 10
    if ("gamma1" \%in\% names(Obj$par)) {
        Lower[grep("gamma1", names(Obj$par))] = -20
        Upper[grep("gamma1", names(Obj$par))] = 20
    }
    if ("gamma2" \%in\% names(Obj$par)) {
        Lower[grep("gamma2", names(Obj$par))] = -20
        Upper[grep("gamma2", names(Obj$par))] = 20
    }
    if ("lambda1" \%in\% names(Obj$par)) {
        Lower[grep("lambda1", names(Obj$par))] = -20
        Upper[grep("lambda1", names(Obj$par))] = 20
    }
    if ("lambda2" \%in\% names(Obj$par)) {
        Lower[grep("lambda2", names(Obj$par))] = -20
        Upper[grep("lambda2", names(Obj$par))] = 20
    }
    Obj$env$inner.control$step.tol <- c(1e-08, 1e-12, 1e-15)[ConvergeTol]
    Obj$env$inner.control$tol10 <- c(1e-06, 1e-08, 1e-12)[ConvergeTol]
    Obj$env$inner.control$grad.tol <- c(1e-08, 1e-12, 1e-15)[ConvergeTol]
    Return = list(Obj = Obj, Upper = Upper, Lower = Lower)
    return(Return)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
