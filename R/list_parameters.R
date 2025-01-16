

#' List fixed and random effects
#'
#' \code{list_parameters} lists all fixed and random effects
#'
#' @param Obj Compiled TMB object
#' @return Return Tagged-list of fixed and random effects (returned invisibly)

#' @export
list_parameters = function( Obj, verbose=TRUE ){
  Return = list()
  Table = data.frame()
  if( length(Obj$env$random)>0 ){
    Return[["Fixed_effects"]] = names(Obj$env$last.par[-Obj$env$random])
    Return[["Random_effects"]] = names(Obj$env$last.par[Obj$env$random])
    Table = data.frame( "Coefficient_name"=names(table(Return[["Fixed_effects"]])),
                        "Number_of_coefficients"=as.numeric(table(Return[["Fixed_effects"]])),
                        "Type"="Fixed")
    Table = rbind( Table,
                   data.frame("Coefficient_name"=names(table(Return[["Random_effects"]])),
                              "Number_of_coefficients"=as.numeric(table(Return[["Random_effects"]])),
                              "Type"="Random"))
  }else{
    Return[["Fixed_effects"]] = names(Obj$env$last.par)
    Table = data.frame( "Coefficient_name"=names(table(Return[["Fixed_effects"]])),
                        "Number_of_coefficients"=as.numeric(table(Return[["Fixed_effects"]])),
                        "Type"="Fixed")
  }
  if( verbose==TRUE ){
    message("List of estimated fixed and random effects:")
    print(Table)
  }
  return( invisible(Table) )
}
