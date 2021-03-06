% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effects.R
\name{Effect.fit_model}
\alias{Effect.fit_model}
\title{Adapts package \code{effects}}
\usage{
Effect.fit_model(focal.predictors, mod, which_formula = "X1", ...)
}
\arguments{
\item{focal.predictors}{a character vector of one or more predictors in the model in any order.}

\item{mod}{an object of the appropriate class.  If no method exists for that class, \code{Effect.default} will be called. }

\item{which_formula}{which formula to use e.g., \code{"X1"}}

\item{...}{arguments to be passed down.}
}
\description{
Calculate effects for plotting
}
