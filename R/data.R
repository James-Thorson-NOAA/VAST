
#' Data to demonstrate and test stream-network options
#'
#' Data sufficient to demonstrate stream-network options in VAST, using \code{Method="Stream_network"},
#'  and also helpful to do integrated testing of this feature
#'
#' \itemize{
#'   \item observations data-frame of biological sampling data
#'   \item network data-frame defining stream-network connectivity, used to define Ornstein-Uhlenbeck covariance function
#'   \item habitat data-frame containing density covariates
#' }
#'
#' @name stream_network_eel_example
#' @docType data
#' @author Merrill Rudd
#' @usage data(stream_network_eel_example)
#' @keywords data
NULL

#' Data to demonstrate and test covariate effects
#'
#' Data sufficient to demonstrate spatially varying coefficient models
#'
#' \itemize{
#'   \item sampling_data data-frame of biological sampling data and associated covariate measurements
#'   \item Region region for model demo
#'   \item strata.limits user-specified stratification of results
#'   \item X_xtp 3D array of covariates, specified at knots given the use of \code{fine_scale=FALSE}
#' }
#'
#' @name GOA_pcod_covariate_example
#' @docType data
#' @author Dave McGowan
#' @usage data(GOA_pcod_covariate_example)
#' @keywords data
NULL

#' Data to demonstrate and test decadal-trends model
#'
#' Data sufficient to demonstrate how spatially varying coefficients (SVCs)
#'      are used to represent decadal trends.  Specifically using public data
#'      from the Alaska Fisheries Science Center bottom trawl survey
#'      (processed to add true zeros) for arrowtooth flounder in the eastern Bering Sea,
#'      accessed Sept. 8, 2022.
#'
#' \itemize{
#'   \item sampling_data data-frame of tow-specific biomass
#' }
#'
#' @name arrowtooth_trends
#' @docType data
#' @author James Thorson
#' @usage data(arrowtooth_trends)
#' @keywords data
NULL

#' Data to demonstrate and test trait-based multispeces SDM
#'
#' Data sufficient to demonstrate how spatially varying coefficients (SVCs)
#'      are used to represent trait effects in a joint (multispecies)
#'      species distribution model.  Specifically using public data
#'      from the Alaska Fisheries Science Center bottom trawl survey
#'      (processed to add true zeros) for numerically abundance fishes in the
#'      Gulf of Alaska, accessed Sept. 8, 2022.  Also including trait estimates
#'      from an initial (unpublished) release of FishLife that
#'      includes morphometric and habitat traits.
#'
#' \itemize{
#'   \item sampling_data data-frame of tow-specific biomass
#'   \item trait_data data-frame of trait-values (for continuous traits) and/or probabilities (for categorical traits)
#' }
#'
#' @name goa_traits
#' @docType data
#' @author James Thorson
#' @usage data(goa_traits)
#' @keywords data
NULL

#' Data to demonstrate and test cohort effects in age-structured SDM
#'
#' Data sufficient to demonstrate how spatially varying coefficients (SVCs)
#'      are used to represent cohort effects in a joint (age-structured)
#'      species distribution model.  Specifically using processed data
#'      from the Alaska Fisheries Science Center bottom trawl survey
#'      that includes a density-dependent correction (based on subsampled acoustics)
#'      and allocated to ages using an age-length-key (ALK), where both were
#'      conducted by Stan Kotwicki in support
#'
#' \itemize{
#'   \item sampling_data data-frame of tow-specific biomass
#' }
#'
#' @references
#' Stevenson, D. E., Kotwicki, S., Thorson, J. T., Correa, G. M., and Buckley, T. 2022.
#' The influence of age and cohort on the distribution of walleye pollock (Gadus chalcogrammus)
#' in the eastern Bering Sea. Canadian Journal of Fisheries and Aquatic Sciences.
#' <https://doi.org/10.1139/cjfas-2021-030>
#'
#' @name pollock_cohorts
#' @docType data
#' @author Stan Kotwicki
#' @usage data(pollock_cohorts)
#' @keywords data
NULL
