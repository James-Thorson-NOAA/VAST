

#' Amend output from VAST for user convenience
#'
#' \code{amend_output} add labels, units, and performs logical operations,
#' e.g., adds zeros as needed, to simplify user and downstream interpretation.
#'
#' @export
amend_output <-
function( fit = NULL,
          TmbData = fit$data_list,
          Report = fit$Report,
          extrapolation_list = fit$extrapolation_list,
          Map = fit$tmb_list$Map,
          Sdreport = fit$parameter_estimates$SD,
          year_labels = fit$year_labels,
          category_names = fit$category_names,
          strata_names = fit$strata_names ){

  # Local functions
  add_dimnames = function( Report, report_names, dimnames ){
    for( i in seq_along(report_names) ){
      if( report_names[i] %in% names(Report) ){
        dimnames(Report[[report_names[i]]]) = dimnames
      }
    }
    return(Report)
  }

  # Defaults
  if( "treat_nonencounter_as_zero" %in% names(TmbData$Options_list$Options) ){
    treat_missing_as_zero = TmbData$Options_list$Options["treat_nonencounter_as_zero"]
  }else{
    treat_missing_as_zero = FALSE
  }

  # Local function
  process_labels = function( labels, prefix, length ){
    if(is.null(labels)){
      labels = paste0( prefix, "_", seq_len(length) )
    }else{
      if(length(labels)!=length) stop("Check labels")
    }
    return(labels)
  }

  # Fill in missing inputs
  if( "D_xt" %in% names(Report)){
    # SpatialDeltaGLMM
    year_labels = process_labels( year_labels, "Time", ncol(Report$D_xt) )
    category_names = "singlespecies"
  }
  if( "D_xct" %in% names(Report)){
    # VAST Version < 2.0.0
    year_labels = process_labels( year_labels, "Time", dim(Report$D_xct)[3] )
    category_names = process_labels( category_names, "Category", dim(Report$D_xct)[2] )
  }
  if( "D_xcy" %in% names(Report)){
    # VAST Version >= 2.0.0
    year_labels = process_labels( year_labels, "Time", dim(Report$D_xcy)[3] )
    category_names = process_labels( category_names, "Category_", dim(Report$D_xcy)[2] )
  }
  if( "D_gcy" %in% names(Report)){
    # VAST Version 8.0.0 through 9.3.0
    year_labels = process_labels( year_labels, "Time", dim(Report$D_gcy)[3] )
    category_names = process_labels( category_names, "Category", dim(Report$D_gcy)[2] )
  }
  if( "D_gct" %in% names(Report)){
    # VAST Version >= 9.4.0
    year_labels = process_labels( year_labels, "Time", dim(Report$D_gct)[3] )
    category_names = process_labels( category_names, "Category", dim(Report$D_gct)[2] )
  }
  if("dhat_ktp" %in% names(Report)){
    # MIST Version <= 14
    year_labels = process_labels( year_labels, "Time", dim(Report$dhat_ktp)[2] )
    category_names = process_labels( category_names, "Category", dim(Report$dhat_ktp)[3] )
  }
  if("dpred_ktp" %in% names(Report)){
    # MIST Version >= 15
    year_labels = process_labels( year_labels, "Time", dim(Report$dpred_ktp)[2] )
    category_names = process_labels( category_names, "Category", dim(Report$dpred_ktp)[3] )
  }
  if("Index_ctl" %in% names(Report)){
    strata_names = process_labels( strata_names, "Stratum", dim(Report$Index_ctl)[3] )
  }

  # Determine year-category pairs with no data
  if( "metadata_ctz" %in% names(TmbData$Options_list) ){
    Num_gct = rep(1,TmbData$n_g) %o% abind::adrop(TmbData$Options_list$metadata_ctz[,,'num_notna',drop=FALSE], drop=3)
    Num_ctl = abind::adrop(TmbData$Options_list$metadata_ctz[,,'num_notna',drop=FALSE], drop=3) %o% rep(1,TmbData$n_l)
    Num_ctm = abind::adrop(TmbData$Options_list$metadata_ctz[,,'num_notna',drop=FALSE], drop=3) %o% rep(1,TmbData$n_m)
    if( treat_missing_as_zero==TRUE ){
      # if treat_missing_as_zero==TRUE, then switch density from year-categories with no data to zero
      if("D_gct"%in%names(Report) & all(dim(Report$D_gct)==dim(Num_gct))) Report$D_gct = ifelse(Num_gct==0, 0, Report$D_gct)
      if("D_gcy"%in%names(Report) & all(dim(Report$D_gct)==dim(Num_gct))) Report$D_gcy = ifelse(Num_gct==0, 0, Report$D_gcy)
      if("Index_ctl"%in%names(Report)) Report$Index_ctl = ifelse(Num_ctl==0, 0, Report$Index_ctl)
      if("Index_cyl"%in%names(Report)) Report$Index_cyl = ifelse(Num_ctl==0, 0, Report$Index_cyl)
    }else{
      # If some intercepts are mapped off, then switch density from year-categories with no data to NA
      if( any(is.na(Map$beta2_ft)) | any(is.na(Map$beta2_ft)) ){
        if("D_gct"%in%names(Report) & all(dim(Report$D_gct)==dim(Num_gct))) Report$D_gct = ifelse(Num_gct==0, NA, Report$D_gct)
        if("D_gcy"%in%names(Report) & all(dim(Report$D_gct)==dim(Num_gct))) Report$D_gcy = ifelse(Num_gct==0, NA, Report$D_gcy)
        if("Index_ctl"%in%names(Report)) Report$Index_ctl = ifelse(Num_ctl==0, NA, Report$Index_ctl)
        if("Index_cyl"%in%names(Report)) Report$Index_cyl = ifelse(Num_ctl==0, NA, Report$Index_cyl)
      }
    }
  }

  # No need to map off spatial statistics mean_Z_ctm or effective_area_ctl in any years
  # These are valid when betas are mapped off, or even when epsilons are mapped off given the potential covariates
  #if( any(is.na(Map$beta2_ft)) | any(is.na(Map$beta2_ft)) ){
  #  if("mean_Z_ctm" %in% names(Report)) Report$mean_Z_ctm = ifelse(Num_ctm==0, NA, Report$mean_Z_ctm)
  #  if("effective_area_ctl" %in% names(Report)) Report$effective_area_ctl = ifelse(Num_ctl==0, NA, Report$effective_area_ctl)
  #}

  # Add labels for all variables plotted using `plot_maps`
  Report = add_dimnames( Report = Report,
                         report_names = c("P1_gct","P2_gct","R1_gct","R2_gct","D_gct","Epsilon1_gct","Epsilon2_gct","eta1_gct","eta2_gct"),
                         dimnames = list("Site"=seq_len(TmbData$n_g), "Category"=category_names, "Time"=year_labels) )
  Report = add_dimnames( Report = Report,
                         report_names = c("Omega1_gc","Omega2_gc"),
                         dimnames = list("Site"=seq_len(TmbData$n_g), "Category"=category_names) )
  Report = add_dimnames( Report = Report,
                         report_names = "Xi1_gcp",
                         dimnames = list("Site"=seq_len(TmbData$n_g), "Category"=category_names, "Covariate"=colnames(TmbData$X1_ip)) )
  Report = add_dimnames( Report = Report,
                         report_names = "Xi2_gcp",
                         dimnames = list("Site"=seq_len(TmbData$n_g), "Category"=category_names, "Covariate"=colnames(TmbData$X2_ip)) )
  Report = add_dimnames( Report = Report,
                         report_names = "Phi1_gk",
                         dimnames = list("Site"=seq_len(TmbData$n_g), "Covariate"=colnames(TmbData$Q1_ik)) )
  Report = add_dimnames( Report = Report,
                         report_names = "Phi2_gk",
                         dimnames = list("Site"=seq_len(TmbData$n_g), "Covariate"=colnames(TmbData$Q2_ik)) )
  Report = add_dimnames( Report = Report,
                         report_names = c("L_omega1_cf","L_omega2_cf","L_beta1_cf","L_beta2_cf","L_epsilon1_cf","L_epsilon2_cf"),
                         dimnames = list("Category"=category_names, NULL) )
  Report = add_dimnames( Report = Report,
                         report_names = c("Ltime_epsilon1_tf","Ltime_epsilon2_tf"),
                         dimnames = list("Time"=year_labels, NULL) )
  Report = add_dimnames( Report = Report,
                         report_names = c("beta1_tc","beta2_tc"),
                         dimnames = list("Time"=year_labels, "Category"=category_names) )

  # Add labels for other useful variables
  Report = add_dimnames( Report = Report,
                         report_names = c("Index_ctl","effective_area_ctl","mean_D_ctl","Bratio_ctl"),
                         dimnames = list("Category"=category_names, "Time"=year_labels, "Stratum"=strata_names) )
  Report = add_dimnames( Report = Report,
                         report_names = c("Fratio_ct"),
                         dimnames = list("Category"=category_names, "Time"=year_labels) )
  Report = add_dimnames( Report = Report,
                         report_names = c("Index_gctl"),
                         dimnames = list("Site"=seq_len(TmbData$n_g), "Category"=category_names, "Time"=year_labels, "Stratum"=strata_names) )
  Report = add_dimnames( Report = Report,
                         report_names = "mean_Z_ctm",
                         dimnames = list("Category"=category_names, "Time"=year_labels, "Spatial_axis"=colnames(TmbData$Z_gm)) )

  # Modify Sdreport
  if( !is.null(Sdreport) ){
    # See plot_biomass_index for how to efficiently use and combine SEs
  }
  #deparse_unit = function(x){
  #  out = units::deparse_unit(x)
  #  if(out=="") out = "1"
  #  return(out)
  #}

  # Add units
  if("Index_ctl" %in% names(Report)){
    # Convert and simplify units ... rescales Report so might be confusing
      converted_units = as_units(1,units(TmbData$b_i)) / as_units(1,units(TmbData$a_i)) * as_units(1,units(extrapolation_list$Area_km2[1]))
      Report$Index_ctl = Report$Index_ctl * converted_units
    # Label units without simplifying them ... clunky plots AND giving error in units version 0.8.5
    #units(Report$Index_ctl) = paste0( deparse_unit(TmbData$b_i), " / ", deparse_unit(TmbData$a_i), " * ", deparse_unit(extrapolation_list$Area_km2[1]) )
    #units(Report$Index_ctl) = deparse_unit(TmbData$b_i)
    #if(deparse_unit(TmbData$a_i) != "") units(Report$Index_ctl) = paste( deparse_unit(Report$Index_ctl), "/", deparse_unit(TmbData$a_i) )
    #if(deparse_unit(extrapolation_list$Area_km2[1]) != "") units(Report$Index_ctl) = paste( deparse_unit(Report$Index_ctl), "*", deparse_unit(extrapolation_list$Area_km2[1]) )
  }
  #if("Index_ctl" %in% names(Report)) units(Report$Index_ctl) = units(TmbData$b_i / TmbData$a_i * extrapolation_list$Area_km2[1])
  if("Fratio_ct" %in% names(Report)) units(Report$Fratio_ct) = units(unitless)
  if("Bratio_ctl" %in% names(Report)) units(Report$Bratio_ctl) = units(unitless)
  if("D_gct" %in% names(Report)){
    # Convert and simplify units ... rescales Report so might be confusing
      converted_units = as_units(1,units(TmbData$b_i)) / as_units(1,units(TmbData$a_i))
      Report$D_gct = Report$D_gct * converted_units
    # Label units without simplifying them ... clunky plots AND giving error in units version 0.8.5
    #units(Report$D_gct) = paste0( deparse_unit(TmbData$b_i), " / ", deparse_unit(TmbData$a_i) )
  }

  # Add units for COG, see: https://github.com/r-quantities/units/issues/291
  # In case of re-running amend_output on objects with existing units
  if("mean_Z_ctm" %in% names(Report)) units(Report$mean_Z_ctm) = as_units("km")
  if("mean_D_ctl" %in% names(Report)) units(Report$mean_D_ctl) = units(TmbData$b_i)
  if("effective_area_ctl" %in% names(Report)) units(Report$effective_area_ctl) = as_units("km^2") #  = paste0( sf::st_crs( extrapolation_list$projargs )$units, "^2" )

  # Check for bad entries
  return( Report )
}

