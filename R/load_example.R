
#' Load an example data set
#'
#' \code{load_example} loads a catch and effort data set from multiple surveys worldwide
#'
#' Current options for \code{data_set} include:
#' \itemize{
#' \item \code{"Chatham_rise_hake"},
#' \item \code{"Iceland_cod"},
#' \item \code{"WCGBTS_canary"},
#' \item \code{"GSL_american_plaice"},
#' \item \code{"BC_pacific_cod"},
#' \item \code{"EBS_pollock"},
#' \item \code{"GOA_Pcod"},
#' \item \code{"GOA_pollock"},
#' \item \code{"GB_spring_haddock"} (with thanks to L. Brooks for assembling the data),
#' \item \code{"GB_fall_haddock"} (with thanks to L. Brooks for assembling the data),
#' \item \code{"SAWC_jacopever"} (with thanks to H. Winker for assembling the data),
#' \item \code{"Aleutian_islands_POP"},
#' \item \code{"GOA_arrowtooth_condition_and_density"} (with thanks to A. Gruss for assembling the data),
#' \item \code{"condition_and_density"} (with thanks to A. Gruss for assembling the data),
#' \item \code{"multimodal_red_snapper"} (with thanks to A. Gruss for assembling the data),
#' \item \code{"lingcod_comp_expansion"} (with thanks to M. Haltuch for assembling the data),
#' \item \code{"covariate_example"} (with thanks to D. McGowan for assembling the data),
#' \item \code{"PESC_example_red_grouper"} (with thanks to A. Gruss for assembling the data),
#' \item \code{"ordination"}, and
#' \item \code{"NWA_yellowtail_seasons"} (with thanks to L. Brooks, C. Adams, and C. Legault for assembling the data).
#' }
#' These examples are used to highlight different functionality for spatio-temporal analysis,
#' as well as during integrated testing to check whether updates are backwards compatible for these examples.
#'
#' @param data_set data set to load
#'
#' @export
load_example = function( data_set="EBS_pollock" ){

  # Determine region for a given data set
  region = switch( tolower(data_set),
                   "chatham_rise_hake" = tolower("New_Zealand"),
                   "wcgbts_canary" = tolower("California_current"),
                   "gsl_american_plaice" = tolower("Gulf_of_St_Lawrence"),
                   "bc_pacific_cod" = tolower("British_Columbia"),
                   "ebs_pollock" = tolower("Eastern_Bering_Sea"),
                   "goa_Pcod" = tolower("Gulf_of_Alaska"),
                   "goa_pollock" = tolower("Gulf_of_Alaska"),
                   "gb_spring_haddock" = tolower("Northwest_Atlantic"),
                   "gb_fall_haddock" = tolower("Northwest_Atlantic"),
                   "sawc_jacopever" = tolower("South_Africa"),
                   "aleutian_islands_pop" = tolower("Aleutian_Islands"),
                   "condition_and_density" = tolower("Eastern_Bering_Sea"),
                   "goa_arrowtooth_condition_and_density" = tolower("Gulf_of_Alaska"),
                   "multimodal_red_snapper" = tolower("Gulf_of_Mexico"),
                   "lingcod_comp_expansion" = tolower("California_current"),
                   "ordination" = tolower("Eastern_Bering_Sea"),
                   "five_species_ordination" = tolower("Eastern_Bering_Sea"),
                   "covariate_example" = tolower("Gulf_of_Alaska"),
                   "goa_pcod_covariate_example" = tolower("Gulf_of_Alaska"),
                   "goa_mice_example" = tolower("Gulf_of_Alaska"),
                   "PESC_example_red_grouper" = tolower("User"),
                   "nwa_yellowtail_seasons" = tolower("Northwest_Atlantic"),
                   tolower("Other") )

  # Initialize all objects as NULL
  Predator_biomass_cath_rate_data = Stomach_content_data =
    input_grid = covariate_data = F_ct = X_xtp = X_gtp = X_itp = Q_ik = NULL

  if( tolower(data_set) %in% tolower("WCGBTS_canary") ){
    data( WCGBTS_Canary_example, package="VAST" )
    Year = as.numeric(sapply(WCGBTS_Canary_example[,'PROJECT_CYCLE'], FUN=function(Char){strsplit(as.character(Char)," ")[[1]][2]}))
    sampling_data = data.frame( "Catch_KG"=WCGBTS_Canary_example[,'HAUL_WT_KG'], "Year"=Year, "Vessel"=paste(WCGBTS_Canary_example[,"VESSEL"],Year,sep="_"), "AreaSwept_km2"=WCGBTS_Canary_example[,"AREA_SWEPT_HA"]/1e2, "Lat"=WCGBTS_Canary_example[,'BEST_LAT_DD'], "Lon"=WCGBTS_Canary_example[,'BEST_LON_DD'], "Pass"=WCGBTS_Canary_example[,'PASS']-1.5)
    strata.limits = data.frame('STRATA'="All_areas")
  }
  if( tolower(data_set) %in% tolower("BC_pacific_cod") ){
    data( BC_pacific_cod_example, package="VAST" )
    sampling_data = data.frame( "Catch_KG"=BC_pacific_cod_example[,'PCOD_WEIGHT'], "Year"=BC_pacific_cod_example[,'Year'], "Vessel"="missing", "AreaSwept_km2"=BC_pacific_cod_example[,'TOW.LENGTH..KM.']/100, "Lat"=BC_pacific_cod_example[,'LAT'], "Lon"=BC_pacific_cod_example[,'LON'], "Pass"=0)
    strata.limits = data.frame('STRATA'="All_areas")
  }
  if( tolower(data_set) %in% tolower("GSL_american_plaice") ){
    data( GSL_american_plaice, package="VAST" )
    Print_Message( "GSL_american_plaice" )
    sampling_data = data.frame( "Year"=GSL_american_plaice[,'year'], "Lat"=GSL_american_plaice[,'latitude'], "Lon"=GSL_american_plaice[,'longitude'], "Vessel"="missing", "AreaSwept_km2"=GSL_american_plaice[,'swept'], "Catch_KG"=GSL_american_plaice[,'biomass']*GSL_american_plaice[,'vstd'] )
    strata.limits = data.frame('STRATA'="All_areas")
  }
  if( tolower(data_set) %in% tolower("EBS_pollock") ){
    data( EBS_pollock_data, package="VAST" )
    covariate_data = EBS_pollock_data$covariate_data
    sampling_data = EBS_pollock_data$sampling_data
    sampling_data = data.frame( "Catch_KG"=sampling_data[,'catch'], "Year"=sampling_data[,'year'], "Vessel"="missing", "AreaSwept_km2"=0.01, "Lat"=sampling_data[,'lat'], "Lon"=sampling_data[,'long'], "Pass"=0)
    strata.limits = data.frame('STRATA'="All_areas")
  }
  if( tolower(data_set) %in% tolower("GOA_Pcod") ){
    data( GOA_pacific_cod , package="VAST")
    sampling_data = data.frame( "Catch_KG"=GOA_pacific_cod[,'catch'], "Year"=GOA_pacific_cod[,'year'], "Vessel"="missing", "AreaSwept_km2"=0.01, "Lat"=GOA_pacific_cod[,'lat'], "Lon"=GOA_pacific_cod[,'lon'], "Pass"=0)
    strata.limits = data.frame('STRATA'="All_areas")
  }
  if( tolower(data_set) %in% tolower("GOA_pollock") ){
    data( GOA_walleye_pollock, package="VAST" )
    sampling_data = data.frame( "Catch_KG"=GOA_walleye_pollock[,'catch'], "Year"=GOA_walleye_pollock[,'year'], "Vessel"="missing", "AreaSwept_km2"=0.01, "Lat"=GOA_walleye_pollock[,'lat'], "Lon"=GOA_walleye_pollock[,'lon'], "Pass"=0)
    strata.limits = data.frame('STRATA'="All_areas")
  }
  if( tolower(data_set) %in% tolower("Aleutian_islands_POP") ){
    data( AI_pacific_ocean_perch, package="VAST" )
    sampling_data = data.frame( "Catch_KG"=AI_pacific_ocean_perch[,'cpue..kg.km.2.'], "Year"=AI_pacific_ocean_perch[,'year'], "Vessel"="missing", "AreaSwept_km2"=1, "Lat"=AI_pacific_ocean_perch[,'start.latitude'], "Lon"=AI_pacific_ocean_perch[,'start.longitude'], "Pass"=0)
    strata.limits = data.frame('STRATA'="All_areas")
  }
  if( tolower(data_set) %in% tolower("GB_spring_haddock") ){
    data( georges_bank_haddock_spring, package="VAST" )
    Print_Message( "GB_haddock" )
    sampling_data = data.frame( "Catch_KG"=georges_bank_haddock_spring[,'CATCH_WT_CAL'], "Year"=georges_bank_haddock_spring[,'YEAR'], "Vessel"="missing", "AreaSwept_km2"=0.0112*1.852^2, "Lat"=georges_bank_haddock_spring[,'LATITUDE'], "Lon"=georges_bank_haddock_spring[,'LONGITUDE'])
    # For NEFSC indices, strata must be specified as a named list of area codes
    strata.limits = list( 'Georges_Bank'=c(1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200, 1210, 1220, 1230, 1240, 1250, 1290, 1300) )
  }
  if( tolower(data_set) %in% tolower("GB_fall_haddock") ){
    data( georges_bank_haddock_fall, package="VAST" )
    Print_Message( "GB_haddock" )
    sampling_data = data.frame( "Catch_KG"=georges_bank_haddock_fall[,'CATCH_WT_CAL'], "Year"=georges_bank_haddock_fall[,'YEAR'], "Vessel"="missing", "AreaSwept_km2"=0.0112*1.852^2, "Lat"=georges_bank_haddock_fall[,'LATITUDE'], "Lon"=georges_bank_haddock_fall[,'LONGITUDE'])
    # For NEFSC indices, strata must be specified as a named list of area codes
    strata.limits = list( 'Georges_Bank'=c(1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200, 1210, 1220, 1230, 1240, 1250, 1290, 1300) )
  }
  if( tolower(data_set) %in% tolower("SAWC_jacopever") ){
    data( south_africa_westcoast_jacopever, package="VAST" )
    sampling_data = data.frame( "Catch_KG"=south_africa_westcoast_jacopever[,'HELDAC'], "Year"=south_africa_westcoast_jacopever[,'Year'], "Vessel"="missing", "AreaSwept_km2"=south_africa_westcoast_jacopever[,'area_swept_nm2']*1.852^2, "Lat"=south_africa_westcoast_jacopever[,'cen_lat'], "Lon"=south_africa_westcoast_jacopever[,'cen_long'])
    strata.limits = data.frame('STRATA'="All_areas")
  }
  if( tolower(data_set) %in% tolower("Iceland_cod") ){
    # WARNING:  This data set has not undergone much evaluation for spatio-temporal analysis
    data( iceland_cod, package="VAST" )
    sampling_data = data.frame( "Catch_KG"=iceland_cod[,'Catch_b'], "Year"=iceland_cod[,'year'], "Vessel"=1, "AreaSwept_km2"=iceland_cod[,'towlength'], "Lat"=iceland_cod[,'lat1'], "Lon"=iceland_cod[,'lon1'])
    strata.limits = data.frame('STRATA'="All_areas")
  }
  if( tolower(data_set) %in% tolower("Chatham_rise_hake") ){
    data( chatham_rise_hake, package="VAST" )
    sampling_data = data.frame( "Catch_KG"=chatham_rise_hake[,'Hake_kg_per_km2'], "Year"=chatham_rise_hake[,'Year'], "Vessel"=1, "AreaSwept_km2"=1, "Lat"=chatham_rise_hake[,'Lat'], "Lon"=chatham_rise_hake[,'Lon'])
    strata.limits = data.frame('STRATA'="All_areas")
  }
  if( tolower(data_set) %in% tolower("condition_and_density") ){
    data( condition_and_density_example, package="VAST" )
    sampling_data = condition_and_density_example
    strata.limits = data.frame('STRATA'="All_areas")
  }
  if( tolower(data_set) %in% tolower("GOA_arrowtooth_condition_and_density") ){
    data( GOA_arrowtooth_condition_and_density, package="VAST" )
    sampling_data = GOA_arrowtooth_condition_and_density
    strata.limits = data.frame('STRATA'="All_areas")
  }
  if( tolower(data_set) %in% tolower("multimodal_red_snapper") ){
    data( multimodal_red_snapper_example, package="VAST" )
    sampling_data = multimodal_red_snapper_example
    strata.limits = data.frame('STRATA'="All_areas")
  }
  if( tolower(data_set) %in% tolower("Lingcod_comp_expansion") ){
    data( comp_expansion_example, package="VAST" )
    sampling_data = comp_expansion_example
    strata.limits = data.frame( 'STRATA' = "ORWA", 'north_border' = 49.0, 'south_border' = 42.0 )
  }
  if( tolower(data_set) %in% tolower(c("ordination","five_species_ordination")) ){
    data( five_species_ordination_example, package="VAST" )
    sampling_data = five_species_ordination_example
    strata.limits = data.frame('STRATA'="All_areas")
  }
  if( tolower(data_set) %in% tolower(c("covariate_example","GOA_pcod_covariate_example")) ){
    data( GOA_pcod_covariate_example, package="VAST" )
    sampling_data = GOA_Pcod_covariate_example$sampling_data
    strata.limits = data.frame('STRATA'="All_areas")
    X_xtp = GOA_Pcod_covariate_example$X_xtp
    covariate_data = GOA_Pcod_covariate_example$covariate_data
  }
  if( tolower(data_set) %in% tolower(c("GOA_MICE_example")) ){
    data( goa_mice_example, package="VAST" )
    sampling_data = goa_mice_example$Data_Geostat
    F_tc = goa_mice_example$F_tc
    sampling_data = sampling_data[ which(sampling_data[,'Year'] %in% unique(F_tc[,1])), ]
    strata.limits = data.frame('STRATA'="All_areas")
    F_ct = t( F_tc[which(F_tc[,'X'] %in% min(sampling_data[,'Year']):max(sampling_data[,'Year'])),-1] )
    colnames(F_ct) = min(sampling_data[,'Year']):max(sampling_data[,'Year'])
  }
  if( tolower(data_set) %in% tolower("PESC_example_red_grouper") ){
    data( PESC_example_red_grouper, package="VAST" )
    Predator_biomass_cath_rate_data = example$Predator_biomass_cath_rate_data
    Stomach_content_data = example$Stomach_content_data
    region = example$Region
    strata.limits = example$strata.limits
    input_grid = example$input_grid
  }
  if( tolower(data_set) %in% tolower("NWA_yellowtail_seasons") ){
    data( NWA_yellowtail_seasons, package="VAST" )
    sampling_data = NWA_yellowtail_seasons
    strata.limits = list( c(1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200, 1210) )
  }
  #sampling_data = na.omit( sampling_data )

  Return = list("sampling_data"=sampling_data, "Region"=region, "strata.limits"=strata.limits)
  if( !is.null(X_xtp)) Return[["X_xtp"]] = X_xtp
  if( !is.null(X_gtp)) Return[["X_gtp"]] = X_gtp
  if( !is.null(X_itp)) Return[["X_itp"]] = X_itp
  if( !is.null(Q_ik)) Return[["Q_ik"]] = Q_ik
  if( !is.null(F_ct)) Return[["F_ct"]] = F_ct
  if( !is.null(covariate_data)) Return[["covariate_data"]] = covariate_data
  if( !is.null(input_grid)) Return[["input_grid"]] = input_grid
  if( !is.null(Predator_biomass_cath_rate_data)) Return[["Predator_biomass_cath_rate_data"]] = Predator_biomass_cath_rate_data
  if( !is.null(Stomach_content_data)) Return[["Stomach_content_data"]] = Stomach_content_data
  
  # return stuff
  return(Return)
}
