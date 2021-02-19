fit_model<- function(settings, Lat_i, Lon_i, t_i, b_i, a_i, c_iz = rep(0, length(b_i)), v_i = rep(0, length(b_i)), working_dir = paste0(getwd(), "/"), X1config_cp = NULL, X2config_cp = NULL, covariate_data, X1_formula = ~0, X2_formula = ~0, Q1config_k = NULL, Q2config_k = NULL, catchability_data, Q1_formula = ~0, Q2_formula = ~0, newtonsteps = 1, silent = TRUE, build_model = TRUE, run_model = TRUE, test_fit = TRUE, ...) {
  
  if(FALSE){
    "settings" = settings
    "Lat_i" = samp_dat[, 'Lat']
    "Lon_i" = samp_dat[, 'Lon']
    "t_i" = samp_dat[, 'Year']
    "c_i" = rep(0, nrow(samp_dat))
    "c_iz" = rep(0, length(b_i))
    "v_i" = rep(0, length(b_i))
    "b_i" = samp_dat[, 'Response']
    "a_i" = samp_dat[, 'Swept']
    "X1config_cp" = X1config_cp_use
    "X2config_cp" = X2config_cp_use
    "covariate_data" = cov_dat
    "X1_formula" = formula_use
    "X2_formula" = formula_use
    "Q1config_k" = NULL
    "Q2config_k" = NULL
    "Q1_formula" = ~0
    "Q2_formula" = ~0
    "newtonsteps" = 1
    "getsd" = TRUE
    "getReportCovariance" = TRUE
    "observations_LL" = cbind("Lat" = samp_dat[, 'Lat'], "Lon" = samp_dat[, 'Lon'])
    "maximum_distance_from_sample" = maximum_distance_from_sample
    "grid_dim_km" = grid_dim_km
    "run_model" = FALSE
    "test_fit" = FALSE
    working_dir = RunDir
    "PredTF_i" = samp_dat[, 'Dummy']
    "Use_REML" = Use_REML
    "getJointPrecision" = FALSE
    extra_args = list()
  }
  extra_args = list(...)
  extra_args = c(extra_args, extra_args$extrapolation_args, extra_args$spatial_args, extra_args$optimize_args, extra_args$model_args)
  data_frame = data.frame(Lat_i = Lat_i, Lon_i = Lon_i, a_i = a_i, v_i = v_i, b_i = b_i, t_i = t_i, c_iz = c_iz)
  year_labels = seq(min(t_i), max(t_i))
  years_to_plot = which(year_labels %in% t_i)
  message("\n### Writing output from `fit_model` in directory: ", working_dir)
  dir.create(working_dir, showWarnings = FALSE, recursive = TRUE)
  capture.output(settings, file = file.path(working_dir, "settings.txt"))
  message("\n### Making extrapolation-grid")
  extrapolation_args_default = list(Region = settings$Region, strata.limits = settings$strata.limits, zone = settings$zone, max_cells = settings$max_cells, DirPath = working_dir)
  extrapolation_args_input = combine_lists(input = extra_args, default = extrapolation_args_default, args_to_use = formalArgs(make_extrapolation_info))
  extrapolation_list = do.call(what = make_extrapolation_info, args = extrapolation_args_input)
  message("\n### Making spatial information")
  spatial_args_default = list(grid_size_km = settings$grid_size_km, n_x = settings$n_x, Method = settings$Method, Lon_i = Lon_i, Lat_i = Lat_i, Extrapolation_List = extrapolation_list, DirPath = working_dir, Save_Results = TRUE, fine_scale = settings$fine_scale, knot_method = settings$knot_method)
  spatial_args_input = combine_lists(input = extra_args, default = spatial_args_default, args_to_use = c(formalArgs(make_spatial_info), formalArgs(INLA::inla.mesh.create)))
  spatial_list = do.call(what = make_spatial_info, args = spatial_args_input)
  message("\n### Making data object")
  if (missing(covariate_data)) 
    covariate_data = NULL
  if (missing(catchability_data)) 
    catchability_data = NULL
  data_args_default = list(Version = settings$Version, FieldConfig = settings$FieldConfig, 
                           OverdispersionConfig = settings$OverdispersionConfig, 
                           RhoConfig = settings$RhoConfig, VamConfig = settings$VamConfig, 
                           ObsModel = settings$ObsModel, c_iz = c_iz, b_i = b_i, 
                           a_i = a_i, v_i = v_i, s_i = spatial_list$knot_i - 1, 
                           t_i = t_i, spatial_list = spatial_list, Options = settings$Options, 
                           Aniso = settings$use_anisotropy, X1config_cp = X1config_cp, 
                           X2config_cp = X2config_cp, covariate_data = covariate_data, 
                           X1_formula = X1_formula, X2_formula = X2_formula, Q1config_k = Q1config_k, 
                           Q2config_k = Q2config_k, catchability_data = catchability_data, 
                           Q1_formula = Q1_formula, Q2_formula = Q2_formula)
  data_args_input = combine_lists(input = extra_args, default = data_args_default)
  data_list = do.call(what = make_data, args = data_args_input)
  
  #####
  ## Breaking into make_data!!!
  #####
  Version = settings$Version
  FieldConfig = settings$FieldConfig
  OverdispersionConfig = settings$OverdispersionConfig
  RhoConfig = settings$RhoConfig
  VamConfig = settings$VamConfig
  ObsModel = settings$ObsModel
  ObsModel_ez = c(PosDist = 1, Link = 0) # In function...
  c_iz = c_iz
  b_i = b_i
  a_i = a_i
  v_i = v_i
  s_i = spatial_list$knot_i - 1
  t_i = t_i
  e_i = c_iz[, 1]
  spatial_list = spatial_list
  Options = settings$Options
  Aniso = settings$use_anisotropy
  X1config_cp = X1config_cp
  X2config_cp = X2config_cp
  covariate_data = covariate_data
  X1_formula = X1_formula
  X2_formula = X2_formula
  Q1config_k = Q1config_k
  Q2config_k = Q2config_k
  catchability_data = catchability_data
  Q1_formula = Q1_formula
  Q2_formula = Q2_formula
  alternate_inputs = list(ObsModel = settings$ObsModel)
  
  MeshList = spatial_list[["MeshList"]]
  GridList = spatial_list[["GridList"]]
  Method = spatial_list[["Method"]]
  a_xl = a_gl = spatial_list[["a_gl"]]
  s_i = spatial_list[["knot_i"]] - 1
  
  Options2use = c(SD_site_density = FALSE, SD_site_logdensity = FALSE, 
                  Calculate_Range = FALSE, SD_observation_density = FALSE, 
                  Calculate_effective_area = FALSE, Calculate_Cov_SE = FALSE, 
                  Calculate_Synchrony = FALSE, Calculate_Coherence = FALSE, 
                  Calculate_proportion = FALSE, normalize_GMRF_in_CPP = TRUE, 
                  Calculate_Fratio = FALSE, Estimate_B0 = FALSE, Project_factors = FALSE, 
                  treat_nonencounter_as_zero = FALSE, simulate_random_effects = TRUE, 
                  observation_error_as_CV = TRUE, report_additional_variables = FALSE, 
                  zerosum_penalty = 0)
  
  for (i in seq_along(Options)) {
    if (tolower(names(Options)[i]) %in% tolower(names(Options2use))) {
      Options2use[[match(tolower(names(Options)[i]), tolower(names(Options2use)))]] = Options[[i]]
    }
  }
  
  if (is.vector(FieldConfig) && length(FieldConfig) == 4) {
    FieldConfig = rbind(matrix(FieldConfig, ncol = 2, dimnames = list(c("Omega", "Epsilon"), c("Component_1", "Component_2"))), Beta = c("IID", "IID"))
  }
  if (is.matrix(FieldConfig) & all(dim(FieldConfig) == c(3, 2))) {
    FieldConfig = rbind(FieldConfig, Epsilon_time = c("Identity", "Identity"))
  }
  if (!is.matrix(FieldConfig) || !all(dim(FieldConfig) ==  c(4, 2))) {
    stop("`FieldConfig` has the wrong dimensions in `make_data`")
  }
  dimnames(FieldConfig) = list(c("Omega", "Epsilon", "Beta", "Epsilon_time"), c("Component_1", "Component_2"))
  tprime_i = t_i - min(t_i, na.rm = TRUE)
  
  c_iz = matrix(c_iz, ncol = 1)
  n_t = max(tprime_i, na.rm = TRUE) + 1
  n_c = max(c_iz, na.rm = TRUE) + 1
  n_e = max(e_i) + 1
  n_v = length(unique(v_i))
  n_i = length(b_i)
  n_x = nrow(a_gl)
  n_l = ncol(a_gl)
  n_g = ifelse(is.null(spatial_list), 1, spatial_list$n_g)
  
  # Weird...
  ObsModel_ez = matrix(ObsModel_ez, ncol = 2, nrow = n_e, byrow = TRUE)
  
  Index = list(factor(c_iz[, 1], levels = 0:max(c_iz[, 1])), factor(tprime_i, levels = 0:max(tprime_i)))
  Num_ct = tapply(b_i, INDEX = Index, FUN = function(vec) {
    sum(vec > 0, na.rm = TRUE)
  })
  Num_ct = ifelse(is.na(Num_ct), 0, Num_ct)
  b_i = ifelse(Num_ct[cbind(as.numeric(Index[[1]]), as.numeric(Index[[2]]))] ==  0, NA, b_i)
  
  X_xtp = alternate_inputs[["X_xtp"]]
  X_gctp = alternate_inputs[["X_gctp"]]
  X1_gctp = alternate_inputs[["X1_gctp"]]
  X2_gctp = alternate_inputs[["X2_gctp"]]
  X_itp = alternate_inputs[["X_itp"]]
  X1_itp = alternate_inputs[["X1_itp"]]
  X2_itp = alternate_inputs[["X2_itp"]]
  
  Covariates_created = FALSE
  if (Covariates_created == FALSE) {
    if (!is.null(covariate_data)) {
      Covariates_created = TRUE
      if (FishStatsUtils::convert_version_name(Version) <= FishStatsUtils::convert_version_name("VAST_v9_4_0")) {
        stop("To use separate formula interface for linear predictors, please use version >= CPP 10.0.0")
      }
      covariate_list = make_covariates_aja(formula = X1_formula, 
                                       covariate_data = covariate_data, Year_i = t_i, 
                                       spatial_list = spatial_list, contrasts = contrasts.arg)
      X1_gtp = covariate_list$X_gtp # 2000, 99, 3
      
      # Why is this 34? 
      
      #####
      ## Breaking into make_covariates
      #####
      formula = X1_formula
      covariate_data = covariate_data
      Year_i = t_i
      spatial_list = spatial_list
    
      Year_Set = min(Year_i):max(Year_i)
      covariate_df = covariate_data[which(!is.na(covariate_data[, "Year"])), ]
    
      for (tI in seq_along(Year_Set)) {
        newrows = covariate_data[which(is.na(covariate_data[, "Year"])), ]
        newrows[, "Year"] = rep(Year_Set[tI], nrow(newrows))
        covariate_df = rbind(covariate_df, newrows)
      }
      Model_matrix = model.matrix(update.formula(formula, ~. +  1), data = covariate_df)
      # Alright, so this is clearly where we are having issues. For example, the first observation is spring 1985. In the model_matrix, though, this observation has the intercept, and SeasonSpring. In this intercept, though is DFOSeason AND Year_Cov1985. Working with contrasts...
      Model_matrix = model.matrix(update.formula(formula, ~. + 1), data = covariate_df, contrasts = list(Season=contrasts(covariate_df$Season, contrasts = FALSE), Year_Cov = contrasts(covariate_df$Year_Cov, contrasts = FALSE)))
      
      Columns_to_keep = which(attr(Model_matrix, "assign") != 0)
      coefficient_names = attr(Model_matrix, "dimnames")[[2]][Columns_to_keep]
      X = Model_matrix[, Columns_to_keep, drop = FALSE]
      dimnames(X) = list(NULL, coefficient_names)
      sample_i = data.frame(Year = Year_i, Lat = spatial_list$latlon_i[, 
                                                                       "Lat"], Lon = spatial_list$latlon_i[, "Lon"])
      latlon_g = spatial_list$latlon_g
      X_gtp = array(NA, dim = c(nrow(latlon_g), length(Year_Set), 
                                ncol(X)), dimnames = list(NULL, Year_Set, colnames(X)))
      X_ip = array(NA, dim = c(nrow(sample_i), ncol(X)), dimnames = list(NULL, 
                                                                         colnames(X)))
      for (tI in seq_along(Year_Set)) {
        tmp_covariate_df = covariate_df[which(Year_Set[tI] == 
                                                covariate_df[, "Year"]), , drop = FALSE]
        tmp_X = X[which(Year_Set[tI] == covariate_df[, "Year"]), 
                  , drop = FALSE]
        if (nrow(tmp_covariate_df) == 0) {
          stop("Year ", Year_Set[tI], " not found in `covariate_data` please specify covariate values for all years")
        }
        Which = which(Year_Set[tI] == sample_i[, "Year"])
        if (length(Which) > 0) {
          NN = RANN::nn2(data = tmp_covariate_df[, c("Lat", 
                                                     "Lon")], query = sample_i[Which, c("Lat", "Lon")], 
                         k = 1)
          X_ip[Which, ] = tmp_X[NN$nn.idx[, 1], , drop = FALSE]
        }
        NN = RANN::nn2(data = tmp_covariate_df[, c("Lat", "Lon")], 
                       query = latlon_g[, c("Lat", "Lon")], k = 1)
        X_gtp[, tI, ] = tmp_X[NN$nn.idx[, 1], , drop = FALSE]
      }
      X_itp = aperm(X_ip %o% rep(1, length(Year_Set)), perm = c(1, 
                                                                3, 2))
      if (any(is.na(X_itp))) 
        stop("Problem with `X_itp` in `make_covariates(.)")
      if (any(is.na(X_gtp))) 
        stop("Problem with `X_gtp` in `make_covariates(.)")
      if (any(apply(X_gtp, MARGIN = 2:3, FUN = sd) > 10 | apply(X_itp, 
                                                                MARGIN = 2:3, FUN = sd) > 10)) {
        warning("The package author recommends that you rescale covariates in `covariate_data` to have mean 0 and standard deviation 1.0")
      
      
      
      
      X1_itp = covariate_list$X_itp # 5358, 99, 3
      X1_gctp = aperm(outer(X1_gtp, rep(1, n_c)), c(1, 4, 2, 3)) # 2000, 1, 99, 3
      covariate_list = make_covariates(formula = X2_formula, 
                                       covariate_data = covariate_data, Year_i = t_i, 
                                       spatial_list = spatial_list)
      X2_gtp = covariate_list$X_gtp
      X2_itp = covariate_list$X_itp
      X2_gctp = aperm(outer(X2_gtp, rep(1, n_c)), c(1, 4, 2, 3))
    }
  }
  
  create_Xconfig = function(Xconfig_cp, n_c, n_p) {
    if (is.null(Xconfig_cp)) {
      Xconfig_cp = array(1, dim = c(n_c, n_p))
    } else {
      if (!is.array(Xconfig_cp) || !(all(dim(Xconfig_cp) == 
                                         c(n_c, n_p)))) {
        stop("`Xconfig_cp` has wrong dimensions")
      }
      if (!all(Xconfig_cp %in% c(-1, 0, 1, 2, 3))) {
        stop("`Xconfig_cp` has some wrong element(s)")
      }
      if (any(Xconfig_cp %in% -1)) {
        warning("Using `Xconfig_cp[] = -1` is unconventional and warrents caution")
      }
    }
    return(Xconfig_cp)
  }
  X1config_cp = create_Xconfig(Xconfig_cp = X1config_cp, n_c = n_c, n_p = dim(X1_gctp)[4])
  X2config_cp = create_Xconfig(Xconfig_cp = X2config_cp, n_c = n_c, n_p = dim(X2_gctp)[4])
  
  # Alright, so clearly something going on here with the number of parameters as this still things we just have 3...
  
  message("\n### Making TMB object")
  model_args_default = list(TmbData = data_list, RunDir = working_dir, 
                            Version = settings$Version, RhoConfig = settings$RhoConfig, 
                            loc_x = spatial_list$loc_x, Method = spatial_list$Method, 
                            build_model = build_model)
  model_args_input = combine_lists(input = extra_args, default = model_args_default, 
                                   args_to_use = formalArgs(make_model))
  tmb_list = do.call(what = make_model, args = model_args_input)
  if (run_model == FALSE | build_model == FALSE) {
    input_args = list(extra_args = extra_args, extrapolation_args_input = extrapolation_args_input, 
                      model_args_input = model_args_input, spatial_args_input = spatial_args_input, 
                      data_args_input = data_args_input)
    Return = list(data_frame = data_frame, extrapolation_list = extrapolation_list, 
                  spatial_list = spatial_list, data_list = data_list, 
                  tmb_list = tmb_list, year_labels = year_labels, 
                  years_to_plot = years_to_plot, settings = settings, 
                  input_args = input_args)
    class(Return) = "fit_model"
    return(Return)
  }
  if (silent == TRUE) 
    tmb_list$Obj$env$beSilent()
  if (test_fit == TRUE) {
    message("\n### Testing model at initial values")
    LogLike0 = tmb_list$Obj$fn(tmb_list$Obj$par)
    Gradient0 = tmb_list$Obj$gr(tmb_list$Obj$par)
    if (any(Gradient0 == 0)) {
      message("\n")
      stop("Please check model structure; some parameter has a gradient of zero at starting values\n", 
           call. = FALSE)
    }
    else {
      message("Looks good: All fixed effects have a nonzero gradient")
    }
  }
  message("\n### Estimating parameters")
  optimize_args_default1 = combine_lists(default = list(lower = tmb_list$Lower, 
                                                        upper = tmb_list$Upper, loopnum = 2), input = extra_args, 
                                         args_to_use = formalArgs(TMBhelper::fit_tmb))
  optimize_args_input1 = list(obj = tmb_list$Obj, savedir = NULL, 
                              newtonsteps = 0, bias.correct = FALSE, control = list(eval.max = 10000, 
                                                                                    iter.max = 10000, trace = 1), quiet = TRUE, getsd = FALSE)
  optimize_args_input1 = combine_lists(default = optimize_args_default1, 
                                       input = optimize_args_input1, args_to_use = formalArgs(TMBhelper::fit_tmb))
  parameter_estimates = do.call(what = TMBhelper::fit_tmb, 
                                args = optimize_args_input1)
  if (exists("check_fit") & test_fit == TRUE) {
    problem_found = VAST::check_fit(parameter_estimates)
    if (problem_found == TRUE) {
      message("\n")
      stop("Please change model structure to avoid problems with parameter estimates and then re-try; see details in `?check_fit`\n", 
           call. = FALSE)
    }
  }
  optimize_args_default2 = list(obj = tmb_list$Obj, lower = tmb_list$Lower, 
                                upper = tmb_list$Upper, savedir = working_dir, bias.correct = settings$bias.correct, 
                                newtonsteps = newtonsteps, bias.correct.control = list(sd = FALSE, 
                                                                                       split = NULL, nsplit = 1, vars_to_correct = settings$vars_to_correct), 
                                control = list(eval.max = 10000, iter.max = 10000, trace = 1), 
                                loopnum = 1)
  optimize_args_input2 = combine_lists(input = extra_args, 
                                       default = optimize_args_default2, args_to_use = formalArgs(TMBhelper::fit_tmb))
  optimize_args_input2 = combine_lists(input = list(startpar = parameter_estimates$par), 
                                       default = optimize_args_input2)
  parameter_estimates = do.call(what = TMBhelper::fit_tmb, 
                                args = optimize_args_input2)
  if ("par" %in% names(parameter_estimates)) {
    Report = tmb_list$Obj$report()
    ParHat = tmb_list$Obj$env$parList(parameter_estimates$par)
  }
  else {
    Report = ParHat = "Model is not converged"
  }
  input_args = list(extra_args = extra_args, extrapolation_args_input = extrapolation_args_input, 
                    model_args_input = model_args_input, spatial_args_input = spatial_args_input, 
                    optimize_args_input1 = optimize_args_input1, optimize_args_input2 = optimize_args_input2, 
                    data_args_input = data_args_input)
  Return = list(data_frame = data_frame, extrapolation_list = extrapolation_list, 
                spatial_list = spatial_list, data_list = data_list, 
                tmb_list = tmb_list, parameter_estimates = parameter_estimates, 
                Report = Report, ParHat = ParHat, year_labels = year_labels, 
                years_to_plot = years_to_plot, settings = settings, 
                input_args = input_args, X1config_cp = X1config_cp, 
                X2config_cp = X2config_cp, covariate_data = covariate_data, 
                X1_formula = X1_formula, X2_formula = X2_formula, Q1config_k = Q1config_k, 
                Q2config_k = Q1config_k, catchability_data = catchability_data, 
                Q1_formula = Q1_formula, Q2_formula = Q2_formula)
  Return$effects = list()
  if (!is.null(catchability_data)) {
    catchability_data_full = data.frame(catchability_data, 
                                        linear_predictor = 0)
    Q1_formula_full = update.formula(Q1_formula, linear_predictor ~ 
                                       . + 0)
    call_Q1 = lm(Q1_formula_full, data = catchability_data_full)$call
    Q2_formula_full = update.formula(Q2_formula, linear_predictor ~ 
                                       . + 0)
    call_Q2 = lm(Q2_formula_full, data = catchability_data_full)$call
    Return$effects = c(Return$effects, list(call_Q1 = call_Q1, 
                                            call_Q2 = call_Q2, catchability_data_full = catchability_data_full))
  }
  if (!is.null(covariate_data)) {
    covariate_data_full = data.frame(covariate_data, linear_predictor = 0)
    X1_formula_full = update.formula(X1_formula, linear_predictor ~ 
                                       . + 0)
    call_X1 = lm(X1_formula_full, data = covariate_data_full)$call
    X2_formula_full = update.formula(X2_formula, linear_predictor ~ 
                                       . + 0)
    call_X2 = lm(X2_formula_full, data = covariate_data_full)$call
    Return$effects = c(Return$effects, list(call_X1 = call_X1, 
                                            call_X2 = call_X2, covariate_data_full = covariate_data_full))
  }
  class(Return) = "fit_model"
  return(Return)
}
