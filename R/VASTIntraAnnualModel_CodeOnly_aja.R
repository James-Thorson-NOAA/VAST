#####
## Fitting seasonal and intra-annual models with VAST
#####

# Original code -----------------------------------------------------------
# Install latest development version of VAST. If you have issues here, please check the VAST GitHub page and then be sure to also look through the issues (open and closed) as many have provided helpful fixes to common installation problems.
install.packages("devtools")
library(devtools)
install_github("james-thorson/VAST@development")
library(VAST)
library(tidyverse)
library(here)

# Older version of TMB to avoid rtweedie issues (https://github.com/James-Thorson-NOAA/VAST/issues/276#issue-804990507)
install_version("TMB", version = "1.7.18", repos = "http://cran.us.r-project.org")
library(TMB)

#####
## Data processing
#####
# Load data and quick exploration of structure 
load(here("R", "VAST_SeasonalExampleData.RData"))

# Quick look at data structure 
str(svdata)
summary(svdata)

# Copying svdata to a new object 
sampling_data<- svdata 

# Set of seasons and years. The DFO spring survey usually occurs before the NOAA NEFSC spring survey, so ordering accordingly.
season_set<- c("DFO", "SPRING", "FALL")
year_set<- sort(unique(sampling_data[,'year']))

# Create a grid with all unique combinations of seasons and years and then combine these into one "season_year" variable
seasonyear_grid<- expand.grid("season" = season_set, "year" = year_set)
t_levels<- apply(seasonyear_grid, MARGIN = 1, FUN = paste, collapse = "_")
t_labels<- round(seasonyear_grid[,'year'] + (as.numeric(factor(seasonyear_grid[,'season'], levels = season_set))-1)/length(season_set), digits=1)

# Similar process, but for the observations
seasonyear_i<- apply(sampling_data[,c("season","year")], MARGIN = 1, FUN = paste, collapse = "_")
seasonyear_i<- factor(seasonyear_i, levels = t_levels)

# Add the season_year factor column to our sampling_data data set
sampling_data <- cbind(sampling_data, "season_year" = seasonyear_i)
sampling_data$season<- factor(sampling_data$season, levels = season_set)

# Some last processing steps
sampling_data<- sampling_data[, c("year", "season", "season_year", "latitude", "longitude", "swept", "weight")]
colnames(sampling_data)<- c("year", "season", "season_year", "latitude", "longitude", "swept", "response")

# Lastly, adding in dummy data for missing season_years
# Make dummy observation for each season-year combination
dummy_data<- data.frame(year = seasonyear_grid[,'year'], season = seasonyear_grid[,'season'], season_year = t_levels, latitude = mean(sampling_data[,'latitude']), longitude = mean(sampling_data[,'longitude']), swept = mean(sampling_data[,'swept']), response = 0, dummy = TRUE)

# Combine with sampling data
sampling_data<- rbind(cbind(sampling_data, dummy = FALSE), dummy_data)

#####
## Initial VAST build
#####
# Extrapolation grid
maximum_distance_from_sample<- 10
grid_dim_km<- c(5, 5)
region_code<- "northwest_atlantic"
strata.limits<- list("Yellowtail" = c(1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200, 1210))

# Observation model -- "Poisson-link delta model" with lognormal positive catch rate distribution
ObsModel<- c(1, 1)

# Number of spatial and spatio-temporal factors to use for each linear predictor. Here, using a single species model, and turning spatial and spatio-temporal variability effects for both linear predictors "on".
FieldConfig<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)

# Setting the structure of parameters across season_year surveys. Here, setting intercepts for the first and second linear predictors (Beta) to be estimated as random effects, which are constant across season_years, and then estimating a first order auto-regressive structure for the spatio-temporal variability (Epsilon) of both linear predictors. This implementation facilitates estimating occurrence at unsampled times/areas.
RhoConfig<- c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 4, "Epsilon2" = 4)
Options<- c('treat_nonencounter_as_zero' = TRUE)

# OutFile
OutFile<- "~/Desktop/VASTIntraAnnualModel/"
if(!file.exists(OutFile)){
  dir.create(OutFile)
}

# Set the run directory
RunDir<- OutFile

# Make settings
settings<- make_settings(n_x = 100, Region = region_code, purpose = "index2", FieldConfig = FieldConfig, RhoConfig = RhoConfig, ObsModel = ObsModel, bias.correct = FALSE, strata.limits = strata.limits, Options = Options)

Use_REML<- FALSE
fit<- fit_model("settings" = settings, "Lat_i" = sampling_data[,'latitude'], "Lon_i" = sampling_data[,'longitude'], "t_i" = as.numeric(sampling_data[,"season_year"])-1, "c_i" = rep(0, nrow(sampling_data)), "b_i" = sampling_data[,'response'], "a_i" = sampling_data[,'swept'], "run_model" = FALSE, "observations_LL" = cbind("Lat" = sampling_data[,'latitude'], "Lon" = sampling_data[,'longitude']), "maximum_distance_from_sample" = maximum_distance_from_sample, "grid_dim_km" = grid_dim_km, "working_dir" = RunDir, "PredTF_i" = sampling_data[,'dummy'], "Use_REML" = Use_REML)

#####
## Modifying covariate matrices for season and year as well as the parameter effects of season and year and rebuilding
#####
## Season covariate matrices
# For observations
# Get the season of each observation (1 = DFO, 2 = SPRING, 3 = FALL)
season_i<- match(sampling_data[,'season'], season_set)

# Now, convert this vector into a design style model matrix -- one row per observation, three columns (one for each season), and a "1" placed in correct season column according to the season of the observation
Xseason_ip<- ThorsonUtilities::vector_to_design_matrix(season_i)
# Expand this to have one of these design matrices per time step in the model. 
Xseason_itp<- aperm(Xseason_ip %o% rep(1, fit$data_list$n_t), c(1,3,2))

# Similar process, but now for the mesh locations
season_t<- match(seasonyear_grid[,'season'], season_set)
Xseason_tp<- ThorsonUtilities::vector_to_design_matrix(season_t)
Xseason_gtp<- rep(1, fit$data_list$n_g) %o% Xseason_tp

## Year covariate matrices
# For observations
year_i<- match(sampling_data[,'year'], year_set)
# Vector to design style model matrix, one row per observation, 33 columns (one for each year), and a "1" placed in correct season column according to the season of the observation
Xyear_ip<- ThorsonUtilities::vector_to_design_matrix(year_i)
# Expand
Xyear_itp<- aperm(Xyear_ip %o% rep(1, fit$data_list$n_t), c(1,3,2))

# For mesh locations
year_t<- match(seasonyear_grid[,'year'], year_set)
Xyear_tp<- ThorsonUtilities::vector_to_design_matrix(year_t)
Xyear_gtp<- rep(1, fit$data_list$n_g) %o% Xyear_tp

## Season and year parameter effects
# First, an empty array for the season effects (rows = linear predictors, columns = categories -- here just 1 with single species model), depth = seasons.
season_Xconfig_zcp<- array(NA, dim = c(2, 1, ncol(Xseason_ip)), dimnames = list(NULL, NULL, colnames(Xseason_ip)))

# Now, specifying the effects of each. See ?fit_model and the X1config_cp description for details.
XiConfig<- c("Xi1_season" = 3, "Xi1_year" = 1, "Xi2_season" = 1, "Xi2_year" = 1) 

# Next, supplying the Xi_season info from earlier, where first linear predictor will have a zero-centered, spatially varying effect and second linear predictor will have just a linear effect for each of the seasons
season_Xconfig_zcp[1,,]<- XiConfig["Xi1_season"]
season_Xconfig_zcp[2,,]<- XiConfig["Xi2_season"]

# Similar process, but now with the year effect
year_Xconfig_zcp<- array(NA, dim = c(2, 1, ncol(Xyear_ip)), dimnames = list(NULL, NULL, colnames(Xyear_ip)) )
year_Xconfig_zcp[1,,]<- XiConfig["Xi1_year"]
year_Xconfig_zcp[2,,]<- XiConfig["Xi2_year"]

# Adjustments. From the original code:
# If Design="Both" and no RhoConfig structure, then drop two DF (corner constraint for season and year main effects, one for corner constraint of main effects + interactions)
# If using RhoConfig structure, drop one less DF. 
# If Design!="Both, drop one less DF
# SOLUTION:  Drop one level from each Season and Year effect
Design<- "Both"
FUN = function(num){
  if(num == 1) return(0)
  if(num == 3) return(2)
  return(num)
}

# Now, applying FUN to first season and year in both linear predictors
season_Xconfig_zcp[1,,1] # 2, indicating zero-centered spatially varying effect of first season on the first linear predictor. 
season_Xconfig_zcp[2,,1] # 0, indicating that there will be no effect of first season on the second linear predictor

# Similar process, but now with year
year_Xconfig_zcp[1,,1]<- FUN(year_Xconfig_zcp[1,,1])
year_Xconfig_zcp[2,,1]<- FUN(year_Xconfig_zcp[2,,1])

# Finally, combine all the information into the three core objects (X_itp, X_gtp and Xconfig_zcp)
# Combine season and year matrices for observations and then mesh locations
Xconfig_zcp = X_itp = X_gtp = NULL
X_itp<- abind::abind(X_itp, Xseason_itp, along = 3)
str(X_itp)
X_itp<- abind::abind(X_itp, Xyear_itp, along = 3)
str(X_itp)

X_gtp<- abind::abind( X_gtp, Xseason_gtp, along=3 )
X_gtp<- abind::abind( X_gtp, Xyear_gtp, along=3 )
str(X_gtp)

Xconfig_zcp<- abind::abind(Xconfig_zcp, season_Xconfig_zcp)
Xconfig_zcp<- abind::abind(Xconfig_zcp, year_Xconfig_zcp)
Xconfig_zcp

# Rebuilding with new X_itp, X_gtp and Xconfig_zcp
fit_seas<- fit_model("settings" = settings, "Lat_i" = sampling_data[,'latitude'], "Lon_i" = sampling_data[,'longitude'], "t_i" = as.numeric(sampling_data[,"season_year"])-1, "c_i" = rep(0, nrow(sampling_data)), "b_i" = sampling_data[,'response'], "a_i" = sampling_data[,'swept'], newtonsteps = 1, getsd = TRUE, "observations_LL" = cbind("Lat" = sampling_data[,'latitude'], "Lon" = sampling_data[,'longitude']),  "maximum_distance_from_sample" = maximum_distance_from_sample, "grid_dim_km" = grid_dim_km, "run_model" = FALSE, "X_itp" = X_itp, "X_gtp" = X_gtp, "Xconfig_zcp" = Xconfig_zcp, "test_fit" = FALSE, working_dir = RunDir, "PredTF_i" = sampling_data[,'dummy'], "Use_REML" = Use_REML)

#####
## Adjusting Map settings for log_sigmaXi and fitting final model
#####
# Apply customization
Map$log_sigmaXi1_cp<- factor(c(rep(Map$log_sigmaXi1_cp[1], dim(Xseason_itp)[3]), rep(Map$log_sigmaXi1_cp[dim(Xseason_itp)[3]+1], dim(Xyear_itp)[3])))
Map$log_sigmaXi2_cp<- factor(c(rep(Map$log_sigmaXi2_cp[1], dim(Xseason_itp)[3]), rep(Map$log_sigmaXi2_cp[dim(Xseason_itp)[3]+1], dim(Xyear_itp)[3])))
str(Map$log_sigmaXi1_cp) # Vector of length 36 (3 seasons, 33 years), [1] = 1, [2] = 1, [3] = 1 and rest  = NA.

# Refit with new mapping argument -- run_model set to FALSE for now as it took a while on my computer to run (could have just been because I had multiple R sessions going)
run_model_use<- FALSE
fit_seas_orig<- fit_model("settings" = settings, "Lat_i" = sampling_data[,'latitude'], "Lon_i" = sampling_data[,'longitude'], "t_i" = as.numeric(sampling_data[,"season_year"])-1, "c_i" = rep(0,nrow(sampling_data)), "b_i" = sampling_data[,'response'], "a_i" = sampling_data[,'swept'], "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "observations_LL" = cbind("Lat" = sampling_data[,'latitude'],"Lon" = sampling_data[,'longitude']), "maximum_distance_from_sample" = maximum_distance_from_sample, "grid_dim_km" = grid_dim_km, "run_model" = run_model_use, "X_itp" = X_itp, "X_gtp" = X_gtp, "Xconfig_zcp" = Xconfig_zcp, "test_fit" = FALSE, "Map" = Map, working_dir = RunDir, "PredTF_i" = sampling_data[,'dummy'],  "Use_REML" = Use_REML, "getJointPrecision" = FALSE)


# New approach, leveraging model formula interface ------------------------
#####
## Beginnings
#####

# Before proceeding, make sure to source the following functions: "fit_model_aja.R", "make_data_aja.R", and "make_covariates_aja.R"
funs_aja<- paste(here("R", c("fit_model_aja.R", "make_data_aja.R", "make_covariates_aja.R")))
sapply(funs_aja, source)

#####
## Data processing
#####
str(sampling_data)

# Creating a sample and covariate data frame. The `make_covariates` function, which uses data supplied to the covariate_data argument in `fit_model` is looking for specific column names.  
# Create sample data
samp_dat<- data.frame("spp" = rep("Yellowtail", nrow(sampling_data)), "Year" = as.numeric(sampling_data$season_year)-1, "Season" = sampling_data$season, "Season_Year" = sampling_data$season_year, "Lat" = sampling_data$latitude, "Lon" = sampling_data$longitude, "Response" = sampling_data$response, "Swept" = sampling_data$swept, "Dummy" = sampling_data$dummy)

# Covariate data. Note here, case sensitive!
cov_dat<- data.frame("spp" = rep("Yellowtail", nrow(sampling_data)), "Year" = as.numeric(sampling_data$season_year)-1, "Year_Cov" = factor(sampling_data$year, levels = year_set), "Season" = sampling_data$season, "Lat" = sampling_data$latitude, "Lon" = sampling_data$longitude)

#####
## Model settings: formula and Xconfig object
#####
# Creating model formula
formula_use<- ~ Season + Year_Cov

# First, let's start off with a simple situation where we want to estimate each season and year main effect as spatially-varying zero-centered in both linear predictors.
season_effect<- rep(2, length(levels(cov_dat$Season)))
year_cov_effect<- rep(2, length(unique(cov_dat$Year_Cov)))

# Alright, since we only have one species/category, this is going to make very simple `Xconfig_cp` matrices for each of the linear predictors with all 2's.
X1config_cp_use<- matrix(data = c(season_effect, year_cov_effect), nrow = 1)
X2config_cp_use<- matrix(data = c(season_effect, year_cov_effect), nrow = 1)

# Force match of Xconfig objects
X1config_cp_use[,]<- fit_seas_orig$data_list$X1config_cp
X2config_cp_use[,]<- fit_seas_orig$data_list$X2config_cp

#####
## Model fit with modified `fit_model_aja` wrapper, which uses `make_data_aja` and `make_covariates_aja` -- all just edited to allow the contrasts list to pass through to the call for model.matrix 
#####

# Fit model -- again run_model set to FALSE here because it took a while
run_model_use<- FALSE
fit_seas_form<- fit_model_aja("settings" = settings, "Lat_i" = samp_dat[, 'Lat'], "Lon_i" = samp_dat[, 'Lon'], "t_i" = samp_dat[, 'Year'], "c_i" = rep(0, nrow(samp_dat)), "b_i" = samp_dat[, 'Response'], "a_i" = samp_dat[, 'Swept'], "X1config_cp" = X1config_cp_use, "X2config_cp" = X2config_cp_use, "covariate_data" = cov_dat, "X1_formula" = formula_use, "X2_formula" = formula_use, contrasts = list(Season=contrasts(covariate_df$Season, contrasts = FALSE), Year_Cov = contrasts(covariate_df$Year_Cov, contrasts = FALSE)), "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "observations_LL" = cbind("Lat" = samp_dat[, 'Lat'], "Lon" = samp_dat[, 'Lon']), "maximum_distance_from_sample" = maximum_distance_from_sample, "grid_dim_km" = grid_dim_km, "run_model" = run_model_use, "test_fit" = FALSE, working_dir = RunDir, "PredTF_i" = samp_dat[, 'Dummy'], "Use_REML" = Use_REML, "getJointPrecision" = FALSE)

# Comparison ------------------------
fit_mods<- list(fit_seas_orig, fit_seas_form)
names(fit_mods)<- c("Original", "New")
params_list<- vector("list", length(fit_mods))
settings_list<- vector("list", length(fit_mods))
Xconfig_comp<- vector("list", length(fit_mods))

for(i in seq_along(fit_mods)){
  mod<- fit_mods[[i]]
  field_config<- mod$data_list$FieldConfig
  rho_config<- mod$data_list$RhoConfig
  settings_list[[i]]<- list(field_config, rho_config)
  Xconfig_comp[[i]]<- mod$data_list$Xconfig_zcp
  params<- mod$parameter_estimates$diagnostics
  params$model<- rep(names(fit_mods)[i], nrow(params))
  params_list[[i]]<- params
}

# Visualizing the differences...want the difference and to check the sign?
params_df<- do.call(rbind.data.frame, params_list)

# More descriptive parameter names?
unique(params_df$Param)
param_namesA<- c("ln_H_input1", "ln_H_input2", "beta1_ft", c(paste("gamma1_cp_", seq(1:34), sep = "")), "L_omega1_z", "L_epsilon1_z", "logkappa1", "Epsilon_rho1_f", "log_sigmaXi1_cp1", "beta2_ft", c(paste("gamma2_cp_", seq(1:34), sep = "")), "L_omega2_z", "L_epsilon2_z", "logkappa2", "Epsilon_rho2_f", "logSigmaM")
param_namesB<- c("ln_H_input1", "ln_H_input2", "beta1_ft", c(paste("gamma1_cp_", seq(1:34), sep = "")), "L_omega1_z", "L_epsilon1_z", "logkappa1", "Epsilon_rho1_f", c(paste("log_sigmaXi1_cp", seq(1:3), sep = "")), "beta2_ft", c(paste("gamma2_cp_", seq(1:34), sep = "")), "L_omega2_z", "L_epsilon2_z", "logkappa2", "Epsilon_rho2_f", "logSigmaM")
params_df$Param_Plot<- c(param_namesA, param_namesB)
params_df$Param_Plot<- factor(params_df$Param_Plot, levels = c("ln_H_input1", "ln_H_input2", "beta1_ft", c(paste("gamma1_cp_", seq(1:34), sep = "")), "L_omega1_z", "L_epsilon1_z", "logkappa1", "Epsilon_rho1_f", c(paste("log_sigmaXi1_cp", seq(1:3), sep = "")), "beta2_ft", c(paste("gamma2_cp_", seq(1:34), sep = "")), "L_omega2_z", "L_epsilon2_z", "logkappa2", "Epsilon_rho2_f", "logSigmaM"))

comp_plot<- ggplot() +
  geom_point(data = params_df, aes(x = Param_Plot, y = MLE, color = model, shape = model), size = 4) +
  scale_shape_manual(name = "Model Code", values = c(21, 3)) +
  scale_color_manual(name = "Model Code", values = c("#1b9e77", "#7570b3")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
comp_plot
