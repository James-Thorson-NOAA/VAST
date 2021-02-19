###### DEVELOPMENT VERSION OF VAST SHOULD BE USED
library(VAST)
library(splines)  # Used to include basis-splines
library(effects)  # Used to visualize covariate effects
library(mgcv)
library(patchwork)

#####
## Exact example from Wiki
#####
# load data set
example = load_example( data_set="EBS_pollock" )
covariate_data = data.frame( "Lat"=0, "Lon"=0, "Year"=example$covariate_data[,'Year'],
                             "CPE"=(example$covariate_data[,'AREA_SUM_KM2_LTE2']-100000)/100000)

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x=100,
                          Region=example$Region,
                          purpose="index2",
                          bias.correct=FALSE)

# Change Beta1 to AR1, to allow linear covariate effect
settings$RhoConfig['Beta1'] = 4

# Define formula.
X1_formula = ~ CPE
X2_formula = ~ 1

#
X1config_cp = array(3,dim=c(1,1))

# Run model
fit = fit_model( "settings" = settings,
                 Lat_i = example$sampling_data[,'Lat'],
                 Lon_i = example$sampling_data[,'Lon'],
                 t_i = example$sampling_data[,'Year'],
                 b_i = example$sampling_data[,'Catch_KG'],
                 a_i = example$sampling_data[,'AreaSwept_km2'],
                 X1_formula = X1_formula,
                 X2_formula = X2_formula,
                 X1config_cp = X1config_cp,
                 covariate_data = covariate_data,
                 test_fit = FALSE )

#####################
# Effects package
#####################

library(effects)  # Used to visualize covariate effects

# Must add data-frames to global environment (hope to fix in future). 
covariate_data_full = fit$effects$covariate_data_full
catchability_data_full = fit$effects$catchability_data_full

# Plot 1st linear predictor, but could use `transformation` to apply link function
pred = Effect.fit_model( fit, 
                         focal.predictors = c("CPE"),
                         which_formula = "X1",
                         xlevels = 100,
                         transformation = list(link=identity, inverse=identity) )
plot(pred)

#####################
# Testing "old" models
#####################
# First, remove the "effects" slot from "fit"
fit_old<- fit
fit_old$effects<- NULL

pred_old = Effect.fit_model( fit_old, 
                         focal.predictors = c("CPE"),
                         which_formula = "X1",
                         xlevels = 100,
                         transformation = list(link=identity, inverse=identity) )
plot(pred_old)

#####
## For comparison with GAMs, leveraging "covariate example" wiki
####
# load data set
# see `?load_example` for list of stocks with example data
# that are installed automatically with `FishStatsUtils`.
example = load_example( data_set="covariate_example" )

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x=100,
                          Region=example$Region,
                          purpose="index2",
                          use_anisotropy=FALSE,
                          bias.correct=FALSE,
                          fine_scale=TRUE )

# Define formula, here using a b-spline basis with the splines library. I have some questions here about the use of "knots = 3." From the help file, the knots define the internal breakpoints of the spline. If NULL (the default), then the b-spline basis will be an ordinary polynomial regression. Typically values are usually the mean or median for one knot or quantiles for more knots. So, I am not sure if that is really what we want? Instead, seems like we would want to set "degree = 3", with a default of 3 for a cubic spline. Going with that for now as I think this will also make it easier for me to replicate things with the mgcv library.  
X1_formula = ~ bs( BOT_DEPTH, degree=3, intercept=FALSE)
X2_formula = ~ bs( BOT_DEPTH, degree=3, intercept=FALSE)

# If all covariates as "static" (not changing among years),
#  then set Year = NA to cause values to be duplicated internally for all values of Year
# If using a mix of static and dynamic covariates,
#  then duplicate rows for static covariates for every value of Year
# Here, all covariates are static, so I'm using the first approach.
example$covariate_data[,'Year'] = NA

# Rescale covariates being used to have an SD >0.1 and <10 (for numerical stability)
# example$covariate_data[,'BOT_DEPTH'] = example$covariate_data[,'BOT_DEPTH'] / 100

# The GAM mgcv package rescales the predictor variables automatically (center and scale). So doing that instead.
example$covariate_data[,'BOT_DEPTH']<- scale(example$covariate_data[,'BOT_DEPTH'])

# Another thing we need to adjust is the default settings and the ObsModel default for "purpose = index2", which uses "ObsModel = c(2, 1)" and specifies the "Poisson" link delta model. I'm not entirely sure how to implement this one in mgcv::gam. So, instead going to use "ObsModel = c(2, 0)", which is the "traditional" two-stage delta log normal GAM.

settings$ObsModel<- c(2, 0)

# Run model
fit = fit_model( "settings" = settings,
                 Lat_i = example$sampling_data[,'Lat'],
                 Lon_i = example$sampling_data[,'Lon'],
                 t_i = example$sampling_data[,'Year'],
                 b_i = example$sampling_data[,'Catch_KG'],
                 a_i = example$sampling_data[,'AreaSwept_km2'],
                 X1_formula = X1_formula,
                 X2_formula = X2_formula,
                 covariate_data = example$covariate_data )

library(effects)  # Used to visualize covariate effects

# Must add data-frames to global environment (hope to fix in future)
covariate_data_full = fit$effects$covariate_data_full
catchability_data_full = fit$effects$catchability_data_full

# Plot 1st linear predictor
pred = Effect.fit_model( fit,
                         focal.predictors = c("BOT_DEPTH"),
                         which_formula = "X1", 
                         xlevels = 100)
plot(pred)


# Now, trying to do the same thing with mgcv. To do this, relying heavily on Gavin Simpson's excellent tutorial here: https://fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/. First, convert "response" to 0 and 1 and then also get the bottom depth variable.
gam_df<- data.frame("Response" = ifelse(example$sampling_data[,'Catch_KG'] > 0, 1, 0), "BOT_DEPTH" = example$covariate_data[,'BOT_DEPTH'])

# Specify and fit the model. According to Gavin's blog post and the R help file, a gam with s(variable, bs = "bs", m = c(3, 2)) is a "cubic B spline with second order penalty" -- so penalty placed on the second derivative of the spline, or its curvature. With the splines::bs(), I am not sure if there is actually any penalty being applied? Going to try having m[2] = 0, then, which should remove penalty from being included with estimating the smooth.
gam_fit<- gam(Response ~ s(BOT_DEPTH, bs = "bs", m = c(3, 0)), data = gam_df, family = "binomial")
plot(gam_fit)
plot(pred)

# Not quite...shape looks a bit similar. Just for kicks, going to try the "second order penalty" and just see if that changes things
gam_fit_2pen<- gam(Response ~ s(BOT_DEPTH, bs = "bs", m = c(3, 2)), data = gam_df, family = "binomial")
plot(gam_fit_2pen)
plot(pred)

# Looks really good!