###### DEVELOPMENT VERSION OF VAST SHOULD BE USED
library(VAST)
library(splines)  # Used to include basis-splines
library(effects)  # Used to visualize covariate effects
library(mgcv)

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

# Define formula, here using a b-spline basis with the splines library. For the cubic regression spline, we want to set "degree = 3".
X1_formula = ~ bs( BOT_DEPTH, degree=3, intercept=FALSE)
X2_formula = ~ bs( BOT_DEPTH, degree=3, intercept=FALSE)

# If all covariates as "static" (not changing among years), then set Year = NA to cause values to be duplicated internally for all values of Year
example$covariate_data[,'Year'] = NA

# In the Wiki example, covariates are rescaled to have an SD >0.1 and <10 (for numerical stability) with the following line. The mgcv::gam function I believe automatically scales/centers any continuous covriates. So, rather than doing the rescaling as in the wiki example, going to scale/center instead for direct comparison between the two approaches.
# example$covariate_data[,'BOT_DEPTH'] = example$covariate_data[,'BOT_DEPTH'] / 100
example$covariate_data[,'BOT_DEPTH']<- scale(example$covariate_data[,'BOT_DEPTH'])

# Another thing we need to adjust is the default settings and the ObsModel default for "purpose = index2", which uses "ObsModel = c(2, 1)" and specifies the "Poisson" link delta model. I'm not entirely sure how to implement that in mgcv::gam. So, instead going to use "ObsModel = c(2, 0)", which is the "traditional" two-stage delta log normal GAM.
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

# Must add data-frames to global environment
covariate_data_full = fit$effects$covariate_data_full
catchability_data_full = fit$effects$catchability_data_full

# Plot 1st linear predictor
pred = Effect.fit_model( fit,
                         focal.predictors = c("BOT_DEPTH"),
                         which_formula = "X1", 
                         xlevels = 100)
plot(pred)


# Now, trying to do the same thing with mgcv. To do this, relies heavily on Gavin Simpson's excellent tutorial here: https://fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/. First, convert "response" to 0 and 1 and then also get the bottom depth covariate
gam_df<- data.frame("Response" = ifelse(example$sampling_data[,'Catch_KG'] > 0, 1, 0), "BOT_DEPTH" = example$covariate_data[,'BOT_DEPTH'])

# Specify and fit the model. According to Gavin's blog post and the R help file, a gam with s(variable, bs = "bs", m = c(3, 2)) is the "garden variety cubic B spline with second order penalty" -- so penalty placed on the second derivative of the spline, or its curvature (as opposed to m[2] = 1, which would penalize the first derivative and departures from a linear function). With splines::bs(), I wasn't sure if there was actually any penalty being applied and tried "m = c(3,0)" first. This gave something similar, but not identical. So, went ahead and tried "m = c(3, 2)" instead. 
gam_fit_2pen<- gam(Response ~ s(BOT_DEPTH, bs = "bs", m = c(3, 2)), data = gam_df, family = "binomial")
plot(gam_fit_2pen, main = "MGCV::GAM, m = c(3, 2)")
plot(pred, main = "VAST SPLINES::BS()")
