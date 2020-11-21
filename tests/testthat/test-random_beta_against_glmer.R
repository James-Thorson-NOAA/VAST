
# Tutorial: http://r-pkgs.had.co.nz/tests.html
# And example see: https://github.com/ss3sim/ss3sim/tree/master/tests/testthat
context("Testing examples")

# Eastern Bering Sea pollcok
test_that("Covariates give identical results to glm(.) ", {
  skip_on_travis()
  library(lme4)

  # load data set
  data( EBS_pollock_data, package="FishStatsUtils" )
  Data = data.frame(EBS_pollock_data,
    "year_factor" = factor(EBS_pollock_data$year,levels=sort(unique(EBS_pollock_data$year))),
    "pres" = ifelse(EBS_pollock_data$catch>0,1,0),
    "AreaSwept_km2" = 0.01 )

  # Make settings (turning off bias.correct to save time for example)
  settings = make_settings( n_x=100,
    Region="Eastern_Bering_Sea",
    purpose="index2",
    use_anisotropy=FALSE,
    strata.limits=example$strata.limits,
    bias.correct=FALSE,
    fine_scale=TRUE,
    FieldConfig=c(0,0,0,0),
    ObsModel=c(1,0),
    RhoConfig=c(1,1,0,0),
    max_cells=Inf )

  # Run model -- Lognormal
  fit = fit_model( settings=settings,
    Lat_i=EBS_pollock_data[,'lat'],
    Lon_i=EBS_pollock_data[,'long'],
    t_i=EBS_pollock_data[,'year'],
    b_i=EBS_pollock_data[,'catch'],
    a_i=rep(0.01,nrow(EBS_pollock_data)),
    working_dir=multispecies_example_path )

  # Function to facilitate comparison
  extract = function(vec, name, remove=FALSE){
    if(remove==FALSE) return( as.vector(vec[grep(name,names(vec))]) )
    if(remove==TRUE) return( as.vector(vec[-grep(name,names(vec))]) )
  }

  # Glm0 -- Bernoulli
  Glm0 = glmer( formula= pres ~ 1 + (1 | year_factor),
    data=Data, family="binomial" )

  # Glm1 -- Lognormal
  Glm1 = lmer( formula= log(catch) ~ 1 + (1 | year_factor),
    data=Data[which(Data$pres==1),], offset=log(AreaSwept_km2) )

  # Comparison with Glm0
  expect_equal( fit$Report$beta1_tc[,1], coef(Glm0)$year_factor[,'(Intercept)'], tolerance=0.001 )

  # Comparison with Glm1
  expect_equal( fit$Report$beta2_tc[,1] - exp(2*fit$parameter_estimates$par['logSigmaM'])/2, coef(Glm1)$year_factor[,'(Intercept)'], tolerance=0.001 )

})

