Description
=============

VAST
* Is an R package for implementing a spatial delta-generalized linear mixed model (delta-GLMM) for multiple categories (species, size, or age classes) when standardizing survey or fishery-dependent data.
* Builds upon a previous R package `SpatialDeltaGLMM` (public available [here](https://github.com/nwfsc-assess/geostatistical_delta-GLMM)), and has unit-testing to automatically confirm that `VAST` and `SpatialDeltaGLMM` give identical results (to the 3rd decimal place for parameter estimates) for several varied real-world case-study examples
* Has built in diagnostic functions and model-comparison tools
* Is intended to improve analysis speed, replicability, peer-review, and interpretation of index standardization methods

Background
* This tool is designed to estimate spatial variation in density using spatially referenced data, with the goal of habitat associations (correlations among species and with habitat) and estimating total abundance for a target species in one or more years.
* The model builds upon spatio-temporal delta-generalized linear mixed modelling techniques (Thorson Shelton Ward Skaug 2015 ICESJMS), which separately models the proportion of tows that catch at least one individual ("encounter probability") and catch rates for tows with at least one individual ("positive catch rates").
* Submodels for encounter probability and positive catch rates by default incorporate variation in density among years (as a fixed effect), and can incorporate variation among sampling vessels (as a random effect, Thorson and Ward 2014) which may be correlated among categories (Thorson Fonner Haltuch Ono Winker In press).
* Spatial and spatiotemporal variation are approximated as Gaussian Markov random fields (Thorson Skaug Kristensen Shelton Ward Harms Banante 2014 Ecology), which imply that correlations in spatial variation decay as a function of distance.

User resources for learning about VAST
=============
There are eleven main resources for learning about VAST:

*  *Model structure*:  Please see the [User Manual](https://github.com/James-Thorson/VAST/blob/master/manual/VAST_model_structure.pdf) for a document listing model equations and relating them to the input/output used in R.
*  *Change documentation*:  Please see [NEWS](https://github.com/James-Thorson/VAST/blob/master/manual/NEWS.pdf) for a document listing major changes in each numbered release; this is useful to track development changes (interface improvements, bug fixes, issues regarding dependencies, etc.) over time. 
*  *Guidance for user decisions*:  Please see [Thorson-2019](https://www.sciencedirect.com/science/article/abs/pii/S0165783618302820) for guidance regarding the 15 major decisions needed in every VAST model
*  *High-level wrapper functions*:  I have recently added high-level wrapper functions, which provide a [gentle introduction](https://github.com/James-Thorson/VAST/wiki/Simple-example) to running `VAST`
*  *Examples*:  The wiki also includes example code to run VAST for many common [use-cases](https://github.com/James-Thorson-NOAA/VAST/wiki), as listed in the right-hand-side toolbar at that link.
*  *R-help documentation*:  Please see the standard R-help documentation, e.g., by typing `?fit_model` or `?make_data` in the R-terminal after installing the package and loading it using `library(VAST)`.
*  *Publications*:  Please see the [publications list](https://github.com/nwfsc-assess/geostatistical_delta-GLMM/wiki/Applications) to identify peer-reviewed publications regarding individual features.  These publications include statistical theory and model testing.
*  *List-serv*: Consider joining the [FishStats listserve](https://groups.google.com/forum/#!forum/fishstats-listserv) for 4-6 updates per year, including training classes.
*  *Slack channel*:  A [slack channel](https://join.slack.com/t/vastsupportgroup/shared_invite/zt-x53wupkg-U5cGF7gsHv9uDKyNMfPpaw) was developed by J. Morano and colleagues to allow real-time, casual discussions among new and longtime users.
*  *Talks available online*:  We post recorded talks and seminars [online](https://www.youtube.com/channel/UCNgFcss1X9Hgox3eWSTRI3Q/)
*  *Issue-tracker*:  Before posting new issues, users should explore the previous issues in the github issue tracker for [VAST](https://github.com/James-Thorson/VAST/issues), [SpatialDeltaGLMM](https://github.com/nwfsc-assess/geostatistical_delta-GLMM/issues), and [FishStatsUtils](https://github.com/james-thorson/FishStatsUtils/issues), including a search for old and closed issues.
*  *Wiki*:  Users should read and are encouraged to actively contribute to the wiki, which is housed at [the github for SpatialDeltaGLMM](https://github.com/nwfsc-assess/geostatistical_delta-GLMM/wiki)

If there are questions that arise after this, please look for a [VAST Point-of-Contact](https://docs.google.com/spreadsheets/d/1YfYeHHTLwHPxh_5jz4_-hhaRd4gbTq-Cvmii84vxWw0/edit) at your institution and consider contacting them prior to posting an issue.

Database
=============

Regions available in the [example script](https://github.com/james-thorson/VAST/blob/master/examples/Example--simple.R):
![alt text](https://github.com/nwfsc-assess/geostatistical_delta-GLMM/raw/master/examples/global_coverage.png "Global data coverage")
and see [FishViz.org](http://www.FishViz.org) for visualization of results for regions with a public API for their data.

Installation Instructions
=============
[![Build Status](https://travis-ci.org/James-Thorson/VAST.svg?branch=master)](https://travis-ci.org/James-Thorson/VAST)
[![DOI](https://zenodo.org/badge/55002718.svg)](https://zenodo.org/badge/55002718.svg)

This function depends on R version >=3.1.1 and a variety of other tools.

First, install the "devtools" package from CRAN

    # Install and load devtools package
    install.packages("devtools")
    library("devtools")

Next, please install the VAST package from this GitHub repository using a function in the "devtools" package.  This may require using the `INSTALL_opts` option depending upon your version of R:

    # Install package
    install_github("james-thorson/VAST", INSTALL_opts="--no-staged-install")
    # Load package
    library(VAST)

If you are having problems with installation, please consider installing dependencies individually, e.g. using:

    # Install TMB from CRAN
    install.packages("TMB")
    # Install INLA using currently recommended method
    install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
    # Install FishStatsUtils from CRAN
    install_github("james-thorson/FishStatsUtils", INSTALL_opts="--no-staged-install")

Finally, please confirm that VAST is installed by running a model, e.g., following the simple example [here](https://github.com/James-Thorson-NOAA/VAST/wiki/Index-standardization).


Known installation/usage issues
=============

1.  If using a NOAA laptop, sometimes the PATH for Rtools is not correctly specified during installation. In those cases, please follow instructions [here](https://github.com/nwfsc-assess/geostatistical_delta-GLMM/issues/50)

2.  Some versions of R are having problems downloading dependencies from GitHub, see details [here](https://github.com/James-Thorson-NOAA/FishStatsUtils/issues/21)

3.  People using R version 3.6.0 or MRAN 3.5.3 are having a problem with changing standards for package namespaces, see details [here](https://github.com/James-Thorson-NOAA/VAST/issues/189), which appears to be particularly a problem with loading INLA due to install issues with that package.

4.  MacOS users have specific install issues and a discussion of potential fixes is [here](https://github.com/James-Thorson-NOAA/VAST/issues/218#issuecomment-587105809)

5.  MacOS users should be aware that significant speed-ups in model fitting can be accomplished by switching the library used for Basic Linear Algebra Subprograms (BLAS) from the default. There are a few BLAS alternatives available, though, the simplest seems to be using the vecLib library, part of Apple's Accelerate Framework and included in most recent R binaries. To switch the BLAS library, run the following lines in the terminal and then confirm the switch with a call to `sessionInfo()` in R.

```
    # Terminal commands to switch R BLAS library to increase speed
    cd /Library/Frameworks/R.framework/Resources/lib
    ln -sf /System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Versions/Current/libBLAS.dylib libRblas.dylib
```
6.  Windows has a speed-limit on the rate that users can access the GitHub API. You can get around this by installing each package locally from a ZIP file.  You'll need to first download a ZIP file for GitHub repositories `TMBhelper` ([here](https://github.com/kaskr/TMB_contrib_R/tree/master/TMBhelper)), then `ThorsonUtilities` ([here](https://github.com/James-Thorson/utilities)), then `FishStatsUtils` ([here](https://github.com/James-Thorson-NOAA/FishStatsUtils/releases)), then `VAST` ([here](https://github.com/James-Thorson-NOAA/VAST/releases)) to your harddrive in a local directory while recording the directory name (which I will reference as `download_dir`), and then install these packages from each ZIP file in the same order. To install each package, please click "clone or download" -> "Download ZIP" -> `devtools::install_local(path=download_dir, dependencies=FALSE)`

References
=============
### Core functionality
* Thorson, J.T., Barnett, L.A.K., 2017. Comparing estimates of abundance trends and distribution shifts using single- and multispecies models of fishes and biogenic habitat. ICES J. Mar. Sci. 74, 1311–1321. https://doi.org/10.1093/icesjms/fsw193
* Thorson, J.T., 2019. Guidance for decisions using the Vector Autoregressive Spatio-Temporal (VAST) package in stock, ecosystem, habitat and climate assessments. Fish. Res. 210, 143–161. https://doi.org/10.1016/j.fishres.2018.10.013

### Correlated spatio-temporal variation among species (a.k.a. "joint species distribution models")
* Thorson, J.T., Ianelli, J.N., Larsen, E., Ries, L., Scheuerell, M.D., Szuwalski, C., and Zipkin, E. 2016. Joint dynamic species distribution models: a tool for community ordination and spatiotemporal monitoring. Glob. Ecol. Biogeogr. 25(9): 1144–1158. doi:10.1111/geb.12464. url: http://onlinelibrary.wiley.com/doi/10.1111/geb.12464/abstract.
* Thorson, J.T., Scheuerell, M.D., Shelton, A.O., See, K.E., Skaug, H.J., and Kristensen, K. 2015. Spatial factor analysis: a new tool for estimating joint species distributions and correlations in species range. Methods Ecol. Evol. 6(6): 627–637. doi:10.1111/2041-210X.12359. url: http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12359/abstract

### Correlated spatio-temporal variation among years (a.k.a. "Empirical Orthogonal functions")
* Thorson, J.T., Ciannelli, L. and Litzow, M. (In press) Defining indices of ecosystem variability using biological samples of fish communities: a generalization of empirical orthogonal functions. Progress In Oceanography.

### Index of abundance
* Thorson, J.T., Shelton, A.O., Ward, E.J., Skaug, H.J., 2015. Geostatistical delta-generalized linear mixed models improve precision for estimated abundance indices for West Coast groundfishes. ICES J. Mar. Sci. J. Cons. 72(5), 1297–1310. doi:10.1093/icesjms/fsu243. URL: http://icesjms.oxfordjournals.org/content/72/5/1297

### Standardizing samples of size/age-composition data
* Thorson, J. T., and Haltuch, M. A. 2018. Spatio-temporal analysis of compositional data: increased precision and improved workflow using model-based inputs to stock assessment. Canadian Journal of Fisheries and Aquatic Sciences. doi:10.1139/cjfas-2018-0015. URL: http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2018-0015#.W0oloTpKiUk

### Range shift metrics
* Thorson, J.T., Pinsky, M.L., Ward, E.J., 2016. Model-based inference for estimating shifts in species distribution, area occupied, and center of gravity. Methods Ecol. Evol. 7(8), 990-1008.  doi:10.1111/2041-210X.12567.  URL: http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12567/full

### Effective area occupied metric
* Thorson, J.T., Rindorf, A., Gao, J., Hanselman, D.H., and Winker, H. 2016. Density-dependent changes in effective area occupied for sea-bottom-associated marine fishes. Proc R Soc B 283(1840): 20161853. doi:10.1098/rspb.2016.1853. URL: http://rspb.royalsocietypublishing.org/content/283/1840/20161853.

### Spatio-temporal statistical methods
* Thorson, J.T., Skaug, H.J., Kristensen, K., Shelton, A.O., Ward, E.J., Harms, J.H., Benante, J.A., 2014. The importance of spatial models for estimating the strength of density dependence. Ecology 96, 1202–1212. doi:10.1890/14-0739.1. URL: http://www.esajournals.org/doi/abs/10.1890/14-0739.1
* Shelton, A.O., Thorson, J.T., Ward, E.J., Feist, B.E., 2014. Spatial semiparametric models improve estimates of species abundance and distribution. Can. J. Fish. Aquat. Sci. 71, 1655–1666. doi:10.1139/cjfas-2013-0508. URL: http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2013-0508#.VMafDf7F_h4

### Accounting for fish shoals using robust observation models
* Thorson, J. T., I. J. Stewart, and A. E. Punt. 2012. Development and application of an agent-based model to evaluate methods for estimating relative abundance indices for shoaling fish such as Pacific rockfish (Sebastes spp.). ICES Journal of Marine Science 69:635–647. URL: http://icesjms.oxfordjournals.org/content/69/4/635
* Thorson, J. T., I. Stewart, and A. Punt. 2011. Accounting for fish shoals in single- and multi-species survey data using mixture distribution models. Canadian Journal of Fisheries and Aquatic Sciences 68:1681–1693. URL: http://www.nrcresearchpress.com/doi/abs/10.1139/f2011-086#.VMafcf7F_h4

### Accounting for variation among vessels
* Helser, T.E., Punt, A.E., Methot, R.D., 2004. A generalized linear mixed model analysis of a multi-vessel fishery resource survey. Fish. Res. 70, 251–264. doi:10.1016/j.fishres.2004.08.007. url: http://www.sciencedirect.com/science/article/pii/S0165783604001705
* Thorson, J.T., Ward, E.J., 2014. Accounting for vessel effects when standardizing catch rates from cooperative surveys. Fish. Res. 155, 168–176. doi:10.1016/j.fishres.2014.02.036.  url: http://www.sciencedirect.com/science/article/pii/S0165783614000836

### Accounting for fisher targetting in fishery-dependent data
* Thorson, J.T., Fonner, R., Haltuch, M., Ono, K., and Winker, H. In press. Accounting for spatiotemporal variation and fisher targeting when estimating abundance from multispecies fishery data. Can. J. Fish. Aquat. Sci. doi:10.1139/cjfas-2015-0598. URL: http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2015-0598
* Dolder, P.J., Thorson, J.T., Minto, C., 2018. Spatial separation of catches in highly mixed fisheries. Sci. Rep. 8, 13886. https://doi.org/10.1038/s41598-018-31881-w

### Bias-correction of estimated indices of abundance
* Thorson, J.T., and Kristensen, K. 2016. Implementing a generic method for bias correction in statistical models using random effects, with spatial and population dynamics examples. Fish. Res. 175: 66–74. doi:10.1016/j.fishres.2015.11.016. url: http://www.sciencedirect.com/science/article/pii/S0165783615301399

### Estimating and attributing variation in size-structured distribution
* Kai, M., Thorson, J. T., Piner, K. R., and Maunder, M. N. 2017. Spatio-temporal variation in size-structured populations using fishery data: an application to shortfin mako (Isurus oxyrinchus) in the Pacific Ocean. Canadian Journal of Fisheries and Aquatic Sciences. doi:10.1139/cjfas-2016-0327. URL: http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2016-0327#.W0olqjpKiUk.
* Thorson, J. T., Ianelli, J. N., and Kotwicki, S. 2018. The relative influence of temperature and size-structure on fish distribution shifts: A case-study on Walleye pollock in the Bering Sea. Fish and Fisheries. doi:10.1111/faf.12225. URL: https://onlinelibrary.wiley.com/doi/abs/10.1111/faf.12225.

### Estimating fishing impacts using spatial surplus production modelling
* Thorson, J. T., Jannot, J., and Somers, K. 2017. Using spatio-temporal models of population growth and movement to monitor overlap between human impacts and fish populations. Journal of Applied Ecology, 54: 577–587.doi:10.1111/1365-2664.12664. URL: https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/1365-2664.12664

### Estimating species interactions using multispecies Gompertz model
* Thorson, J. T., Munch, S. B., and Swain, D. P. 2017. Estimating partial regulation in spatiotemporal models of community dynamics. Ecology, 98: 1277–1289. doi:10.1002/ecy.1760. URL: https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.1760
* Thorson, J.T., Adams, G., Holsman, K., In press. Spatio-temporal models of intermediate complexity for ecosystem assessments: A new tool for spatial fisheries management. Fish and Fisheries. https://doi.org/10.1111/faf.12398

### Estimating synchrony among species and locations as measure of risk-exposure
* Thorson, J.T., Scheuerell, M.D., Olden, J.D., Schindler, D.E., 2018. Spatial heterogeneity contributes more to portfolio effects than species variability in bottom-associated marine fishes. Proc R Soc B 285, 20180915. https://doi.org/10.1098/rspb.2018.0915

### Forecasting future changes in distribution or abundance
* Thorson, 2019. Forecast skill for predicting distribution shifts:  A retrospective experiment for marine fishes in the Eastern Bering Sea. Fish Fish. 20(1): 159-173. https://doi.org/10.1111/faf.12330

### Combining multiple types of data (e.g., biomass, count, encounter)
* Grüss, A. and Thorson, J.T. (2019) Developing spatio-temporal models using multiple data types for evaluating population trends and habitat usage. ICES Journal of Marine Science 76, 1748–1761. doi:10.1093/icesjms/fsz075.

### Spatially varying coefficient models and their use for fisheries oceanography
* Thorson, J.T., 2019. Measuring the impact of oceanographic indices on species distribution shifts: The spatially varying effect of cold-pool extent in the eastern Bering Sea. Limnol. Oceanogr. 64, 2632–2645. https://doi.org/10.1002/lno.11238

Funding and support for the tool
=============
* Ongoing:  Support from Fisheries Resource Analysis and Monitoring Division (FRAM), Northwest Fisheries Science Center, National Marine Fisheries Service, NOAA
* Ongoing:  Support from Danish Technical University (in particular Kasper Kristensen) for  development of Template Model Builder software, URL: https://www.jstatsoft.org/article/view/v070i05
* Generous support from people knowledgeable about each region and survey! Specific contributions are listed [here](https://github.com/nwfsc-assess/geostatistical_delta-GLMM/blob/master/shiny/Acknowledgements_for_regional_inputs.csv).
* Hoff, G., Thorson, J., and Punt, A.  2018. Spatio-temporal dynamics of groundfish availability to EBS bottom trawl surveys and abundance estimate uncertainties.  North Pacific Research Board (NPRB) 2018 RFP.   
* Rooper, C., Thorson, J., and Boldt, J.  2017. Detecting changes in life history traits and distribution shifts in eastern Bering Sea fishes in response to climate change.  Habitat Information for Stock Assessments 2016 RFP.  
* Ianelli, J., Thorson, J., Kotwicki, S. 2017. Combining acoustic and bottom-trawl data in a spatio-temporal index standardization model for Eastern Bering Sea pollock.  Improve a Stock Assessment 2017 RFP. 
* Thorson, J., Ianelli, J., and O’Brien, L.  2015.  Distribution and application of a new geostatistical index standardization and habitat modeling tool for stock assessments and essential fish habitat designation in Alaska and Northwest Atlantic regions.  Habitat Assessment Improvement Plan 2014 RFP.  URL: https://www.st.nmfs.noaa.gov/ecosystems/habitat/funding/projects/project15-027

=============
### Disclaimer

“The United States Department of Commerce (DOC) GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. DOC has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information. Any claims against the Department of Commerce stemming from the use of its GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.”
