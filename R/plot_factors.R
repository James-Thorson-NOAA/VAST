
#' Plot factor-decomposition of covariance
#'
#' \code{plot_factors} plots factor loadings, average spatial factors, and spatio-temporal factors
#'
#' @inheritParams plot_overdispersion
#' @inheritParams summarize_covariance
#' @inheritParams rotate_factors
#' @inheritParams plot_maps
#' @param mapdetails_list output from \code{\link{make_map_info}}
#' @param Dim_year Plotting dimension (row,column) for plot of years (default: square with sufficient size for number of years)
#' @param Dim_species Plotting dimension (row,column) for plot of categories (default: square with sufficient size for number of categories)
#' @param plotdir directory for saving plots
#' @param land_color color for filling in land (use \code{land_color=rgb(0,0,0,alpha=0)} for transparent land)
#' @param factor_samples number of samples from joint precision used to estimate precision for derived quantities
#' @param ... additional arguments passed to \code{\link{plot_maps}} and/or \code{\link{plot_variable}} when plotting factor-values on a map
#'
#' @references For details regarding spatial factor analysis see \url{https://doi.org/10.1111/2041-210X.12359}
#' @references For details regarding multi-species ordination see \url{https://doi.org/10.1016/j.fishres.2018.10.013}
#' @export
plot_factors <-
function( fit,
          Report = fit$Report,
          ParHat = fit$ParHat,
          Data = fit$data_list,
          Obj = fit$tmb_list$Obj,
          SD = fit$parameter_estimates$SD,
          year_labels = fit$year_labels,
          category_names = fit$category_names,
          RotationMethod = "PCA",
          mapdetails_list = NULL,
          Dim_year = NULL,
          Dim_species = NULL,
          projargs = '+proj=longlat',
          plotdir = getwd(),
          land_color = "grey",
          zlim = NA,
          testcutoff = 1e-4,
          factor_samples = 100,
          plot_value = "estimate",
          ... ){

  #
  if(is.null(mapdetails_list)) message( "`plot_factors(.) skipping plots because argument `mapdetails_list` is missing")

  # Extract Options and Options_vec (depends upon version)
  if( all(c("Options","Options_vec") %in% names(Data)) ){
    Options_vec = Data$Options_vec
    Options = Data$Options
  }
  if( "Options_list" %in% names(Data) ){
    Options_vec = Data$Options_list$Options_vec
    Options = Data$Options_list$Options
  }

  #### Deals with backwards compatibility for FieldConfig
  # Converts from 4-vector to 3-by-2 matrix
  if( is.vector(Data$FieldConfig) && length(Data$FieldConfig)==4 ){
    Data$FieldConfig = rbind( matrix(Data$FieldConfig,ncol=2,dimnames=list(c("Omega","Epsilon"),c("Component_1","Component_2"))), "Beta"=c(-2,-2) )
  }
  # Converts from 3-by-2 matrix to 4-by-2 matrix
  if( is.matrix(Data$FieldConfig) & all(dim(Data$FieldConfig)==c(3,2)) ){
    Data$FieldConfig = rbind( Data$FieldConfig, "Epsilon_time"=c(-3,-3) )
  }
  # Checks for errors
  if( !is.matrix(Data$FieldConfig) || !all(dim(Data$FieldConfig)==c(4,2)) ){
    stop("`FieldConfig` has the wrong dimensions in `plot_factors`")
  }
  # Renames
  dimnames(Data$FieldConfig) = list( c("Omega","Epsilon","Beta","Epsilon_time"), c("Component_1","Component_2") )

  # Fill in missing inputs
  # COULD REPLACE WITH amend_output, except missing `extrapolation_list`
  #if( "D_gct" %in% names(Report) ){
  #  if( is.null(year_labels) ) year_labels = 1:dim(Report$D_gct)[3]
  #  if( is.null(category_names) ) category_names = 1:dim(Report$D_gct)[2]
  #}
  #if( "D_xcy" %in% names(Report) ){
  #  if( is.null(year_labels) ) year_labels = 1:dim(Report$D_xcy)[3]
  #  if( is.null(category_names) ) category_names = 1:dim(Report$D_xcy)[2]
  #}
  #if( "D_gcy" %in% names(Report) ){
  #  if( is.null(year_labels) ) year_labels = 1:dim(Report$D_gcy)[3]
  #  if( is.null(category_names) ) category_names = 1:dim(Report$D_gcy)[2]
  #}
  # Overwrite labels using run-time user inputs if provided
  Report = amend_output( fit = fit,
                         TmbData = Data,
                         Report = Report,
                         Sdreport = SD,
                         year_labels = year_labels,
                         category_names = category_names )

  # Dimensions for plotting
  Dim = function( num ) c(ceiling(sqrt(num)), ceiling(num/ceiling(sqrt(num))) )
  Dim_year = Dim(length(year_labels))
  Dim_species = Dim(length(category_names))

  # Extract loadings matrices (more numerically stable than extracting covariances, and then re-creating Cholesky)
  Psi2prime_list = Psiprime_list = Psi2prime_SE_list = Lprime_SE_list = Hinv_list = Lprime_list = L_list = vector("list", length=8)    # Add names at end so that NULL doesn't interfere

  # Loop through
  for(i in 1:8){

    # Variable names
    Par_name = c("Omega1", "Epsilon1", "Beta1", "EpsilonTime1", "Omega2", "Epsilon2", "Beta2", "EpsilonTime2")[i]
    Lpar_name = c("L_omega1_z", "L_epsilon1_z", "L_beta1_z", "Ltime_epsilon1_z", "L_omega2_z", "L_epsilon2_z", "L_beta2_z", "Ltime_epsilon2_z")[i]

    # Backwards compatible loading of variables and names
    if(Par_name == "Omega1"){ Var_name = "Omegainput1_sf"; Var2_name = "Omegainput1_gf"; L_name = "L_omega1_cf" }
    if(Par_name == "Epsilon1"){ Var_name = "Epsiloninput1_sft"; Var2_name = "Epsiloninput1_gft"; L_name = "L_epsilon1_cf" }
    if(Par_name == "Beta1"){ Var_name = "beta1_ft"; Var2_name = "missing"; L_name = "L_beta1_cf" }
    if(Par_name == "EpsilonTime1"){ Var_name = "Epsiloninput1_sff"; Var2_name = "Epsiloninput1_gff"; L_name = "Ltime_epsilon1_tf" }
    if(Par_name == "Omega2"){ Var_name = "Omegainput2_sf"; Var2_name = "Omegainput2_gf"; L_name = "L_omega2_cf" }
    if(Par_name == "Epsilon2"){ Var_name = "Epsiloninput2_sft"; Var2_name = "Epsiloninput2_gft"; L_name = "L_epsilon2_cf" }
    if(Par_name == "Beta2"){ Var_name = "beta2_ft"; Var2_name = "missing"; L_name = "L_beta2_cf" }
    if(Par_name == "EpsilonTime2"){ Var_name = "Epsiloninput2_sff"; Var2_name = "Epsiloninput2_gff"; L_name = "Ltime_epsilon2_tf" }

    # Continue if component is included
    if( as.vector(Data[["FieldConfig"]])[i] > 1 ){

      # Get loadings matrix
      if( L_name %in% names(Report) ){
        L_list[[i]] = Report[[L_name]]
      }else{
        L_list[[i]] = calc_cov( L_z=ParHat[[Lpar_name]], n_f=as.vector(Data[["FieldConfig"]])[i], n_c=Data$n_c, returntype="loadings_matrix" )
      }
      if( !Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
        rownames(L_list[[i]]) = category_names
      }

      # Get covariance
      if(Var_name%in%names(ParHat)) Psi_sjt = ParHat[[Var_name]]
      if(Var_name%in%names(Report)) Psi_sjt = Report[[Var_name]]
      Psi_gjt = Report[[Var2_name]]
      # Get tau
      logkappa = unlist(ParHat[c('logkappa1','logkappa2')])[c(1,1,1,1,2,2,2,2)[i]]
      if(Options_vec[8]==0){
        tau = 1 / (exp(logkappa) * sqrt(4*pi));
      }else if(Options_vec[8]==1){
        tau = 1 / sqrt(1-exp(logkappa*2));
      }else stop("Check 'Options_vec[8]' for allowable entries")
      ## the betas and EpsilonTimeare transposed compared to others so fix that here
      if( Par_name %in% c("Beta1","Beta2") ){
        Psi_sjt = t(Psi_sjt)
        tau = 1
      }
      if( Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
        Psi_sjt = aperm( Psi_sjt, c(1,3,2) )
        Psi_gjt = aperm( Psi_gjt, c(1,3,2) )
      }
      if(is.null(Psi_sjt)){
        stop(paste("Covariance is empty for parameter", Var_name))
      }

      # Rotate stuff
      Var_rot = rotate_factors( L_pj = L_list[[i]],
                Psi_sjt = Psi_sjt/tau,
                RotationMethod = RotationMethod,
                testcutoff = testcutoff,
                quiet = TRUE )
      if( Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
        Var_rot$Psi_rot = aperm( Var_rot$Psi_rot, c(1,3,2) )
      }
      Report_tmp = list("D_gct"=Var_rot$Psi_rot, "Epsilon1_gct"=Var_rot$Psi_rot, "Epsilon2_gct"=Var_rot$Psi_rot)
      Lprime_list[[i]] = Var_rot$L_pj_rot
      if( !Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
        rownames(Lprime_list[[i]]) = category_names
      }
      Psiprime_list[[i]] = Var_rot$Psi_rot
      Hinv_list[[i]] = Var_rot$Hinv

      # Extract SEs if available
      # Could also edit to extract L_SE_list and Psi2_SE_list
      if( !is.null(Obj) && class(SD)=="sdreport" ){
        L_cfr = sample_variable( Sdreport=SD, Obj=Obj, variable_name=L_name, n_samples=factor_samples, sample_fixed=TRUE, seed=123456 )
        Psi_gjtr = sample_variable( Sdreport=SD, Obj=Obj, variable_name=Var2_name, n_samples=factor_samples, sample_fixed=TRUE, seed=123456 )
        if( Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
          Psi_gjtr = aperm( Psi_gjtr, c(1,3,2,4) )
        }
        Lprime_cfr = array(NA, dim=dim(L_cfr) )
        Psiprime_gjtr = array(NA, dim=dim(Psi_gjtr) )
        for( rI in seq_len(factor_samples) ){
          tmplist = rotate_factors( L_pj = array(L_cfr[,,rI],dim=dim(L_cfr)[1:2]),
            Psi_sjt = array(Psi_gjtr[,,,rI],dim=dim(Psi_gjtr)[1:3])/tau,
            RotationMethod = RotationMethod,
            testcutoff = testcutoff,
            quiet = TRUE )
          Lprime_cfr[,,rI] = tmplist$L_pj_rot
          Psiprime_gjtr[,,,rI] = tmplist$Psi_rot
        }
        if( Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
          Psiprime_gjtr = aperm( Psiprime_gjtr, c(1,3,2,4) )
        }
        Lprime_SE_list[[i]] = apply(Lprime_cfr, MARGIN=1:2, FUN=sd)
        if( nrow(Lprime_SE_list[[i]])==length(category_names) ){
          rownames(Lprime_SE_list[[i]]) = category_names
        }
        Psi2prime_SE_list[[i]] = apply(Psiprime_gjtr, MARGIN=1:3, FUN=sd)
      }

      # Extract projected factors is available
      if( !is.null(Psi_gjt) ){
        Var2_rot = rotate_factors( L_pj = L_list[[i]],
          Psi_sjt = Psi_gjt/tau,
          RotationMethod = RotationMethod,
          testcutoff = testcutoff,
          quiet = TRUE )
        if( Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
          Var2_rot$Psi_rot = aperm( Var2_rot$Psi_rot, c(1,3,2) )
        }
        Report2_tmp = list("D_gct"=Var2_rot$Psi_rot, "Epsilon1_gct"=Var2_rot$Psi_rot, "Epsilon2_gct"=Var2_rot$Psi_rot)
        Psi2prime_list[[i]] = Var2_rot$Psi_rot
      }else{
        Report2_tmp = NULL
      }

      # Plot loadings
      Dim_factor = Dim( as.vector(Data[["FieldConfig"]])[i] )
      png( file=file.path(plotdir,paste0("Factor_loadings--",Par_name,".png")), width=Dim_factor[2]*4, height=Dim_factor[1]*4, units="in", res=200 )
        par( mfrow=Dim_factor, mar=c(2,2,1,0), oma=c(0,0,0,0), mgp=c(2,0.5,0), tck=-0.02 )
        for( cI in 1:as.vector(Data[["FieldConfig"]])[i] ){
          if( Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
            if(any(is.na(as.numeric(year_labels)))) stop("Check `At` in `plot_factors(.)`")
            plot_loadings( L_pj=Lprime_list[[i]], Lsd_pj=Lprime_SE_list[[i]], whichfactor=cI, At=as.numeric(year_labels), LabelPosition="Side" )
          }else{
            plot_loadings( L_pj=Lprime_list[[i]], Lsd_pj=Lprime_SE_list[[i]], whichfactor=cI, At=1:nrow(Var_rot$L_pj_rot) )
          }
        }
      dev.off()

      # Plot Beta
      if( Par_name %in% c("Beta1","Beta2") ){
        png(file=file.path(plotdir,paste0("factor_values--",Par_name,".png")),
            width=4, height=4, res=200, units='in')
        matplot( x = as.numeric(year_labels),
                 y = array(Psiprime_list[[i]][,,1],dim=dim(Psiprime_list[[i]])[1:2]),
                 type = "l",
                 col = rainbow(ncol(Psiprime_list[[i]])),
                 lty = "solid",
                 xlab = "Time",
                 ylab = "" )
        legend( "topleft", bty = "n", fill = rainbow(ncol(Psiprime_list[[i]])), legend = 1:ncol(Psiprime_list[[i]]))
        dev.off()
      }

      # Plot factors
      if( !is.null(mapdetails_list) & !is.null(Report2_tmp) ){

        # Plot Epsilon
        # Use plot_maps to automatically make one figure per factor
        if( Par_name %in% c("Epsilon1","Epsilon2") ){
          factor_setting = Data$FieldConfig['Epsilon',match(Par_name,c("Epsilon1","Epsilon2"))]
          if( factor_setting==(-3) & dim(Var_rot$Psi_rot)[2]==length(category_names) ){
            factor_names = category_names
          }else{
            factor_names = paste0("Factor_",1:dim(Var_rot$Psi_rot)[2])
          }
          plot_maps( plot_set=c(6,6,NA,6,7,7,NA,7)[i],
                     fit = fit,
                     Report = Report2_tmp,
                     PlotDF = mapdetails_list[["PlotDF"]],
                     MapSizeRatio = mapdetails_list[["MapSizeRatio"]],
                     working_dir = plotdir,
                     year_labels = year_labels,
                     category_names = factor_names,
                     legend_x = mapdetails_list[["Legend"]]$x/100,
                     legend_y = mapdetails_list[["Legend"]]$y/100,
                     zlim = zlim,
                     land_color = land_color,
                     n_cells = NULL,
                     projargs = projargs,
                    # ...,
                     plot_value = "estimate" )
        }  #

        # Plot Omega
        # Use plot_variable to plot all factors on single figure
        if( Par_name %in% c("Omega1", "Omega2")){
          plot_variable( Y_gt = array(Report2_tmp$D_gct[,,1],dim = dim(Report2_tmp$D_gct)[1:2]),
                         map_list = mapdetails_list,
                         working_dir = plotdir,
                         panel_labels = paste0("Factor_",1:dim(Var_rot$Psi_rot)[2]),
                         file_name = paste0("Factor_maps--",Par_name),
                         land_color = land_color,
                        # ...,
                         projargs = projargs )
        }

        ## Doesn't make sense to make maps of beta factors since they aren't spatial

        # Plot EpsilonTime
        if( Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
          if( dim(Var_rot$Psi_rot)[2] == length(category_names) ){
            factor_names = category_names
          }else{
            factor_names = paste0("Factor_",1:dim(Var_rot$Psi_rot)[2])
          }
          tmp_names = paste0("Factor_",1:dim(Var_rot$Psi_rot)[3])
          plot_maps( plot_set=c(6,6,NA,6,7,7,NA,7)[i],
                     fit = fit,
                     Report = Report2_tmp,
                     PlotDF = mapdetails_list[["PlotDF"]],
                     MapSizeRatio = mapdetails_list[["MapSizeRatio"]],
                     working_dir = plotdir,
                     category_names = factor_names,
                     #year_labels = tmp_names,
                     Panel = "Year",
                     legend_x = mapdetails_list[["Legend"]]$x/100,
                     legend_y = mapdetails_list[["Legend"]]$y/100,
                     zlim = zlim,
                     n_cells = NULL,
                     land_color = land_color,
                     projargs = projargs,
                    # ...,
                     plot_value = "estimate" )
        }  #
      }
    }else{
      Psi2prime_SE_list[[i]] = Psiprime_list[[i]] = Psi2prime_list[[i]] = Lprime_SE_list[[i]] = Hinv_list[[i]] = Lprime_list[[i]] = L_list[[i]] = "Element not estimated, and therefore empty"
    }
  }

  # Return stuff invisibly
  names(Hinv_list) = names(Psi2prime_list) = names(Psiprime_list) = names(Lprime_SE_list) = names(Lprime_list) = names(L_list) = c("Omega1", "Epsilon1", "Beta1", "EpsilonTime1", "Omega2", "Epsilon2", "Beta2", "EpsilonTime2")
  Return = list("Loadings"=L_list, "Rotated_loadings"=Lprime_list, "Rotated_factors"=Psiprime_list, "Rotated_projected_factors"=Psi2prime_list, "Rotation_matrices"=Hinv_list)
  if( !is.null(Obj) && class(SD)=="sdreport" ){
    Return[["Rotated_loadings_SE"]] = Lprime_SE_list
    Return[["Rotated_projected_factors_SE"]] = Psi2prime_SE_list
  }
  return( invisible(Return) )
}
