#' Get information for adequate mapping for regions
#'
#' Function provides information to be used when plotting results
#'
#' @export
make_map_info <-
function( Region,
          Extrapolation_List,
          spatial_list = NULL,
          NN_Extrap = spatial_list$PolygonList$NN_Extrap,
          fine_scale = spatial_list$fine_scale,
          Include ){

  # Fix defaults
  if( is.null(fine_scale) ) fine_scale = FALSE
  if( is.null(spatial_list) ){
    if( fine_scale==FALSE ){
      warning("Consider updating inputs to `make_map_info` to enable future use of feature `fine_scale=TRUE`")
    }else{
      stop("Must update inputs to `make_map_info` to enable feature `fine_scale=TRUE`")
    }
  }
  if( missing(Include) ){
    Include = strip_units(Extrapolation_List[["Area_km2_x"]])>0 & strip_units(rowSums(Extrapolation_List[["a_el"]]))>0
  }

  # Initialize
  PlotDF = NULL

  # Loop through regions
  if( tolower(Region)[1]== "california_current" ){
    PlotDF = cbind( Extrapolation_List[["Data_Extrap"]][,c('Lat','Lon')], 'x2i'=NA, 'Include'=Include)
    MappingDetails = list("state", c("alabama","arizona","arkansas","california","colorado","connecticut","delaware","district of columbia","florida","georgia","idaho","illinois","indiana","iowa","kansas","kentucky","louisiana","maine","maryland","massachusetts:martha's vineyard","massachusetts:main","massachusetts:nantucket","michigan:north","michigan:south","minnesota","mississippi","missouri","montana","nebraska","nevada","new hampshire","new jersey","new mexico","new york:manhattan","new york:main","new york:statenisland","new york:longisland","north carolina:knotts","north carolina:main","north carolina:spit","north dakota","ohio","oklahoma","oregon","pennsylvania","rhode island","south carolina","south dakota","tennessee","texas","utah","vermont","virginia:chesapeake","virginia:chincoteague","virginia:main","washington:san juan island","washington:lopez island","washington:orcas island","washington:whidbey island","washington:main","west virginia","wisconsin","wyoming"))
    Xlim=c(-126,-117)
    Ylim=c(32,49)
    #MapSizeRatio = c("Height(in)"=4,"Width(in)"=1.55)
    Rotate = 20     # Degrees counter-clockwise
    Cex = 0.01
    Legend = list(use=TRUE, x=c(65,75), y=c(35,65))
  }
  if( tolower(Region)[1] == "british_columbia" ){
    PlotDF = cbind( Extrapolation_List[["Data_Extrap"]][,c('Lat','Lon')], 'x2i'=NA, 'Include'=Include )
    MappingDetails = list("worldHires", NULL)
    Xlim=c(-133,-126)
    Ylim=c(50,55)
    #MapSizeRatio = c("Height(in)"=2,"Width(in)"=2)
    Rotate = 0
    Cex = 0.1
    Legend = list(use=FALSE,x=c(5,10),y=c(5,45))
  }
  if( tolower(Region)[1] == "eastern_bering_sea" ){
    PlotDF = cbind( Extrapolation_List[["Data_Extrap"]][,c('Lat','Lon')], 'x2i'=NA, 'Include'=Include )
    MappingDetails = list("worldHires", NULL)
    Xlim = c(-180,-158)
    Ylim=c(54,63)
    #MapSizeRatio = c("Height(in)"=4,"Width(in)"=5)
    Rotate = 0
    Cex = 0.01
    Legend = list(use=TRUE,x=c(81,86),y=c(48,88))
  }
  if( tolower(Region)[1] == "aleutian_islands" ){
    PlotDF = cbind( Extrapolation_List[["Data_Extrap"]][,c('Lat','Lon')], 'x2i'=NA, 'Include'=Include )
    PlotDF[,'Lon'] = PlotDF[,'Lon'] %% 360 # Change units to match world2Hires
    MappingDetails = list("world2Hires", NULL)
    Xlim = c(170,195)
    Ylim=c(51,55)
    #MapSizeRatio = c("Height(in)"=2,"Width(in)"=5)
    Rotate = 0
    Cex = 0.2
    Legend = list(use=FALSE,x=c(5,10),y=c(5,45))
  }
  if( tolower(Region)[1] == "gulf_of_alaska" ){
    PlotDF = cbind( Extrapolation_List[["Data_Extrap"]][,c('Lat','Lon')], 'x2i'=NA, 'Include'=Include )
    MappingDetails = list("world", NULL)
    Xlim = c(-171,-132)
    Ylim=c(52,61)
    #MapSizeRatio = c("Height(in)"=2.5,"Width(in)"=6)
    Rotate = 0
    Cex = 0.01
    Legend = list(use=TRUE,x=c(5,10),y=c(30,65))
  }
  if( tolower(Region)[1] == "northwest_atlantic" ){
    PlotDF = cbind( Extrapolation_List[["Data_Extrap"]][,c('Lat','Lon')], 'x2i'=NA, 'Include'=Include )
    #MappingDetails = list("world", NULL)
    MappingDetails = list("state", c("alabama","arizona","arkansas","california","colorado","connecticut","delaware","district of columbia","florida","georgia","idaho","illinois","indiana","iowa","kansas","kentucky","louisiana","maine","maryland","massachusetts:martha's vineyard","massachusetts:main","massachusetts:nantucket","michigan:north","michigan:south","minnesota","mississippi","missouri","montana","nebraska","nevada","new hampshire","new jersey","new mexico","new york:manhattan","new york:main","new york:statenisland","new york:longisland","north carolina:knotts","north carolina:main","north carolina:spit","north dakota","ohio","oklahoma","oregon","pennsylvania","rhode island","south carolina","south dakota","tennessee","texas","utah","vermont","virginia:chesapeake","virginia:chincoteague","virginia:main","washington:san juan island","washington:lopez island","washington:orcas island","washington:whidbey island","washington:main","west virginia","wisconsin","wyoming"))
    Xlim = c(-80,-65)
    Ylim=c(32,45)
    #MapSizeRatio = c("Height(in)"=4,"Width(in)"=3)
    Rotate = 0
    Cex = 0.01
    Legend = list(use=TRUE,x=c(75,80),y=c(5,35))
  }
  if( tolower(Region)[1] == "south_africa" ){
    PlotDF = cbind( Extrapolation_List[["Data_Extrap"]][,c('Lat','Lon')], 'x2i'=NA, 'Include'=Include )
    MappingDetails = list("worldHires", NULL )
    Xlim = c(14,26)
    Ylim=c(-37,-28)
    #MapSizeRatio = c("Height(in)"=4,"Width(in)"=3)
    Rotate = 0
    Cex = 0.1
    Legend = list(use=FALSE,x=c(5,10),y=c(4,45))
  }
  if( tolower(Region)[1] == "gulf_of_st_lawrence" ){
    PlotDF = cbind( Extrapolation_List[["Data_Extrap"]][,c('Lat','Lon')], 'x2i'=NA, 'Include'=Include )
    MappingDetails = list("worldHires", "Canada" )
    Xlim = range(Extrapolation_List[["Data_Extrap"]][which(strip_units(Extrapolation_List[["Area_km2_x"]])>0),'Lon'])
    Ylim = range(Extrapolation_List[["Data_Extrap"]][which(strip_units(Extrapolation_List[["Area_km2_x"]])>0),'Lat'])
    #MapSizeRatio = c("Height(in)"=4,"Width(in)"=4)
    Rotate = 0
    Cex = 1.0
    Legend = list(use=FALSE,x=c(5,10),y=c(4,45))
  }
  if( tolower(Region)[1] == "new_zealand" ){
    PlotDF = cbind( Extrapolation_List[["Data_Extrap"]][,c('Lat','Lon')], 'x2i'=NA, 'Include'=Include )
    MappingDetails = list("worldHires", NULL )
    Xlim=c(172,187)
    Ylim=c(-46,-41)
    #MapSizeRatio = c("Height(in)"=2,"Width(in)"=5)
    Rotate = 0     # Degrees counter-clockwise
    Cex = 0.01
    Legend = list(use=FALSE,x=c(5,10),y=c(5,45))
  }
  if( tolower(Region)[1] == "habcam" ){
    PlotDF = cbind( Extrapolation_List[["Data_Extrap"]][,c('Lat','Lon')], 'x2i'=NA, 'Include'=Include)
    #MappingDetails = list("state", c("alabama","arizona","arkansas","california","colorado","connecticut","delaware","district of columbia","florida","georgia","idaho","illinois","indiana","iowa","kansas","kentucky","louisiana","maine","maryland","massachusetts:martha's vineyard","massachusetts:main","massachusetts:nantucket","michigan:north","michigan:south","minnesota","mississippi","missouri","montana","nebraska","nevada","new hampshire","new jersey","new mexico","new york:manhattan","new york:main","new york:statenisland","new york:longisland","north carolina:knotts","north carolina:main","north carolina:spit","north dakota","ohio","oklahoma","oregon","pennsylvania","rhode island","south carolina","south dakota","tennessee","texas","utah","vermont","virginia:chesapeake","virginia:chincoteague","virginia:main","washington:san juan island","washington:lopez island","washington:orcas island","washington:whidbey island","washington:main","west virginia","wisconsin","wyoming"))
    MappingDetails = list("worldHires", NULL )
    Xlim = range(Extrapolation_List[["Data_Extrap"]][which(strip_units(Extrapolation_List[["Area_km2_x"]])>0),'Lon'])
    Ylim = range(Extrapolation_List[["Data_Extrap"]][which(strip_units(Extrapolation_List[["Area_km2_x"]])>0),'Lat'])
    #MapSizeRatio = c("Height(in)"=4,"Width(in)"=3)
    Rotate = 20     # Degrees counter-clockwise
    Cex = 0.01
    Legend = list(use=TRUE,x=c(70,90),y=c(5,35))
  }
  if( tolower(Region)[1] == "gulf_of_mexico" ){
    PlotDF = cbind( Extrapolation_List[["Data_Extrap"]][,c('Lat','Lon')], 'x2i'=NA, 'Include'=Include)
    #MappingDetails = list("state", c("alabama","arizona","arkansas","california","colorado","connecticut","delaware","district of columbia","florida","georgia","idaho","illinois","indiana","iowa","kansas","kentucky","louisiana","maine","maryland","massachusetts:martha's vineyard","massachusetts:main","massachusetts:nantucket","michigan:north","michigan:south","minnesota","mississippi","missouri","montana","nebraska","nevada","new hampshire","new jersey","new mexico","new york:manhattan","new york:main","new york:statenisland","new york:longisland","north carolina:knotts","north carolina:main","north carolina:spit","north dakota","ohio","oklahoma","oregon","pennsylvania","rhode island","south carolina","south dakota","tennessee","texas","utah","vermont","virginia:chesapeake","virginia:chincoteague","virginia:main","washington:san juan island","washington:lopez island","washington:orcas island","washington:whidbey island","washington:main","west virginia","wisconsin","wyoming"))
    MappingDetails = list("worldHires", NULL )
    Xlim = range(Extrapolation_List[["Data_Extrap"]][which(strip_units(Extrapolation_List[["Area_km2_x"]])>0),'Lon'])
    Ylim = range(Extrapolation_List[["Data_Extrap"]][which(strip_units(Extrapolation_List[["Area_km2_x"]])>0),'Lat'])
    #MapSizeRatio = c("Height(in)"=4,"Width(in)"=3)
    Rotate = 0     # Degrees counter-clockwise
    Cex = 0.60
    Legend = list(use=TRUE,x=c(81,86),y=c(48,83))
  }
  if( is.null(PlotDF) ){
    PlotDF = cbind( Extrapolation_List[["Data_Extrap"]][,c('Lat','Lon')], 'x2i'=NA, 'Include'=Include )
    MappingDetails = list("worldHires", NULL )
    Xlim = range(Extrapolation_List[["Data_Extrap"]][which(strip_units(Extrapolation_List[["Area_km2_x"]])>0),'Lon'])
    Ylim = range(Extrapolation_List[["Data_Extrap"]][which(strip_units(Extrapolation_List[["Area_km2_x"]])>0),'Lat'])
    Rotate = 0
    Cex = 0.1
    Legend = list(use=FALSE,x=c(5,10),y=c(5,45))
  }

  #
  #if( fine_scale==TRUE | spatial_list$Method=="Stream_network" ){
  #  PlotDF[,'x2i'] = NA
  #  PlotDF[ which(Extrapolation_List[["Area_km2_x"]]>0),'x2i'] = 1:length(which(Extrapolation_List[["Area_km2_x"]]>0))
  #}else{
  #  PlotDF[,'x2i'] = NN_Extrap$nn.idx[,1]
  #}
  ## Exceptions
  #if( tolower(Region)[1] == "eastern_bering_sea" ){
  #  PlotDF = PlotDF[which(PlotDF[,'Lon']<0),]
  #}
  PlotDF[,'x2i'] = spatial_list$g_e

  # Plotting zone has to match zone used for Extrapolation_List, for COG estimates (in UTM via Z_xm) to match map UTM
  if( is.numeric(Extrapolation_List$zone) ){
    Zone = Extrapolation_List[["zone"]] - ifelse(Extrapolation_List$flip_around_dateline==TRUE, 30, 0)
  }else{
    Zone = Extrapolation_List$zone
  }

  # Determine map size (equal distance along x-axis and y-axis)
  if( all(!is.na(Extrapolation_List$Data_Extrap[,c('N_km','E_km')])) ){
    MapSizeRatio = c("Height(in)"=diff(range(Extrapolation_List$Data_Extrap[,'N_km'])) , "Width(in)"=diff(range(Extrapolation_List$Data_Extrap[,'E_km'])) )
  }else{
    MapSizeRatio = c("Height(in)"=diff(range(Extrapolation_List$Data_Extrap[,'Lat'])) , "Width(in)"=diff(range(Extrapolation_List$Data_Extrap[,'Lon'])) )
  }
  MapSizeRatio = MapSizeRatio / sqrt(prod(MapSizeRatio)) * 4  # 4 square-inches

  # bundle and return
  mapdetails_list = list("PlotDF"=PlotDF, "MappingDetails"=MappingDetails, "Xlim"=Xlim, "Ylim"=Ylim,
    "MapSizeRatio"=MapSizeRatio, "Rotate"=Rotate, "Cex"=Cex, "Legend"=Legend, "Zone"=Zone, "fine_scale"=fine_scale )
  return( mapdetails_list )
}
