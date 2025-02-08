
#' Make mesh for distances among points
#'
#' \code{make_mesh} builds a tagged list representing distances for isotropic or geometric anisotropic triangulated mesh
#'
#' @param loc_x location (eastings and northings in kilometers, UTM) for each sample or knot
#' @param Method spatial method determines ("Mesh" and "Grid" give
#' @param anisotropic_mesh OPTIONAL, anisotropic mesh (if missing, its recalculated from loc_x)
#' @param ... Arguments passed to \code{INLA::inla.mesh.create}

#' @return Tagged list containing distance metrics

#' @export
make_mesh <-
function( loc_x,
          loc_g,
          loc_i,
          Method,
          Extrapolation_List,
          anisotropic_mesh = NULL,
          fine_scale = FALSE,
          map_data,
          mesh_package = c("INLA","fmesher"),
          ...){

  #######################
  # Create the anisotropic SPDE mesh using 2D coordinates
  #######################
  mesh_package = match.arg( mesh_package )

  # 2D coordinates SPDE
  if( is.null(anisotropic_mesh)){
    if( fine_scale==FALSE ){
      if( mesh_package=="INLA" ){
        if(!require(INLA)) stop("Please install INLA or switch to `settings$mesh_package=`fmesher`")
        anisotropic_mesh = INLA::inla.mesh.create( loc_x, plot.delay=NULL, ...)
      }else{
        anisotropic_mesh = fm_mesh_2d( loc_x, ... )
      }
    }else{
      if( mesh_package=="INLA" ){
        loc_z = rbind( loc_x, loc_g, loc_i )
        if(!require(INLA)) stop("Please install INLA or switch to `settings$mesh_package=`fmesher`")
        outer_hull = INLA::inla.nonconvex.hull( as.matrix(loc_z), convex = -0.05, concave = -0.05)
        anisotropic_mesh = INLA::inla.mesh.create( loc_x, plot.delay=NULL, boundary=outer_hull, ...)
      }else{
        #outer_hull = fm_nonconvex_hull( as.matrix(loc_z), convex = -0.05, concave = -0.05 )
        anisotropic_mesh = fm_mesh_2d( loc_x, ... ) # , boundary=outer_hull, ... )
      }
    }
  }

  #anisotropic_spde = INLA::inla.spde2.matern(anisotropic_mesh, alpha=2)
  anisotropic_spde = fm_fem( anisotropic_mesh, order=2 )
  anisotropic_spde$param.inla = list( "M0" = anisotropic_spde$c0,
                                      "M1" = anisotropic_spde$g1,
                                      "M2" = anisotropic_spde$g2)
  anisotropic_spde$mesh = anisotropic_mesh
  anisotropic_spde$n.spde = anisotropic_mesh$n

  # Pre-processing in R for anisotropy
  Dset = 1:2
  # Triangle info
  TV = anisotropic_mesh$graph$tv       # Triangle to vertex indexing
  V0 = anisotropic_mesh$loc[TV[,1],Dset]   # V = vertices for each triangle
  V1 = anisotropic_mesh$loc[TV[,2],Dset]
  V2 = anisotropic_mesh$loc[TV[,3],Dset]
  E0 = V2 - V1                      # E = edge for each triangle
  E1 = V0 - V2
  E2 = V1 - V0
  
  # Pre-processing for barriers
  # Barriers don't affect projection matrix A
  # Obtain polygon for water
  if( missing(map_data) ){
    map_data = rnaturalearth::ne_countries( scale=50, returnclass="sf" )
    #map_data = sf::as_Spatial( map_data )
    #attr(map_data,"proj4string") = sp::CRS("+proj=longlat +datum=WGS84")
    #proj4string(map_data) = sp::CRS("+proj=longlat +datum=WGS84")
  }

  # Calculate centroid of each triangle in mesh and convert to SpatialPoints
  n_triangles = length(anisotropic_mesh$graph$tv[,1])
  posTri = matrix(NA, nrow=n_triangles, ncol=2)
  for(tri_index in 1:n_triangles){
    temp = anisotropic_mesh$loc[ anisotropic_mesh$graph$tv[tri_index,], ]
    posTri[tri_index,] = colMeans(temp)[c(1,2)]
  }
  #posTri = sp::SpatialPoints(posTri, proj4string=sp::CRS(Extrapolation_List$projargs) )
  #posTri = sp::spTransform(posTri, CRSobj=sp::CRS(map_data) )
  posTri = sf::st_sfc( sf::st_multipoint(posTri), crs=sf::st_crs(Extrapolation_List$projargs) )
  posTri = sf::st_transform( posTri, crs=sf::st_crs(map_data) )

  # Calculate set of triangles barrier.triangles with centroid over land
  if( Method == "Barrier" ){
    #anisotropic_mesh_triangles_over_land = unlist(sp::over(map_data, posTri, returnList=TRUE))
    anisotropic_mesh_triangles_over_land = unlist(sf::st_intersects(map_data, posTri))
  }else{
    anisotropic_mesh_triangles_over_land = vector()
  }
  #
  #plot( x=posTri@coords[,1], y=posTri@coords[,2], col=ifelse(1:n_triangles%in%triangles_over_land,"black","red") )

  # Create Barrier object if requested
    # Don't do this unless necessary, because it sometimes throws an error
  #Diagnose issues:  assign("anisotropic_mesh", anisotropic_mesh, envir = .GlobalEnv)
  barrier_finite_elements = inla.barrier.fem.copy(mesh=anisotropic_mesh,
    barrier.triangles=anisotropic_mesh_triangles_over_land)
  barrier_list = list(C0 = barrier_finite_elements$C[[1]],
    C1 = barrier_finite_elements$C[[2]],
    D0 = barrier_finite_elements$D[[1]],
    D1 = barrier_finite_elements$D[[2]],
    I = barrier_finite_elements$I )
  # sp::plot( INLA::inla.barrier.polygon(anisotropic_mesh, triangles_over_land) )

  # Calculate Areas 
  crossprod_fn = function(Vec1,Vec2) abs(det( rbind(Vec1,Vec2) ))
  Tri_Area = rep(NA, nrow(E0))
  for(i in 1:length(Tri_Area)) Tri_Area[i] = crossprod_fn( E0[i,],E1[i,] )/2   # T = area of each triangle

  ################
  # Add the isotropic SPDE mesh for spherical or 2D projection, depending upon `Method` input
  ################

  # Mesh and SPDE for different inputs
  if(Method %in% c("Mesh","Grid","Stream_network","Barrier")){
    loc_isotropic_mesh = loc_x
    isotropic_mesh = anisotropic_mesh
  }
  if(Method %in% c("Spherical_mesh")){
    #loc_isotropic_mesh = INLA::inla.mesh.map(loc_x, projection="longlat", inverse=TRUE) # Project from lat/long to mesh coordinates
    #isotropic_mesh = INLA::inla.mesh.create( loc_isotropic_mesh, plot.delay=NULL, ...)
    stop("Method `Spherical_mesh` no longer supported")
  }
  #isotropic_spde = INLA::inla.spde2.matern(isotropic_mesh, alpha=2)
  isotropic_spde = fm_fem(isotropic_mesh, order=2)
  isotropic_spde$param.inla = list( "M0" = isotropic_spde$c0,
                                    "M1" = isotropic_spde$g1,
                                    "M2" = isotropic_spde$g2)
  isotropic_spde$mesh = isotropic_mesh
  isotropic_spde$n.spde = isotropic_mesh$n

  ####################
  # Return stuff
  ####################
  #if( isotropic_mesh$n != anisotropic_mesh$n ) stop("Check `Calc_Anisotropic_Mesh` for problem")

  Return = list( "loc_x"=loc_x,
                "loc_isotropic_mesh"=loc_isotropic_mesh,
                "isotropic_mesh"=isotropic_mesh,
                "isotropic_spde"=isotropic_spde,
                "anisotropic_mesh"=anisotropic_mesh,
                "anisotropic_spde"=anisotropic_spde,
                "Tri_Area"=Tri_Area,
                "TV"=TV, "E0"=E0, "E1"=E1, "E2"=E2,
                "anisotropic_mesh_triangles_over_land"=anisotropic_mesh_triangles_over_land,
                "barrier_list"=barrier_list )
  return(Return)
}
