

latlon_to_SpatialPolygons = function( lat, lon, alpha=1 ){

  #==================================
  # Convert alpha shapes to polygons
  #==================================
  #alphashape_1 = ashape(xym, alpha=alpha)
  #ashape2poly = function(ashape){
  #  # Convert node numbers into characters
  #  ashape$edges[,1] = as.character(ashape$edges[,1])
  #  ashape_graph = graph_from_edgelist(ashape$edges[,1:2], directed = FALSE)
  #  if (!is.connected(ashape_graph)) {
  #    stop("Graph not connected")
  #  }
  #  if (any(degree(ashape_graph) != 2)) {
  #    stop("Graph not circular")
  #  }
  #  if (clusters(ashape_graph)$no > 1) {
  #    stop("Graph composed of more than one circle")
  #  }
  #  # Delete one edge to create a chain
  #  cut_graph = ashape_graph - E(ashape_graph)[1]
  #  # Find chain end points
  #  ends = names(which(degree(cut_graph) == 1))
  #  path = get.shortest.paths(cut_graph, ends[1], ends[2])$vpath[[1]]
  #  # this is an index into the points
  #  pathX = as.numeric(V(ashape_graph)[path]$name)
  #  # join the ends
  #  pathX = c(pathX, pathX[1])
  #  return(pathX)
  #}
  #alphapoly_1 = ashape2poly(alphashape_1)

  #========================
  # Convert arcs to lines
  #========================



  #========================
  # Convert arcs to lines
  #========================

  # function to convert an arc into line segments
  # given the center of the arc, the radius, the vector, and the angle (radians)
  arc2line = function(center, r, vector, theta, npoints = 100) {
    # Get the angles at the extremes of the arcs
    angles = alphahull::anglesArc(vector, theta)
    # Generate sequence of angles along the arc to determine the points
    seqang = seq(angles[1], angles[2], length = npoints)
    # Generate x coordinates for points along the arc
    x = center[1] + r * cos(seqang)
    # Generate y coordinates for points along the arc
    y = center[2] + r * sin(seqang)
    coords.xy = cbind(x,y)
    line = sp::Line(coords = coords.xy)
    return(line)
  }

  #==========================================
  #Convert a-hull into a SpatialLines object
  #==========================================
  ahull2lines = function(hull){
    arclist = hull$arcs
    lines = list()
    for (i in 1:nrow(arclist)) {
      # Extract the attributes of arc i
      center_i = arclist[i, 1:2]
      radius_i = arclist[i, 3]
      vector_i = arclist[i, 4:5]
      theta_i = arclist[i, 6]
      # Convert arc i into a Line object
      line_i = arc2line(center = center_i, r = radius_i, vector = vector_i, theta = theta_i)
      list_length = length(lines)
      if(list_length > 0){
        # If a line has already been added to the list of lines
        # Define last_line_coords as the coordinates of the last line added to the list before the ith line
        last_line_coords = lines[[list_length]]@coords
      }
      if(i == 1){
        # Add the first line to the list of lines
        lines[[i]] = line_i
      } else if(all.equal(line_i@coords[1,], last_line_coords[nrow(last_line_coords),])){
        # If the first coordinate in the ith line is equal to the last coordinate in the previous line
        # then those lines should be connected
        # Row bind the coordinates for the ith line to the coordinates of the previous line in the list
        lines[[list_length]]@coords = rbind(last_line_coords, line_i@coords[2:nrow(line_i@coords),])
      } else {
        # If the first coordinate in the ith line does not match the last coordinate in the previous line
        # then the ith line represents a new line
        # Add the ith line to the list as a new element
        lines[[length(lines) + 1]] = line_i
      }
    }

    # Convert the list of lines to a Line object
    lines = sp::Lines(lines, ID = 'l')
    # Convert the Line object to a SpatialLines object
    sp_lines = sp::SpatialLines(list(lines))
    return(sp_lines)
  }

  #===============================================
  # Convert spatial lines object to spatialpolygon
  #===============================================
  spLines2poly = function(sp_lines){
    # Extract the lines slot
    lines_slot = sp_lines@lines[[1]]
    # Create a list of booleans indicating whether a given Line represents a polygon
    poly_bool = sapply(lines_slot@Lines, function(x){
      coords = lines_slot@Lines[[1]]@coords
      # Check if the first coordinate in the line is the same as the last
      all.equal(coords[1,], coords[nrow(coords),])
    })
    # Pull out the lines that form polygons
    poly_lines = sp_lines[poly_bool]
    poly_lines_slot = poly_lines@lines
    # Create SpatialPolygons
    sp_polys = sp::SpatialPolygons(list(sp::Polygons(lapply(poly_lines_slot, function(x) {
      sp::Polygon(slot(slot(x, "Lines")[[1]], "coords"))
    }), ID = "1")))
    return(sp_polys)
  }

  #===============================================
  # Bundle pieces
  #===============================================
  ahull2poly = function(hull){
    # Convert the alpha hull to SpatialLines
    hull2SpatialLines = ahull2lines(hull)
    # Convert SpatialLines to SpatialPolygon
    SpatialLines2SpatialPolygon = spLines2poly(hull2SpatialLines)
    return(SpatialLines2SpatialPolygon)
  }

  lonlat = cbind( lon, lat )
  lonlat = unique(lonlat)
  alphahull = alphahull::ahull(lonlat, alpha=alpha )
  SpatialPolygons = ahull2poly(alphahull)
  return(SpatialPolygons)
}

#setwd("C:/Users/James.Thorson/Desktop/Work files/AFSC/2020-06 -- Question from Ellen")
#temp2 = read.csv("EBS_stations.csv") # Feb 2020 northern Bering Sea 6 jellyfish species
#lat = temp2$lat
#lon = temp2$long
#
#SpatialPolygons = latlon_to_SpatialPolygons( lat=lat, lon=lon, alpha=3 )
#sp::plot(SpatialPolygons)
