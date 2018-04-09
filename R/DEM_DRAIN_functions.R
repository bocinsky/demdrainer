# The master function of the DEM draining.
# Given a digital elevation model (DEM), return
# a version with estimates for elevations under
# reservoirs and dams.
# Save downloaded data along the way.
drainDEM <- function(dem, raw.dir="./Input/NHD", out.dir="Output", dsn.vectors="Output/vectors", force.redo=FALSE) {
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(dsn.vectors, showWarnings = FALSE, recursive = TRUE)
  suppressWarnings(writeOGR(SPDFfromPolygon(polygonFromExtent(extent(dem), projection(dem))), dsn.vectors, "sim_poly", "ESRI Shapefile", overwrite_layer = TRUE))

  res.elevs.rast <- dem

  getNHD(dem)

  getNRCS_Dams(dem)

  Streams <- readOGR(dsn.vectors, "Streams")
  Reservoirs <- readOGR(dsn.vectors, "Reservoirs")
  Dams <- readOGR(dsn.vectors, "Dams")

  Dams <- gBuffer(Dams, width = 1 * res(dem)[1])

  reservoirsWithDams <- gWithinDistance(Reservoirs, Dams, dist = 2 * res(dem)[1], byid = T)

  reservoirsWithDams.vector <- vector("logical", length(Reservoirs))
  for (i in 1:length(Reservoirs)) {
    reservoirsWithDams.vector[i] <- any(reservoirsWithDams[, i])
  }

  Reservoirs <- Reservoirs[reservoirsWithDams.vector]

  damsWithReservoirs <- gWithinDistance(Dams, Reservoirs, dist = 100, byid = T)

  damsWithReservoirs.vector <- vector("logical", length(Dams))
  for (i in 1:length(Dams)) {
    damsWithReservoirs.vector[i] <- any(damsWithReservoirs[, i])
  }

  Dams <- Dams[damsWithReservoirs.vector, ]

  Streams.named <- Streams[!is.na(Streams$GNIS_Nm), ]
  Streams.unnamed <- Streams[is.na(Streams$GNIS_Nm), ]
  Streams.named.merged <- gLineMerge(Streams.named, byid = T, id = Streams.named$GNIS_Nm)
  Streams.unnamed.merged <- gLineMerge(Streams.unnamed)
  Streams <- spRbind(Streams.named.merged, Streams.unnamed.merged)
  Streams <- mergeStreams(Streams, projection(dem))

  # Join the dam and water polygons, and rasterize
  Dams <- spTransform(Dams, Reservoirs@proj4string)
  bad.data.vect <- gUnion(Dams, Reservoirs)
  bad.data.rast <- extract(dem, bad.data.vect, cellnumbers = T)

  res.elevs.rast[bad.data.rast[[1]][, 1]] <- NA


  # Find which streams have missing elevation data
  bad.data.vect <- spTransform(bad.data.vect, CRS(projection(dem)))
  Streams.gapped <- Streams[c(Streams %over% gBuffer(bad.data.vect, width = 500)) == 1]

  Streams.gapped <- forceMergeStreams(Streams.gapped, CRS(projection(dem)))

  res.elevs.rast <- bootstrapDrainDEM(res.elevs.rast, Streams.gapped, Reservoirs, dem)

  dem.final <- fillIDW(res.elevs.rast)

  # A 5x5 mean smooth
  dem.final.smoothed <- focal(dem.final, w = matrix(1, nrow = 5, ncol = 5), fun = mean, na.rm = T, pad = T)

  return(dem.final.smoothed)
}

# A method for removing modern features from the National Hydrography Dataset.
removeModernFeatures <- function(shapes, layers) {
  cat("Removing modern or unwanted features\n")
  # Only keep flowlines that represent Streams/Rivers (46006/46003), or
  # "Artificial Path," which are generally historic river courses under modern reservoirs
  if (!is.na(match("NHDFlowline", layers))) {
    index <- match("NHDFlowline", layers)
    shapes[[index]] <- shapes[[index]][shapes[[index]]$FCode %in% c(46003, 46006, 55800), ]
  }

  # Only keep waterbodies that are Lakes/ponds (39004/39001/39009), and reservoirs (39010)
  if (!is.na(match("NHDWaterbody", layers))) {
    index <- match("NHDWaterbody", layers)
    shapes[[index]] <- shapes[[index]][shapes[[index]]$FCode %in% c(39004, 39001, 39009, 39010), ]
  }
  return(shapes)
}

# A method for splitting features in the National Hydrography Dataset by their feature type.
splitAndExport <- function(shapes, layers, dsn.vectors="Output/vectors") {
  # Split shapefiles by feature types.
  cat("\nJoining shapefiles by feature types.\n")
  index <- match("NHDFlowline", layers)
  Streams <- shapes[[index]][shapes[[index]]$FCode %in% c(55800, 46006, 46003), ]

  index <- match("NHDWaterbody", layers)
  Lakes <- shapes[[index]][shapes[[index]]$FCode %in% c(39004, 39009, 39010, 39001), ]
  Lakes$ReachCd <- NULL
  # "Areas" are large Bureau of Reclamation reservoirs
  index <- match("NHDArea", layers)
  Areas <- shapes[[index]]
  # Join waterbodies and reservoir layers
  Reservoirs <- spRbind(Lakes, Areas)

  # Remove holes, and small lakes (< 400x400 meters)
  Reservoirs <- Reservoirs[Reservoirs$AreSqKm > 0.16, ]
  for (i in 1:length(Reservoirs)) {
    Reservoirs@polygons[[i]] <- remove.holes(Reservoirs@polygons[[i]])
  }

  # Remove water bodies that touch no streams
  Reservoirs <- spTransform(Reservoirs, CRS(projection(Streams)))
  #   Reservoirs <- Reservoirs[!is.na(c(Reservoirs %over% Streams)),]
  Reservoirs <- Reservoirs[!is.na(c(as(Reservoirs, "SpatialPolygons") %over% as(Streams, "SpatialLines"))), ]

  # Export final vector datasets for the study area.
  cat("\nExporting final vector and raster data.\n")
  suppressWarnings(writeOGR(Reservoirs, dsn.vectors, "Reservoirs", "ESRI Shapefile", overwrite_layer = TRUE))
  suppressWarnings(writeOGR(Streams, dsn.vectors, "Streams", "ESRI Shapefile", overwrite_layer = TRUE))
}

# A wrapper for the RGEOS function "gLineMerge" specifically to merge
# Stream segments of the National Hydrography Dataset
mergeStreams <- function(x, proj4string) {
  Streams.final <- vector("list")
  for (i in 1:length(x)) {
    temp.merge <- gLineMerge(SpatialLines(list(x[i, ]@lines[[1]]), proj4string = master.proj))
    for (j in 1:length(temp.merge@lines)) {
      for (k in 1:length(temp.merge@lines[[j]]@Lines)) {
        Streams.final <- c(Streams.final, temp.merge@lines[[j]]@Lines[[k]])
      }
    }
  }
  Streams.lines.list <- lapply(Streams.final, function(...) {
    Lines(..., ID = "lines")
  })
  for (i in 1:length(Streams.lines.list)) {
    Streams.lines.list[[i]]@ID <- paste("Stream", i)
  }
  Streams <- SpatialLines(Streams.lines.list, proj4string = CRS(proj4string))
  return(Streams)
}

# Strong-arm the joining of stream segments that remain unmerged by gLineMerge.
forceMergeStreams <- function(streams, proj4string) {
  # Do one more loop of gLineMerge to allow RGEOS to attempt to merge more segments
  done <- FALSE
  while (!done) {
    start.length <- length(streams)
    streams <- gLineMerge(streams, byid = TRUE)
    Streams.lines.list <- lapply(streams@lines[[1]]@Lines, function(...) {
      Lines(..., ID = "lines")
    })
    for (i in 1:length(Streams.lines.list)) {
      Streams.lines.list[[i]]@ID <- paste("Stream", i)
    }
    streams <- SpatialLines(Streams.lines.list, proj4string = proj4string)
    if (length(streams) == start.length) {
      done <- TRUE
    }
  }

  # RGEOS doesn't seem to do a complete job for some reason,
  # so this loop strong-arms the joining of stream segments that remain.
  done <- FALSE
  while (!done) {
    start.length <- length(streams)
    gapped.list <- vector("list")
    joined <- vector("logical", length(streams))
    counter <- 1
    for (i in 1:length(streams)) {
      for (j in 1:length(streams)) {
        if (joined[i]) break
        if (joined[j]) next
        if (i == j) next
        if (all(tail(streams@lines[[i]]@Lines[[1]]@coords, n = 1) == head(streams@lines[[j]]@Lines[[1]]@coords, n = 1))) {
          gapped.list <- c(gapped.list, Lines(Line(rbind(streams@lines[[i]]@Lines[[1]]@coords, streams@lines[[j]]@Lines[[1]]@coords)), ID = as.character(counter)))
          counter <- counter + 1
          joined[c(i, j)] <- TRUE
          break
        }
      }
    }

    for (i in 1:length(streams)) {
      if (!joined[i]) gapped.list <- c(gapped.list, streams@lines[[i]])
    }

    for (i in 1:length(gapped.list)) {
      gapped.list[[i]]@ID <- paste("Stream ", i, sep = "")
    }

    streams <- SpatialLines(gapped.list, proj4string = proj4string)

    if (length(streams) == start.length) done <- TRUE
  }

  return(streams)
}

# This loop is the work-horse of the algorithm. It first calculates the
# elevations of cells along each stream. If a stream BOTH enters and exits
# a reservoir, then it linearly interpolates between the entry and exit
# elevations for all elevations in between. It then assigns those elevations
# to the "holey" DEM along the stream, and also at a 1-cell buffer (this is dependent
# on the resolution of one's DEM).
#
# The algorithm then loops back and extracts the elevations along streams again.
# This time, streams that end at one of the streams that get interpolated in the
# first loop will have "ending" elevations, and may be interpolated themselves, and
# interpolated values added to the DEM.
#
# This process is repeated until no more streams may be interpolated.
bootstrapDrainDEM <- function(gappedDEM, streams, reservoirs, orig.dem) {
  done <- FALSE
  last <- FALSE
  extracts <- vector("list", length(streams))
  finished.lines <- vector("logical", length(extracts))
  counter <- 1
  while (!done | !last) {
    if (done) {
      last <- TRUE
    }
    cat("\nBootstrap-draining DEM, iteration", counter)
    counter <- counter + 1
    # A variable used to assess our progress
    start.length <- length(finished.lines[finished.lines == TRUE])

    # To make things faster, first calculate the elevations at the line-ends.
    # Only compute lines that have both ends.
    ends <- vector("list", length(streams))
    end.logic <- vector("logical", length(streams))
    for (i in 1:length(streams)) {
      if (finished.lines[i]) next
      ends[[i]] <- extract(gappedDEM, matrix(c(head(streams[i]@lines[[1]]@Lines[[1]]@coords, n = 1), tail(streams[i]@lines[[1]]@Lines[[1]]@coords, n = 1)), nrow = 2, byrow = T), cellnumbers = T)

      if (is.na(ends[[i]][1, 2])) {
        minLocalValue <- suppressWarnings(min(gappedDEM[adjacent(gappedDEM, ends[[i]][1, 1], directions = 8, include = T)[, 2]], na.rm = T))
        if (is.finite(minLocalValue)) {
          gappedDEM[ends[[i]][1, 1]] <- minLocalValue
          ends[[i]][1, 2] <- minLocalValue
        }
      }

      if (is.na(ends[[i]][2, 2])) {
        minLocalValue <- suppressWarnings(min(gappedDEM[adjacent(gappedDEM, ends[[i]][2, 1], directions = 8, include = T)[, 2]], na.rm = T))
        if (is.finite(minLocalValue)) {
          gappedDEM[ends[[i]][2, 1]] <- minLocalValue
          ends[[i]][2, 2] <- minLocalValue
        }
      }
      end.logic[i] <- all(!is.na(ends[[i]]))
    }

    # Calculare the elevations of the cells along each stream.
    # Don't re-calculate if a stream has already been interpolated.
    for (i in 1:length(streams)) {
      if (finished.lines[i] | !end.logic[i]) next
      extracts[[i]] <- extract(gappedDEM, streams[i], along = T, cellnumbers = T)[[1]]
    }

    if (last) {
      ends <- vector("list", length(streams))
      for (i in 1:length(streams)) {
        if (finished.lines[i]) next

        ends[[i]] <- extract(gappedDEM, matrix(c(head(streams[i]@lines[[1]]@Lines[[1]]@coords, n = 1), tail(streams[i]@lines[[1]]@Lines[[1]]@coords, n = 1)), nrow = 2, byrow = T), cellnumbers = TRUE)

        if (length((ends[[i]][, 1])[is.na(ends[[i]][, 2])]) > 0) {
          minimum <- getNearestMinimum(gappedDEM, (ends[[i]][, 1])[is.na(ends[[i]][, 2])])
          gappedDEM[(ends[[i]][, 1])[is.na(ends[[i]][, 2])]] <- minimum
        }
        extracts[[i]] <- extract(gappedDEM, streams[i], along = T, cellnumbers = T)[[1]]
      }
    }

    # For each extracted stream...
    for (i in 1:length(extracts)) {
      # Don't re-process if already finished
      if (finished.lines[i] | is.null(extracts[[i]])) {
        next
      }

      # Don't process if there are less than 4 total cells covered, or
      # if there are no NA cells. Consider these streams "finished"

      if (length(extracts[[i]][, 2]) < 4 | all(!is.na(extracts[[i]][, 2]))) {
        finished.lines[i] <- TRUE
        next
      }


      temp <- data.frame(1:length(extracts[[i]][, 2]), extracts[[i]][, 2], extracts[[i]][, 1])
      names(temp) <- c("index", "elevation", "cell")
      if (all(is.na(head(temp$elevation, n = 6))) | all(is.na(tail(temp$elevation, n = 6)))) next

      temp$elevation <- make.monotonic(temp$elevation)

      gaps <- vector("list")
      temp$group <- 0
      group.counter <- 0
      for (j in 2:nrow(temp)) {
        # Enters a batch of NAs
        if (!is.na(temp$elevation[j]) & is.na(temp$elevation[j + 1])) {
          if (!is.na(temp$elevation[j - 1]) & j != nrow(temp)) {
            group.counter <- group.counter + 1
          }
          temp$group[j] <- group.counter
        } else if (!is.na(temp$elevation[j]) & is.na(temp$elevation[j - 1])) {
          temp$group[j] <- group.counter
          if (!is.na(temp$elevation[j + 1]) & !is.na(temp$elevation[j + 2])) {
            group.counter <- group.counter + 1
          }
        } else {
          temp$group[j] <- group.counter
        }
      }

      temp.na <- temp[is.na(temp$elevation), ]

      for (j in unique(temp$group)) {
        if (any(is.na(temp[temp$group == j, ]$elevation)) & length(temp[temp$group == j, ]$elevation[!is.na(temp[temp$group == j, ]$elevation)]) >= 2) {
          temp[temp$group == j, ]$elevation <- approx(temp[temp$group == j, ]$index, temp[temp$group == j, ]$elevation, xout = temp[temp$group == j, ]$index)$y
        }
      }

      # Update the raster
      gappedDEM[temp$cell] <- temp$elevation

      local <- adjacent(gappedDEM, temp.na$cell, pairs = F, directions = 16)

      gappedDEM[local] <- gappedDEM[cellFromXY(gappedDEM, t(apply(xyFromCell(gappedDEM, local), 1, function(x) {
        nearestPointOnLine(streams[i]@lines[[1]]@Lines[[1]]@coords, x)
      })))]

      gappedDEM <- setMinMax(gappedDEM)

      finished.lines[i] <- TRUE
    }
    if (length(finished.lines[finished.lines == TRUE]) == start.length) done <- TRUE
  }

  if (any(!finished.lines)) {
    # Lakes where no interpolated values could be determined (usually
    # shallow stock-ponds) are re-filled with their original elevations.
    reservoirs <- spTransform(reservoirs, CRS(projection(gappedDEM)))
    unfinished.streams <- spTransform(Streams.gapped[!finished.lines], CRS(projection(gappedDEM)))
    if (any(is.na(c(reservoirs %over% unfinished.streams)))) {
      reservoirs.no.fill <- reservoirs[is.na(c(reservoirs %over% unfinished.streams)), ]
      reservoirs.no.fill.rast <- extract(orig.dem, reservoirs.no.fill, cellnumbers = T)
      for (i in 1:length(reservoirs.no.fill.rast)) {
        gappedDEM[reservoirs.no.fill.rast[[i]][, 1]] <- reservoirs.no.fill.rast[[i]][, 2]
      }
    }
  }

  return(gappedDEM)
}

# A function to get the nearest minimum value from a raster, given a cell location.
getNearestMinimum <- function(raster, cell) {
  done <- FALSE
  radius <- 1
  while (!done) {
    nRowCols <- (radius * 2) + 1
    neighbors <- matrix(nrow = nRowCols, ncol = nRowCols, data = 1)
    neighbors[radius + 1, radius + 1] <- 0
    minimum <- suppressWarnings(min(raster[adjacent(raster, cell, directions = neighbors)[, 2]], na.rm = T))
    radius <- radius + 1
    if (is.finite(minimum)) done <- TRUE
  }
  return(minimum)
}

# A function that fills missing values in a raster using Inverse Distance Weighted interpolation
# with an IDW exponent of 0.5.
fillIDW <- function(gappedDEM) {
  # Finally, interpolate the missing values using IDW interpolation!
  # This uses the known elevations to derive an Inverse Distance Weighted model
  # of spatial prediction, then predicts elevation values at the still-missing points.
  gappedDEM.df <- data.frame(coordinates(gappedDEM), values(gappedDEM))
  names(gappedDEM.df) <- c("x", "y", "ELEVATION")
  gappedDEM.df.known <- gappedDEM.df[!is.na(gappedDEM.df$ELEVATION), ]
  gappedDEM.df.unknown <- gappedDEM.df[is.na(gappedDEM.df$ELEVATION), ]
  gappedDEM.idw.model <- gstat(id = "ELEVATION", formula = ELEVATION ~ 1, locations = ~ x + y, data = gappedDEM.df.known, nmax = 7, set = list(idp = 0.5))
  gappedDEM.df.new <- predict.gstat(gappedDEM.idw.model, newdata = gappedDEM.df.unknown)
  gappedDEM.df.new$CELL <- cellFromXY(gappedDEM, as.matrix(gappedDEM.df.new[, 1:2]))
  gappedDEM.final <- gappedDEM
  gappedDEM.final[gappedDEM.df.new$CELL] <- gappedDEM.df.new$ELEVATION.pred
  return(gappedDEM.final)
}

# A function that loads the National Hydrography Dataset for a provided study area defined by "x."
getNHD <- function(x, raw.dir="./Input/NHD", remove.modern=TRUE, out.dir="Output", dsn.vectors="Output/vectors", force.redo=FALSE) {
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(dsn.vectors, showWarnings = FALSE, recursive = TRUE)
  suppressWarnings(writeOGR(SPDFfromPolygon(polygonFromExtent(extent(x), projection(x))), dsn.vectors, "sim_poly", "ESRI Shapefile", overwrite_layer = TRUE))

  # Set the layers to be extracted
  layers <- c("NHDFlowline", "NHDWaterbody", "NHDArea")

  HUC4 <- loadHUC4(x, raw.dir = raw.dir, out.dir = out.dir, dsn.vectors = dsn.vectors, force.redo = force.redo)

  area.list <- getAreaList(HUC4)

  shapes <- loadNHDSubregions(x, layers = layers, area.list = area.list, raw.dir = raw.dir, out.dir = out.dir, dsn.vectors = dsn.vectors, force.redo = force.redo)

  splitAndExport(shapes = shapes, layers = layers, dsn.vectors = dsn.vectors)

  #   system(paste('ogr2ogr -overwrite -f "FileGDB" ',vep.output.dir,'/NHD_vectors.gdb ',dsn, sep=''))
  #   unlink(paste(vep.output.dir,"/NHD_vectors",sep=''),recursive = T, force = T)

  cat("\nNHD processed for study area defined by x")
}

# A function that loads the Natural Resources Conservation Service dam polygons
getNRCS_Dams <- function(x, aggFactor=1, raw.dir="./Input/NRCS", layers=NULL, out.dir="Output", dsn.vectors="Output/vectors", force.redo=FALSE) {
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(dsn.vectors, showWarnings = FALSE, recursive = TRUE)
  suppressWarnings(writeOGR(SPDFfromPolygon(polygonFromExtent(extent(x), projection(x))), dsn.vectors, "sim_poly", "ESRI Shapefile", overwrite_layer = TRUE))


  # Load the NRCS study areas
  cat("\nLoading map of NRCS survey areas\n")
  NRCS.areas <- loadNRCSStudyAreas(x, raw.dir = raw.dir, dsn.vectors = dsn.vectors, force.redo = force.redo)

  # Load the NRCS mapunit polygons
  cat("\nLoading NRCS mapunit polygons\n")
  NRCS.polys <- loadNRCSMapUnitPolygons(x = x, raw.dir = raw.dir, dsn.vectors = dsn.vectors, force.redo = force.redo)

  # Load tabular soil data for study area
  getSoilData(x, raw.dir = raw.dir, dsn.vectors = dsn.vectors, out.dir = out.dir, force.redo = force.redo)

  # Create a final NRCS soils polygon file
  NRCS.final <- NRCS.polys[, c(4)]

  soils <- suppressWarnings(readOGR(paste(dsn.vectors, sep = ""), layer = "soils"))
  mapunits <- read.csv(paste(out.dir, "/mapunit.csv", sep = ""))
  dam.mukeys <- mapunits[mapunits$muname == "Dam", ]$mukey
  dams <- soils[soils$MUKEY %in% dam.mukeys, ]
  dams <- spTransform(dams, CRS(projection(x)))

  writeOGR(dams, dsn.vectors, "Dams", "ESRI Shapefile", overwrite_layer = TRUE)
}

# A function that loads a map of NRCS study areas (the individual soil surveys)
loadNRCSStudyAreas <- function(x, raw.dir="./Input/NRCS", dsn.vectors="Output/vectors", force.redo=F) {
  # Get NRCS study areas, and save them to disk
  if (!("NRCS_SurveyAreas" %in% ogrListLayers(dsn.vectors)) | force.redo) {
    NRCS.areas <- getNRCSStudyAreas(x, raw.dir = raw.dir)
    suppressWarnings(writeOGR(NRCS.areas, dsn.vectors, "NRCS_SurveyAreas", "ESRI Shapefile", overwrite_layer = TRUE))
  } else {
    NRCS.areas <- readOGR(dsn.vectors, "NRCS_SurveyAreas")
  }

  # Check to see if all survey areas are available
  if (0 %in% NRCS.areas@data$iscomplete) {
    cat("WARNING! Some of the soil surveys in your area are unavailable.\n")
    cat("Soils and productivity data will have holes.\n")
    cat("Missing areas:\n")
    cat(as.vector(NRCS.areas@data[NRCS.areas@data$iscomplete == 0, ]$areasymbol))
    cat("\n\n")
    cat("Continuing with processing available soils.\n\n")
    NRCS.areas <- NRCS.areas[NRCS.areas@data$iscomplete != 0, ]
  }

  return(NRCS.areas)
}

# A function that downloads the maps of NRCS study areas, and crops it to a given area.
getNRCSStudyAreas <- function(x, raw.dir="./Input/NRCS") {
  # Import the shapefile of NRCS study areas.
  # This is available at
  # http://soildatamart.sc.egov.usda.gov/download/StatusMaps/soilsa_a_nrcs.zip
  if (url.exists("http://websoilsurvey.sc.egov.usda.gov/DataAvailability/SoilDataAvailabilityShapefile.zip")) {
    f <- CFILE(paste(raw.dir, "/SoilDataAvailabilityShapefile.zip", sep = ""), mode = "wb")
    curlPerform(url = "http://websoilsurvey.sc.egov.usda.gov/DataAvailability/SoilDataAvailabilityShapefile.zip", writedata = f@ref)
    close(f)
  } else {
    stop("Unable to download NRCS study area shapefile! Please find a copy of the NRCS study areas shapefile, available from the USGS.")
  }
  unzip(paste(raw.dir, "/SoilDataAvailabilityShapefile.zip", sep = ""), exdir = paste(raw.dir, "/SoilDataAvailabilityShapefile", sep = ""))
  NRCS.areas <- suppressWarnings(readOGR(paste(raw.dir, "/SoilDataAvailabilityShapefile", sep = ""), layer = "soilsa_a_nrcs", verbose = FALSE))

  # Get a list of NHD subregions within the project study area
  extent.poly <- polygonFromExtent(extent(x), projection(x))
  extent.poly.z15 <- spTransform(extent.poly, CRS(projection(NRCS.areas)))
  NRCS.areas <- crop.to.studyArea(NRCS.areas, extent.poly.z15)
  NRCS.areas <- spTransform(NRCS.areas, CRS(projection(x)))

  return(NRCS.areas)
}

# A function that loads NRCS Map Unit polygons given a study area.
loadNRCSMapUnitPolygons <- function(x, raw.dir="./Input/NRCS", dsn.vectors="Output/vectors", force.redo=F) {
  if (!("soils" %in% ogrListLayers(dsn.vectors)) | force.redo) {
    NRCS.areas <- loadNRCSStudyAreas(x, raw.dir = raw.dir, dsn.vectors = dsn.vectors, force.redo = F)
    NRCS.polys <- getNRCSMapUnitPolygons(x, NRCS.areas, raw.dir = raw.dir)
    # Export final vector dataset for the study area.
    cat("Exporting final vector dataset for the study area.\n")
    writeOGR(NRCS.polys, dsn.vectors, "soils", "ESRI Shapefile", overwrite_layer = TRUE)
  } else {
    NRCS.polys <- readOGR(dsn.vectors, "soils")
  }
  return(NRCS.polys)
}

# A function that downloads NRCS map unit polygons, and transforms and crops them given a study area.
getNRCSMapUnitPolygons <- function(x, NRCS.areas, raw.dir="./Input/NRCS") {
  # Load raw NRCS Map Unit polygons from the regions specified.
  cat("Loading the NRCS soil survey data for each region.\n")
  NRCS.polys <- vector("list", length(NRCS.areas))
  for (i in 1:length(NRCS.areas)) {
    cat("\n(Down)Loading soils data for", as.character(NRCS.areas@data$areasymbol[i]), "\n")
    if (!file.exists(paste(raw.dir, "/wss_SSA_", NRCS.areas@data$areasymbol[i], "_[", as.Date(NRCS.areas@data$saverest[i]), "].zip", sep = ""))) {
      f <- CFILE(paste(raw.dir, "/wss_SSA_", NRCS.areas@data$areasymbol[i], "_[", as.Date(NRCS.areas@data$saverest[i]), "].zip", sep = ""), mode = "wb")
      curlPerform(url = paste("http://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/wss_SSA_", NRCS.areas@data$areasymbol[i], "_[", as.Date(NRCS.areas@data$saverest[i]), "].zip", sep = ""), writedata = f@ref)
      close(f)
    }

    unzip(paste(raw.dir, "/wss_SSA_", NRCS.areas@data$areasymbol[i], "_[", as.Date(NRCS.areas@data$saverest[i]), "].zip", sep = ""), exdir = raw.dir)

    NRCS.polys[[i]] <- readOGR(paste(raw.dir, "/wss_SSA_", NRCS.areas@data$areasymbol[i], "_[", as.Date(NRCS.areas@data$saverest[i]), "]/spatial", sep = ""), layer = paste("soilmu_a_", tolower(NRCS.areas@data$areasymbol[i]), sep = ""))
    # Change all spatial IDs to prepare from merging
    NRCS.polys[[i]] <- spChFIDs(NRCS.polys[[i]], as.character(paste(NRCS.areas@data$areasymbol[i], "_", row.names(NRCS.polys[[i]]@data), sep = "")))
  }

  # Merging all NRCS Map Unit polygons
  NRCS.polys <- do.call("rbind", NRCS.polys)

  # Transform NRCS Map Unit polygons to CRS of x
  cat("Transforming to the CRS of x\n")
  NRCS.polys <- spTransform(NRCS.polys, CRS(projection(x)))

  # Crop to area of x
  cat("Cropping NRCS Map Unit polygons to the extent of x\n")
  extent.poly <- polygonFromExtent(extent(x), projection(x))
  NRCS.polys <- crop.to.studyArea(NRCS.polys, extent.poly)

  # Make MUKEY column numeric, and not factor
  NRCS.polys$MUKEY <- as.numeric(as.character(NRCS.polys$MUKEY))

  return(NRCS.polys)
}

# A function that extracts a particular table ("x") from the NRCS soils database, for a
# list of study areas, and aggregates them together.
loadAndAggregateSoilTable <- function(x, NRCS.areas, raw.dir="./Input/NRCS") {
  tables <- vector("list", length(NRCS.areas))
  for (i in 1:length(NRCS.areas)) {
    if (length(readLines(paste(raw.dir, "/wss_SSA_", NRCS.areas@data$areasymbol[i], "_[", as.Date(NRCS.areas@data$saverest[i]), "]/tabular/", x, ".txt", sep = ""))) > 0) {
      tables[[i]] <- read.delim(paste(raw.dir, "/wss_SSA_", NRCS.areas@data$areasymbol[i], "_[", as.Date(NRCS.areas@data$saverest[i]), "]/tabular/", x, ".txt", sep = ""), header = F, sep = "|")
    }
  }
  table <- do.call("rbind", tables)
}

# A function that extracts the "mapunit" table from a set of sudy areas, and aggregates it.
getSoilData <- function(x, raw.dir="./Input/NRCS", dsn.vectors="Output/vectors", out.dir="Output", force.redo=F) {
  # Load the survey area and soils data
  NRCS.areas <- loadNRCSStudyAreas(x = x, raw.dir = raw.dir, dsn.vectors = dsn.vectors, force.redo = F)
  NRCS.polys <- loadNRCSMapUnitPolygons(x = x, raw.dir = raw.dir, dsn.vectors = dsn.vectors, force.redo = F)

  # Load the primary mapunit data
  if (!file.exists(paste(out.dir, "/mapunit.csv", sep = "")) | force.redo) {
    NRCS.mapunit <- loadAndAggregateSoilTable("mapunit", NRCS.areas, raw.dir = raw.dir)[, c(1:5, 12:24)]
    names(NRCS.mapunit) <- c("musym", "muname", "mukind", "mustatus", "muacres", "farmlndcl", "muhelcl", "muwathelcl", "muwndhelcl", "interpfocus", "invesintens", "iacornsr", "nhiforsoigrp", "nhspiagr", "vtsepticsyscl", "mucertstat", "lkey", "mukey")
    NRCS.mapunit <- NRCS.mapunit[NRCS.mapunit$mukey %in% unique(NRCS.polys$MUKEY), ]
    write.csv(NRCS.mapunit, paste(out.dir, "/mapunit.csv", sep = ""), row.names = F)
  }
}


# A method to crop vector shapefiles to a given study area.
crop.to.studyArea <- function(x, y) {
  gI <- gIntersects(x, y, byid = TRUE)
  out <- vector(mode = "list", length = length(which(gI)))
  ii <- 1
  if (length(out) == 0) return(NULL)
  for (i in seq(along = gI)) {
    if (gI[i]) {
      out[[ii]] <- gIntersection(x[i, ], y)
      row.names(out[[ii]]) <- row.names(x)[i]
      ii <- ii + 1
    }
  }
  out_class <- sapply(out, class)

  if (class(x)[1] == "SpatialLinesDataFrame") {
    ri <- do.call("rbind", out[out_class == "SpatialLines"])
    ri <- SpatialLinesDataFrame(ri, as.data.frame(x))
    return(ri)
  } else if (class(x)[1] == "SpatialPolygonsDataFrame") {
    ri <- do.call("rbind", out[out_class == "SpatialPolygons"])
    ri <- SpatialPolygonsDataFrame(ri, as.data.frame(x)[row.names(ri), ])
    return(ri)
  } else if (class(x)[1] == "SpatialPointsDataFrame") {
    ri <- do.call("rbind", out[out_class == "SpatialPoints"])
    ri <- SpatialPointsDataFrame(ri, as.data.frame(x)[row.names(ri), ], match.ID = TRUE)
    return(ri)
  }
}

# A method to create a spatial polygon from the extent of a spatial* or raster* object
polygonFromExtent <- function(x, proj4string) {
  extent.matrix <- rbind(c(x@xmin, x@ymin), c(x@xmin, x@ymax), c(x@xmax, x@ymax), c(x@xmax, x@ymin), c(x@xmin, x@ymin)) # clockwise, 5 points to close it
  extent.SP <- SpatialPolygons(list(Polygons(list(Polygon(extent.matrix)), "extent")), proj4string = CRS(proj4string))
  return(extent.SP)
}

# A method to create a SpatialPolygonsDataFrame from a single polygon.
SPDFfromPolygon <- function(x) {
  IDs <- sapply(slot(x, "polygons"), function(x) slot(x, "ID"))
  df <- data.frame(rep(0, length(IDs)), row.names = IDs)
  x <- SpatialPolygonsDataFrame(x, df)
  return(x)
}

# A function that takes a sometimes erratic stream profile
# and makes it monotonic (strictly downward or upward sloping)
# for the purposes of interpolating across gaps.
# The monotonic values do not replace the existing values in the DEM.
make.monotonic <- function(x) {
  for (i in 2:length(x)) {
    if (!is.na(x[i]) & !is.na(x[i - 1])) {
      if (!is.na(x[i - 1]) & x[i] > x[i - 1]) {
        x[i] <- x[i - 1]
      }
    }
  }
  return(x)
}

# A method to plot rasters as a PNG
raster.png <- function(x, zlim, palette, file.name, title, legend.lab, extra.vector.plots.function, extra.label.plots.function) {
  image.width <- 7
  image.height <- 7
  x <- setMinMax(x)

  if (is.na(zlim[1])) {
    zlim <- c(min(getValues(x)), max(getValues(x)))
  }

  ratio.raster <- ncol(x) / nrow(x)

  if (ratio.raster >= 1) {
    image.height <- max(image.width * (1 / ratio.raster), 4)
  } else {
    image.width <- image.height * ratio.raster
  }

  colors <- palette
  quartz(file = file.name, width = image.width, height = image.height, antialias = FALSE, bg = "white", type = "png", family = "Gulim", pointsize = 1, dpi = 1200)
  if (ratio.raster >= 1) {
    plot.width <- 6
    plot.height <- plot.width * (1 / ratio.raster)
    inch <- (extent(x)@xmax - extent(x)@xmin) / (image.width - 1)
    par(mai = c((image.height - plot.height) / 2, 0.5, (image.height - plot.height) / 2, 0.5), omi = c(0, 0, 0, 0), pin = c(plot.width, plot.height))
  } else {
    plot.height <- 6
    plot.width <- plot.height * ratio.raster
    inch <- (extent(x)@ymax - extent(x)@ymin) / (image.height - 1)
    par(mai = c(0.5, (image.width - plot.width) / 2, 0.5, (image.width - plot.width) / 2), omi = c(0, 0, 0, 0), pin = c(plot.width, plot.height))
  }
  par(bg = "white", fg = "black", col.lab = "black", col.main = "black", col.axis = "black", family = "Helvetica Bold", lend = 2, ljoin = 1)
  plot(x, maxpixels = ncell(x), xlim = c(extent(x)@xmin, extent(x)@xmax), ylim = c(extent(x)@ymin, extent(x)@ymax), zlim = zlim, xlab = "", ylab = "", axes = FALSE, main = "", asp = 1, col = colors, useRaster = FALSE, legend = FALSE)

  plot(sim.poly, add = TRUE)

  if (!is.na(extra.vector.plots.function)) {
    FUN <- match.fun(extra.vector.plots.function)
    FUN()
  }

  # plot(x,add=T,col=colors, useRaster=FALSE, legend=FALSE)
  xseq <- c(seq(extent(x)@xmin, extent(x)@xmax, by = round_any((extent(x)@xmax - extent(x)@xmin) / 10, 1000)), extent(x)@xmax)
  yseq <- c(seq(extent(x)@ymin, extent(x)@ymax, by = round_any((extent(x)@xmax - extent(x)@xmin) / 10, 1000)), extent(x)@ymax)

  if (abs(diff(tail(xseq, n = 2))) < ((extent(x)@xmax - extent(x)@xmin) / (length(xseq) - 1))) {
    xseq <- xseq[-(length(xseq) - 1)]
  }

  if (abs(diff(tail(yseq, n = 2))) < ((extent(x)@xmax - extent(x)@xmin) / (length(yseq) - 1))) {
    yseq <- yseq[-(length(yseq) - 1)]
  }

  #   par(tck=(-1*(0.01)), family="Helvetica", lend=2)
  #   axis(3, at=xseq, cex.axis=7, lwd=-1, lwd.ticks=1, padj=-.4, pos=(extent(x)@ymax), labels=FALSE)
  #   text(xseq, extent(x)@ymax+(0.08*inch), labels = c(xseq[-length(xseq)],''), srt = 45, pos = 4, xpd = TRUE, cex=7)
  #   axis(4, at=yseq, cex.axis=7, lwd=-1, lwd.ticks=1, padj=-.4, pos=(extent(x)@xmax), labels=FALSE)
  #   text(extent(x)@xmax+(0.08*inch), yseq, labels = c(yseq[-length(yseq)],''), srt = 45, pos = 4, xpd = TRUE, cex=7)
  #
  #   text(extent(x)@xmax, extent(x)@ymax+(0.08*inch), labels = extent(x)@xmax, srt = 90, pos = 4, xpd = TRUE, cex=7)
  #   text(extent(x)@xmax+(0.08*inch), extent(x)@ymax,labels = extent(x)@ymax, srt = 0, pos = 4, xpd = TRUE, cex=7)
  #
  plot(sim.poly, add = TRUE)

  par(bty = "o", tck = -(1 / 4))

  if (!is.na(extra.label.plots.function)) {
    FUN <- match.fun(extra.label.plots.function)
    FUN(x, inch)
  } else {
    plot(x, maxpixels = ncell(x), add = TRUE, legend.only = TRUE, col = colors, smallplot = c(0.5, (((image.width - plot.width) / 2) + plot.width) / image.width, (((image.height - plot.height) / 2) - 0.15) / image.height, (((image.height - plot.height) / 2) - 0.05) / image.height), horiz = TRUE, axis.args = list(cex.axis = 7, lwd = 0, lwd.ticks = 0.15 * 6, pos = -0.25, padj = 1.1, bty = "o"), zlim = zlim, legend.args = list(text = legend.lab, side = 1, font = 2, line = 16, cex = 8))
  }

  text(extent(x)@xmin, (extent(x)@ymin - (0.14 * inch)), labels = title, pos = 4, xpd = TRUE, cex = 14, font = 2)

  #   text(extent(x)@xmin,(extent(x)@ymin-(0.21*inch)), labels='R. Kyle Bocinsky', pos = 4, xpd = TRUE, cex=4,font=2)
  #   text(extent(x)@xmin,(extent(x)@ymin-(0.28*inch)), labels=format(Sys.Date(), "%d %B %Y"), pos = 4, xpd = TRUE, cex=4,font=2)
  #   text(extent(x)@xmin,(extent(x)@ymin-(0.35*inch)), labels=paste('PROJ.4 String: ',projection(x),sep=''), pos = 4, xpd = TRUE, cex=4,font=2)

  dev.off()
}

# A method to plot multiband (RBG) rasters as a PNG
rasterRGB.png <- function(x, file.name, title, extra.vector.plots.function, extra.label.plots.function) {
  image.width <- 7
  image.height <- 7
  x <- setMinMax(x)

  ratio.raster <- ncol(x) / nrow(x)

  if (ratio.raster >= 1) {
    image.height <- max(image.width * (1 / ratio.raster), 4)
  } else {
    image.width <- image.height * ratio.raster
  }

  colors <- palette
  quartz(file = file.name, width = image.width, height = image.height, antialias = FALSE, bg = "white", type = "png", family = "Gulim", pointsize = 1, dpi = 1200)
  if (ratio.raster >= 1) {
    plot.width <- 6
    plot.height <- plot.width * (1 / ratio.raster)
    inch <- (extent(x)@xmax - extent(x)@xmin) / (image.width - 1)
    par(mai = c((image.height - plot.height) / 2, 0.5, (image.height - plot.height) / 2, 0.5), omi = c(0, 0, 0, 0), pin = c(plot.width, plot.height))
  } else {
    plot.height <- 6
    plot.width <- plot.height * ratio.raster
    inch <- (extent(x)@ymax - extent(x)@ymin) / (image.height - 1)
    par(mai = c(0.5, (image.width - plot.width) / 2, 0.5, (image.width - plot.width) / 2), omi = c(0, 0, 0, 0), pin = c(plot.width, plot.height))
  }
  par(bg = "white", fg = "black", col.lab = "black", col.main = "black", col.axis = "black", family = "Helvetica Bold", lend = 2, ljoin = 1)
  sim.raster.new <- sim.raster
  sim.raster.new[] <- 0
  plot(sim.raster.new, xlab = "", ylab = "", axes = FALSE, main = "", asp = 1, col = "white", useRaster = FALSE, legend = FALSE)
  plotRGB(x, maxpixels = ncell(x), add = T)

  # plot(x,add=T,col=colors, useRaster=FALSE, legend=FALSE)
  xseq <- c(seq(extent(x)@xmin, extent(x)@xmax, by = round_any((extent(x)@xmax - extent(x)@xmin) / 10, 1000)), extent(x)@xmax)
  yseq <- c(seq(extent(x)@ymin, extent(x)@ymax, by = round_any((extent(x)@xmax - extent(x)@xmin) / 10, 1000)), extent(x)@ymax)

  if (abs(diff(tail(xseq, n = 2))) < round_any((extent(x)@xmax - extent(x)@xmin) / 15, 1000)) {
    xseq <- xseq[-(length(xseq) - 1)]
  }

  if (abs(diff(tail(yseq, n = 2))) < round_any((extent(x)@xmax - extent(x)@xmin) / 15, 1000)) {
    yseq <- yseq[-(length(yseq) - 1)]
  }

  #   par(tck=(-1*(0.01)), family="Helvetica", lend=2)
  #   axis(3, at=xseq, cex.axis=7, lwd=-1, lwd.ticks=1, padj=-.4, pos=(extent(x)@ymax), labels=FALSE)
  #   text(xseq, extent(x)@ymax+(0.08*inch), labels = c(xseq[-length(xseq)],''), srt = 45, pos = 4, xpd = TRUE, cex=7)
  #   axis(4, at=yseq, cex.axis=7, lwd=-1, lwd.ticks=1, padj=-.4, pos=(extent(x)@xmax), labels=FALSE)
  #   text(extent(x)@xmax+(0.08*inch), yseq, labels = c(yseq[-length(yseq)],''), srt = 45, pos = 4, xpd = TRUE, cex=7)
  #
  #   text(extent(x)@xmax, extent(x)@ymax+(0.08*inch), labels = extent(x)@xmax, srt = 90, pos = 4, xpd = TRUE, cex=7)
  #   text(extent(x)@xmax+(0.08*inch), extent(x)@ymax,labels = extent(x)@ymax, srt = 0, pos = 4, xpd = TRUE, cex=7)

  plot(sim.poly, add = TRUE)

  if (!is.na(extra.vector.plots.function)) {
    FUN <- match.fun(extra.vector.plots.function)
    FUN()
  }

  plot(sim.poly, add = TRUE)

  par(bty = "o", tck = -(1 / 4))

  text(extent(x)@xmin, (extent(x)@ymin - (0.14 * inch)), labels = title, pos = 4, xpd = TRUE, cex = 14, font = 2)

  #   text(extent(x)@xmin,(extent(x)@ymin-(0.21*inch)), labels='R. Kyle Bocinsky', pos = 4, xpd = TRUE, cex=4,font=2)
  #   text(extent(x)@xmin,(extent(x)@ymin-(0.28*inch)), labels=format(Sys.Date(), "%d %B %Y"), pos = 4, xpd = TRUE, cex=4,font=2)
  #   text(extent(x)@xmin,(extent(x)@ymin-(0.35*inch)), labels=paste('PROJ.4 String: ',projection(x),sep=''), pos = 4, xpd = TRUE, cex=4,font=2)

  dev.off()
}
