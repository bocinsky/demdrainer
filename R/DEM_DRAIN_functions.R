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
      ends[[i]] <- raster::extract(gappedDEM,
        matrix(c(head(streams[i, ]@lines[[1]]@Lines[[1]]@coords,
          n = 1
        ), tail(streams[i, ]@lines[[1]]@Lines[[1]]@coords,
          n = 1
        )),
        nrow = 2,
        byrow = T
        ),
        cellnumbers = T
      )

      if (is.na(ends[[i]][1, 2])) {
        minLocalValue <- suppressWarnings(
          min(
            gappedDEM[raster::adjacent(gappedDEM, ends[[i]][1, 1], directions = 8, include = T)[, 2]],
            na.rm = T
          )
        )

        if (is.finite(minLocalValue)) {
          gappedDEM[ends[[i]][1, 1]] <- minLocalValue
          ends[[i]][1, 2] <- minLocalValue
        }
      }

      if (is.na(ends[[i]][2, 2])) {
        minLocalValue <- suppressWarnings(
          min(
            gappedDEM[raster::adjacent(gappedDEM, ends[[i]][2, 1], directions = 8, include = T)[, 2]],
            na.rm = T
          )
        )
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
      extracts[[i]] <- raster::extract(gappedDEM, streams[i, ], along = T, cellnumbers = T)[[1]]
    }

    if (last) {
      ends <- vector("list", length(streams))
      for (i in 1:length(streams)) {
        if (finished.lines[i]) next

        ends[[i]] <- raster::extract(gappedDEM,
          matrix(c(
            head(streams[i, ]@lines[[1]]@Lines[[1]]@coords, n = 1),
            tail(streams[i, ]@lines[[1]]@Lines[[1]]@coords, n = 1)
          ),
          nrow = 2,
          byrow = T
          ),
          cellnumbers = TRUE
        )

        if (length((ends[[i]][, 1])[is.na(ends[[i]][, 2])]) > 0) {
          minimum <- getNearestMinimum(gappedDEM, (ends[[i]][, 1])[is.na(ends[[i]][, 2])])
          gappedDEM[(ends[[i]][, 1])[is.na(ends[[i]][, 2])]] <- minimum
        }
        extracts[[i]] <- raster::extract(gappedDEM, streams[i, ], along = T, cellnumbers = T)[[1]]
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
        if (
          any(is.na(temp[temp$group == j, ]$elevation)) &
            length(temp[temp$group == j, ]$elevation[!is.na(temp[temp$group == j, ]$elevation)]) >= 2) {
          temp[temp$group == j, ]$elevation <-
            stats::approx(temp[temp$group == j, ]$index,
              temp[temp$group == j, ]$elevation,
              xout = temp[temp$group == j, ]$index
            )$y
        }
      }

      # Update the raster
      gappedDEM[temp$cell] <- temp$elevation

      local <- raster::adjacent(gappedDEM, temp.na$cell, pairs = F, directions = 16)

      gappedDEM[local] <-
        gappedDEM[
          raster::cellFromXY(
            gappedDEM,
            t(
              apply(
                raster::xyFromCell(gappedDEM, local),
                1,
                function(x) {
                  maptools::nearestPointOnLine(streams[i, ]@lines[[1]]@Lines[[1]]@coords, x)
                }
              )
            )
          )
        ]

      gappedDEM <- raster::setMinMax(gappedDEM)

      finished.lines[i] <- TRUE
    }
    if (length(finished.lines[finished.lines == TRUE]) == start.length) done <- TRUE
  }

  if (any(!finished.lines)) {
    # Lakes where no interpolated values could be determined (usually
    # shallow stock-ponds) are re-filled with their original elevations.
    reservoirs <- sp::spTransform(reservoirs, CRS(projection(gappedDEM)))
    unfinished.streams <- sp::spTransform(Streams.gapped[!finished.lines], CRS(projection(gappedDEM)))
    if (any(is.na(c(reservoirs %over% unfinished.streams)))) {
      reservoirs.no.fill <- reservoirs[is.na(c(reservoirs %over% unfinished.streams)), ]
      reservoirs.no.fill.rast <- raster::extract(orig.dem, reservoirs.no.fill, cellnumbers = T)
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
    minimum <- suppressWarnings(min(raster[raster::adjacent(raster, cell, directions = neighbors)[, 2]], na.rm = T))
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
  gappedDEM.df <- data.frame(raster::coordinates(gappedDEM), raster::values(gappedDEM))
  names(gappedDEM.df) <- c("x", "y", "ELEVATION")
  gappedDEM.df.known <- gappedDEM.df[!is.na(gappedDEM.df$ELEVATION), ]
  gappedDEM.df.unknown <- gappedDEM.df[is.na(gappedDEM.df$ELEVATION), ]
  gappedDEM.idw.model <- gstat::gstat(
    id = "ELEVATION",
    formula = ELEVATION ~ 1,
    locations = ~ x + y,
    data = gappedDEM.df.known,
    nmax = 7,
    set = list(idp = 0.5)
  )
  gappedDEM.df.new <- predict(gappedDEM.idw.model, newdata = gappedDEM.df.unknown)
  gappedDEM.df.new$CELL <- raster::cellFromXY(gappedDEM, as.matrix(gappedDEM.df.new[, 1:2]))
  gappedDEM.final <- gappedDEM
  gappedDEM.final[gappedDEM.df.new$CELL] <- gappedDEM.df.new$ELEVATION.pred
  return(gappedDEM.final)
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
