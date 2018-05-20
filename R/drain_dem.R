#' Drain Your Dammed DEM
#'
#' Given a digital elevation model (DEM) with human-made dammed reservoirs,
#' returns a version with estimates for elevations under reservoirs and dams.
#'
#' @param x A [raster::raster] DEM. Alternatively, a [sf::st_bbox] object for which a DEM will
#' be downloaded
#' @param ... Arguments passed on to various [FedData] functions.
#'
#' @return A drained [raster::raster] DEM.
#' @importFrom magrittr %>% %<>% %$%
#' @export
#'
#' @examples
#' library(sf)
#' library(demdrainer)
#' x <- mcor::mt_counties %>% dplyr::filter(County == "Missoula") %>% sf::st_bbox()
drain_dem <-
  function(x,
             label,
             ...) {
    if (missing(label)) {
      label <- "test"
    }

    if (class(x) == "bbox") {
      x %<>%
        FedData::get_ned(
          template = x %>%
            sf::st_as_sfc() %>%
            as("Spatial"),
          label = label
        )
    }

    nhd <- FedData::get_nhd(x,
                            label = label)

    ssurgo <- FedData::get_ssurgo(x,
                                  label = label)


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
