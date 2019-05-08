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
                            label = label,
                            ...)

    ssurgo <- FedData::get_ssurgo(x,
                                  label = label,
                                  ...)

    Dams <- ssurgo$spatial %>%
      sf::st_as_sf() %>%
      dplyr::left_join(ssurgo$tabular$mapunit %>%
                         dplyr::mutate(MUKEY = as.factor(mukey))) %>%
      dplyr::filter(muname == "Dam") %>%
      sf::st_buffer(dist = 1 * res(x)[1]) %>%
      sf::st_transform(projection(x))

    Streams <- nhd$`_Flowline` %>%
      sf::st_as_sf() %>%
      dplyr::filter(FCode %in% c(55800, 46006, 46003)) %>%
      sf::st_transform(projection(x))

    Reservoirs <- nhd$`_Waterbody` %>%
      sf::st_as_sf() %>%
      dplyr::filter(FCode %in% c(39004, 39009, 39010, 39001, 39012)) %>%
      # tibble::as_tibble() %>%
      # dplyr::bind_rows(nhd$`_Area` %>%
      #                    sf::st_as_sf() %>%
      #                    tibble::as_tibble()) %>%
      # sf::st_as_sf() %>%
      dplyr::filter(AreSqKm > 0.16) %>%
      lwgeom::st_make_valid() %>%
      sf::st_transform(projection(x))

    # Dams %<>%
    #   as("Spatial")
    #
    # Streams %<>%
    #   as("Spatial")
    #
    # Reservoirs %<>%
    #   as("Spatial")

    reservoirsWithDams <- st_is_within_distance(Reservoirs,
                                                Dams,
                                                dist = 2 * res(x)[1],
                                                sparse = FALSE) %>%
      t()

    reservoirsWithDams.vector <- vector("logical", nrow(Reservoirs))
    for (i in 1:nrow(Reservoirs)) {
      reservoirsWithDams.vector[i] <- any(reservoirsWithDams[, i])
    }

    # Reservoirs %<>%
    #   dplyr::filter(reservoirsWithDams.vector)

    damsWithReservoirs <- st_is_within_distance(Dams,
                                                Reservoirs,
                                                dist = 100,
                                                sparse = FALSE) %>%
      t()

    damsWithReservoirs.vector <- vector("logical", nrow(Dams))
    for (i in 1:nrow(Dams)) {
      damsWithReservoirs.vector[i] <- any(damsWithReservoirs[, i])
    }

    Dams %<>%
      dplyr::filter(damsWithReservoirs.vector)

    Streams %<>%
      dplyr::select(GNIS_Nm)

    Streams.named.merged <- Streams %>%
      dplyr::filter(!is.na(GNIS_Nm),
                    GNIS_Nm != "<Null>") %>%
      dplyr::group_by(GNIS_Nm) %>%
      dplyr::summarise() %>%
      dplyr::ungroup() %>%
      sf::st_cast() %>%
      sf::st_line_merge()

    Streams.unnamed.merged <- Streams %>%
      dplyr::filter(is.na(GNIS_Nm) |
                    GNIS_Nm == "<Null>")  %>%
      dplyr::mutate(GNIS_Nm = "Unnamed") %>%
      dplyr::group_by(GNIS_Nm) %>%
      dplyr::summarise() %>%
      dplyr::ungroup() %>%
      sf::st_cast() %>%
      sf::st_line_merge()

    Streams <- rbind(Streams.named.merged,Streams.unnamed.merged)

    # Streams <- spRbind(Streams.named.merged, Streams.unnamed.merged)
    # Streams <- mergeStreams(Streams, projection(x))
    #
    # Streams %<>%
    #   dplyr::select(GNIS_Nm) %>%
    #   dplyr::group_by(GNIS_Nm) %>%
    #   dplyr::summarise() %>%
    #   dplyr::ungroup() %>%
    #   sf::st_cast() %>%
    #   sf::st_line_merge()
    #
    # %>%
    #   as("Spatial") %>%
    #   mergeStreams(projection(x))
    #
    #
    # Streams.named <- Streams[!is.na(Streams$GNIS_Nm), ]
    # Streams.unnamed <- Streams[is.na(Streams$GNIS_Nm), ]
    # Streams.named.merged <- gLineMerge(Streams.named, byid = T, id = Streams.named$GNIS_Nm)
    # Streams.unnamed.merged <- gLineMerge(Streams.unnamed)
    # Streams <- spRbind(Streams.named.merged, Streams.unnamed.merged)
    # Streams <- mergeStreams(Streams, projection(x))

    # Join the dam and water polygons, and rasterize
    bad.data.vect <- sf::st_union(Dams, Reservoirs) %>%
      sf::st_cast() %>%
      dplyr::select()

    res.elevs.rast <- raster::mask(x,
                                   fasterize::fasterize(bad.data.vect, x),
                                   inverse = TRUE)

    bad.data.vect %<>%
      sf::st_union()

    # Find which streams have missing elevation data
    Streams.gapped <- Streams %>%
      sf::st_intersection(bad.data.vect %>%
                            sf::st_buffer(dist = 0.002))

    Streams.gapped %<>%
      dplyr::filter(GNIS_Nm == "Unnamed") %>%
      st_cast('LINESTRING', do_split = TRUE) %>%
      rbind(Streams.gapped %>%
              dplyr::filter(GNIS_Nm != "Unnamed"))

    Streams.gapped %<>%
      st_cast('MULTILINESTRING') %>%
      st_cast('LINESTRING', do_split = TRUE)

    res.elevs.rast <- bootstrapDrainDEM(gappedDEM = res.elevs.rast,
                                        streams = Streams.gapped %>%
                                          as("Spatial"),
                                        reservoirs = Reservoirs %>%
                                          as("Spatial"),
                                        orig.dem = x)

    dem.final <- fillIDW(res.elevs.rast)

    # A 5x5 mean smooth
    # dem.final.smoothed <- focal(dem.final, w = matrix(1, nrow = 5, ncol = 5), fun = mean, na.rm = T, pad = T)

    # return(dem.final.smoothed)
    return(dem.final)
  }
