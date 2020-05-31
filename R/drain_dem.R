#' Drain Your Dammed DEM
#'
#' Given a digital elevation model (DEM) with human-made dammed reservoirs,
#' returns a version with estimates for elevations under reservoirs and dams.
#'
#' @param x A [raster::raster] DEM. Alternatively, a [sf::st_bbox] object for which a DEM will
#' be downloaded
#' @param label A character string naming the study area.
#' @param ... Arguments passed on to various `FedData` functions.
#'
#' @return A drained [raster::raster] DEM.
#' @importFrom magrittr %>% %<>% %$%
#' @export
#'
#' @examples
#' library(magrittr)
#' library(demdrainer)
#'
#' dem <-
#'   FedData::get_ned(demdrainer::mcphee,
#'     label = "McPhee"
#'   )
#'
#' drained_dem <-
#'   drain_dem(dem,
#'     label = "McPhee"
#'   )
#'
#' raster::plot(dem)
#' raster::plot(drained_dem)
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
            methods::as("Spatial"),
          label = label
        )
    }

    nhd <- FedData::get_nhd(x,
                            label = label,
                            ...
    )

    ssurgo <- FedData::get_ssurgo(x,
                                  label = label,
                                  ...
    )

    Dams <-
      ssurgo$spatial %>%
      sf::st_as_sf() %>%
      dplyr::left_join(ssurgo$tabular$mapunit %>%
                         dplyr::mutate(MUKEY = as.factor(mukey))) %>%
      dplyr::filter(muname == "Dam") %>%
      sf::st_buffer(dist = 1 * raster::res(x)[1]) %>%
      sf::st_transform(raster::projection(x))

    Streams <-
      nhd$`Flowline` %>%
      sf::st_as_sf() %>%
      dplyr::filter(FCode %in% c(55800, 46006, 46003)) %>%
      sf::st_transform(raster::projection(x))

    Reservoirs <-
      nhd$`Waterbody` %>%
      sf::st_as_sf() %>%
      dplyr::filter(FCode %in% c(39004, 39009, 39010, 39001, 39012)) %>%
      dplyr::filter(AreaSqKm > 0.16) %>%
      sf::st_make_valid() %>%
      sf::st_transform(raster::projection(x))

    reservoirsWithDams <- sf::st_is_within_distance(Reservoirs,
                                                Dams,
                                                dist = 2 * raster::res(x)[1],
                                                sparse = FALSE
    ) %>%
      t()

    reservoirsWithDams.vector <- vector("logical", nrow(Reservoirs))
    for (i in 1:nrow(Reservoirs)) {
      reservoirsWithDams.vector[i] <- any(reservoirsWithDams[, i])
    }

    damsWithReservoirs <- sf::st_is_within_distance(Dams,
                                                Reservoirs,
                                                dist = 100,
                                                sparse = FALSE
    ) %>%
      t()

    damsWithReservoirs.vector <- vector("logical", nrow(Dams))
    for (i in 1:nrow(Dams)) {
      damsWithReservoirs.vector[i] <- any(damsWithReservoirs[, i])
    }

    Dams %<>%
      dplyr::filter(damsWithReservoirs.vector)

    Streams %<>%
      dplyr::select(GNIS_Name)

    Streams.named.merged <- Streams %>%
      dplyr::filter(
        !is.na(GNIS_Name),
        GNIS_Name != "<Null>"
      ) %>%
      dplyr::group_by(GNIS_Name) %>%
      dplyr::summarise() %>%
      dplyr::ungroup() %>%
      sf::st_cast() %>%
      sf::st_line_merge()

    Streams.unnamed.merged <- Streams %>%
      dplyr::filter(is.na(GNIS_Name) |
                      GNIS_Name == "<Null>") %>%
      dplyr::mutate(GNIS_Name = "Unnamed") %>%
      dplyr::group_by(GNIS_Name) %>%
      dplyr::summarise() %>%
      dplyr::ungroup() %>%
      sf::st_cast() %>%
      sf::st_line_merge()

    Streams <- rbind(Streams.named.merged, Streams.unnamed.merged)

    # Join the dam and water polygons, and rasterize
    bad.data.vect <- sf::st_union(Dams, Reservoirs) %>%
      sf::st_cast() %>%
      dplyr::select()

    res.elevs.rast <- raster::mask(x,
                                   fasterize::fasterize(bad.data.vect, x),
                                   inverse = TRUE
    )

    bad.data.vect %<>%
      sf::st_union()

    # Find which streams have missing elevation data
    Streams.gapped <- Streams %>%
      sf::st_intersection(bad.data.vect %>%
                            sf::st_buffer(dist = 0.002))

    Streams.gapped %<>%
      dplyr::filter(GNIS_Name == "Unnamed") %>%
      sf::st_cast("LINESTRING", do_split = TRUE) %>%
      rbind(Streams.gapped %>%
              dplyr::filter(GNIS_Name != "Unnamed"))

    Streams.gapped %<>%
      sf::st_cast("MULTILINESTRING") %>%
      sf::st_cast("LINESTRING", do_split = TRUE)

    res.elevs.rast <- bootstrapDrainDEM(
      gappedDEM = res.elevs.rast,
      streams = Streams.gapped %>%
        methods::as("Spatial"),
      reservoirs = Reservoirs %>%
        methods::as("Spatial"),
      orig.dem = x
    )

    dem.final <- fillIDW(res.elevs.rast)

    dem.final[dem.final > x] <- x[dem.final > x]

    return(dem.final)
  }
