library(magrittr)

mcphee <-
  sf::st_bbox(c(xmin = 707000,
                xmax = 725000,
                ymin = 4147000,
                ymax = 4165000),
              crs = 26912) %>%
  sf::st_as_sfc() %>%
  as("Spatial")

usethis::use_data(mcphee)
