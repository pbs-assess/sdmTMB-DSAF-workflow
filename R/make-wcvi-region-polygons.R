library(dplyr)
library(ggplot2)
library(sf)

source(here::here("R", "00-settings.R"))

# WCVI survey blocks define the extent to split - keeps the two rectangles
# tight to the actual survey grid instead of an arbitrary bounding box.
sb <- gfdata::load_survey_blocks(type = "polygon")
wcvi_blocks <- sb |> filter(survey_abbrev == "SYN WCVI", active_block)

bbox <- sf::st_bbox(wcvi_blocks)

# Buffer the outer extent so the polygons extend beyond grid.
# Units (m, UTM 9N) match the blocks' native CRS.
buffer_m <- 10000
xmin <- bbox["xmin"] - buffer_m
xmax <- bbox["xmax"] + buffer_m
ymin <- bbox["ymin"] - buffer_m
ymax <- bbox["ymax"] + buffer_m

# Arbitrary N/S split
split_northing <- 5420000

wcvi_south <- sf::st_polygon(list(rbind(
  c(xmin, ymin),
  c(xmax, ymin),
  c(xmax, split_northing),
  c(xmin, split_northing),
  c(xmin, ymin)
)))

wcvi_north <- sf::st_polygon(list(rbind(
  c(xmin, split_northing),
  c(xmax, split_northing),
  c(xmax, ymax),
  c(xmin, ymax),
  c(xmin, split_northing)
)))

wcvi_regions <- sf::st_sf(
  region = c("WCVI N", "WCVI S"),
  geometry = sf::st_sfc(wcvi_north, wcvi_south, crs = sf::st_crs(wcvi_blocks))
)

ggplot() +
  geom_sf(data = wcvi_blocks, fill = NA, colour = "grey70") +
  geom_sf(data = wcvi_regions, aes(fill = region), alpha = 0.3) +
  gfplot::theme_pbs() +
  gfplot::coord_sf_auto(wcvi_blocks)

# Shapefiles are a multi-file format (.shp/.shx/.dbf/.prj) - give this one its
# own directory so the files stay bundled together in data-raw/.
# Best mimics data we would get.
wcvi_regions_dir <- file.path(data_raw_dir, "wcvi-region-polygons")
dir.create(wcvi_regions_dir, showWarnings = FALSE, recursive = TRUE)

sf::st_write(
  wcvi_regions,
  file.path(wcvi_regions_dir, "wcvi-region-polygons.shp"),
  delete_dsn = TRUE,
  quiet = TRUE
)
