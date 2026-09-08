library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(terra)

# @TODO: will abstract this as a function that takes a standardised (to be defined)
# - NetCDF input and survey grid OR survey data, to make it easy to swap/test
#   different covariate data sources in model fitting / downstream.
# - time-banding (not just year) in NetCDF
# - account for static (e.g., depth) and time-varying covariates (temperature, oxygen, ...)

source(here::here("R", "00-settings.R"))

# This demo attaches a covariate (bottom temperature, a NetCDF written by
# demo-00-get-GLORYS.R) onto a spatial prediction grid (from
# demo-00-get-survey-grid.R), for every grid cell x year combination.
#
# Extract with `terra::extract()`, matches the approach in 01-data-prep.qmd

# 1. Load prediction grid (2km cells over the dogfish survey domain, UTM 9N)
dog_grid <- readRDS(file.path(data_gen_dir, "dogfish-poly-grid.rds")) |>
  mutate(cell_id = dplyr::row_number())
st_crs(dog_grid)

grid_centroids <- sf::st_centroid(dog_grid)

# 2. Load bottom temperature raster (longitude x latitude x year bands, WGS84)
# Written as a NetCDF by demo-00-get-GLORYS.R (year as time step)
bt_stack <- terra::rast(file.path(data_gen_dir, "glorys-bottom-temp-june-yearly.nc"))
# NetCDF round-trips don't reliably preserve arbitrary layer `names()` -- derive
# them from the time metadata instead, which we set explicitly on write.
names(bt_stack) <- format(terra::time(bt_stack), "%Y")
terra::crs(bt_stack)

plot(bt_stack, col = cmocean::cmocean("thermal")(100))

# 3. Extract bottom temp at each grid cell centroid, for every year
grid_vect <- grid_centroids |>
  sf::st_transform(terra::crs(bt_stack)) |>
  terra::vect()

bt_extracted <-
  terra::extract(bt_stack, grid_vect, method = "simple") |> # value under cell returned
  # terra::extract(bt_stack, grid_vect, method = "bilinear") |> # value under cell bilinearly interpolated
  mutate(cell_id = grid_centroids$cell_id[ID]) |>
  select(-ID) |>
  tidyr::pivot_longer(-cell_id, names_to = "year", values_to = "bottom_temp") |>
  mutate(year = as.integer(year))

# 4. Join back onto the grid geometry to build the prediction grid
prediction_grid <- dog_grid |>
  inner_join(bt_extracted, by = "cell_id")

# 5. Visual check: does the extracted covariate line up spatially?
ggplot(prediction_grid |> filter(year %in% c(2006, 2016, 2026))) +
  geom_sf(aes(fill = bottom_temp), linewidth = 0.08, colour = "grey80") +
  facet_wrap(~year, ncol = 3) +
  scale_fill_gradientn(
    colours = cmocean::cmocean("thermal")(100), na.value = NA,
    guide = guide_colorbar(
      position = "top",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(10, "cm"),
      barheight = unit(0.4, "cm"),
      direction = "horizontal"
    )
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  labs(fill = "Bottom temperature (°C)") +
  gfplot::coord_sf_auto(dog_grid)

saveRDS(prediction_grid, file.path(data_gen_dir, "dogfish-prediction-grid-GLORYS.rds"))

# Comparison with terra::extract() using bilinear interpolation:
bt_extracted_bi <-
  terra::extract(bt_stack, grid_vect, method = "bilinear") |> # value under cell bilinearly interpolated
  mutate(cell_id = grid_centroids$cell_id[ID]) |>
  select(-ID) |>
  tidyr::pivot_longer(-cell_id, names_to = "year", values_to = "bottom_temp") |>
  mutate(year = as.integer(year))

prediction_grid_bi <- dog_grid |>
  inner_join(bt_extracted_bi, by = "cell_id")

bind_rows(
   "no interpolation" = prediction_grid,
   "bilinear interpolation" = prediction_grid_bi,
   .id = "interpolation"
 ) |>
   filter(year == 2026) |>
  mutate(interpolation = forcats::fct_inorder(interpolation)) |>
ggplot() +
  geom_sf(aes(fill = bottom_temp), linewidth = 0.08, colour = "grey80") +
  facet_wrap(~interpolation, ncol = 2) +
  scale_fill_gradientn(
    colours = cmocean::cmocean("thermal")(100), na.value = NA,
    guide = guide_colorbar(
      position = "top",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(10, "cm"),
      barheight = unit(0.4, "cm"),
      direction = "horizontal"
    )
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  gfplot::coord_sf_auto(dog_grid)
