library(dplyr)
library(ggplot2)
library(sf)

source(here::here("R", "00-settings.R"))

# @Check - are these current?
data("major", package = "PBSdata", envir = environment())
data("minor", package = "PBSdata", envir = environment())

data("locality.plus", package = "PBSdata", envir = environment())
data("locality", package = "PBSdata", envir = environment())

sb <- gfdata::load_survey_blocks(type = "polygon")

major_meta <- attributes(major)$PolyData
major_sf <- major |>
  dplyr::group_by(PID) |>
  dplyr::summarise(geometry = st_sfc(
    st_polygon(list(rbind(cbind(X, Y), c(X[1], Y[1]))))
  ), .groups = 'drop') |>
  st_as_sf(crs = st_crs(4326)) |>
  left_join(major_meta)
# Rowan already specifies labels in decent positions
major_centroids <- major_sf |>
  st_transform(st_crs(3156)) |> # no need to change s2 (sphere) to r2 (plane) i.e., `sf_use_s2(FALSE)`
  st_centroid()

minor_meta <- attributes(minor)$PolyData
minor_sf <- minor |>
  dplyr::group_by(PID) |>
  dplyr::summarise(
    geometry = st_sfc(
      st_polygon(list(rbind(cbind(X, Y), c(X[1], Y[1]))))
    ),
    .groups = "drop"
  ) |>
  st_as_sf(crs = st_crs(4326)) |>
  left_join(minor_meta)

ggplot() +
  geom_sf(data = sb, aes(fill = survey_abbrev), colour = NA, alpha = 0.5) +
  geom_sf(data = pacea::bc_coast) +
  geom_sf(data = major_sf, fill = NA) +
  geom_sf_label(data = major_sf, aes(label = label)) +
  gfplot::theme_pbs() +
  gfplot::coord_sf_auto(sb)

major_centroids

saveRDS(major_sf, file.path(data_raw_dir, "pfma.rds"))
