library(dplyr)
library(ggplot2)

source(here::here("R", "00-settings.R"))

# 1. Load survey data
dog <- sdmTMB::dogfish
utm_crs <- sdmTMB::get_crs(dog, ll_names = c("longitude", "latitude")) # get UTM crs for data region

# 2. Convert dataframe to sf object if not already; use equal distance projection (needed??)
dog_sf <-
  dog |>
  dplyr::select(year, longitude, latitude) |>
  dplyr::filter(!(is.na(longitude) | is.na(latitude))) |>
  sf::st_as_sf(coords = c('longitude', 'latitude'), crs = 'WGS84') |>
  sf::st_transform(crs = utm_crs)

# Optional: Can make survey grid using recent / modern sampling domain
# dog_sf <- dog_sf |> filter(year >= 2010)

ggplot(data = dog_sf) +
  geom_sf() +
  theme_minimal()

# 3. Make and overlay grid onto survey domain
# Use 2x2 km grid size in units of our polygon shape file
grid_spacing <- 2000 # maybe use 7km since that is ~3.5 nautical miles - I think used on east coast?

full_grid <- dog_sf |>
  sf::st_make_grid(cellsize = c(grid_spacing, grid_spacing)) |>
  sf::st_as_sf() |>
  dplyr::rename(geometry = x)


# 4. Overlay grid with sampling domain - 2 options:
  # a) intersect against sampled data points only
  # b) polygon covering full survey domain; land removed

# 4a) intersect against sampled data points only  ------------------------------
intersect_points <- sf::st_intersects(full_grid, dog_sf)
point_intersect_id <- intersect_points |> purrr::map_dbl(length) > 0

point_grid <- full_grid[point_intersect_id, ]

ggplot(data = point_grid) +
  geom_sf() +
  geom_sf(data = dog_sf, alpha = 0.5) +
  theme_minimal()

# 4b) polygon covering full survey domain; land removed ------------------------
dog_poly <- dog_sf |>
  sf::st_union() |>
  sf::st_concave_hull(ratio = 0.3, allow_holes = FALSE) |> # ratio adjust detail/corners (1 = convex hull)
  # sf::st_concave_hull(ratio = 0.01, allow_holes = TRUE) |> # Allows holes, more holes as ratio decreases
  sf::st_buffer(1000) # buffer so that grid extends 1km beyond points

intersect_poly <- sf::st_intersects(full_grid, dog_poly)
poly_intersect_id <- intersect_poly |> purrr::map_dbl(length) > 0

poly_grid <- full_grid[poly_intersect_id, ]

ggplot() +
  geom_sf(data = poly_grid) +
  geom_sf(data = dog_sf) +
  geom_sf(data = dog_poly, fill = NA) +
  theme_minimal()

# Remove land footprint. Here using `bcmaps::bc_bound()`
# bcmaps::available_layers()
# bcmaps::bc_bound_hres() # (if detailed coastal features needed; but very large object)
land_bc0 <- bcmaps::bc_bound() |>
  sf::st_transform(crs = utm_crs) |>
  sf::st_union()
land_bc <- land_bc0 |>
  sf::st_buffer(1000) # could apply buffer to land

dog_poly_land_clipped <- sf::st_difference(dog_poly, land_bc)

# st_intersection() cuts each cell's geometry down to the overwater area,
poly_grid_land_clipped <- suppressWarnings(
  full_grid |>
    sf::st_intersection(dog_poly_land_clipped) |>
    sf::st_collection_extract("POLYGON") # drop stray point/line slivers from edge-touching intersections
)

# Note that we now have grid cells with different areas
poly_grid_land_clipped |> mutate(area = sf::st_area(geometry)) |>
  arrange(area)

ggplot() +
  geom_sf(data = land_bc0) +
  geom_sf(data = poly_grid_land_clipped) +
  geom_sf(data = dog_sf) +
  geom_sf(data = dog_poly_land_clipped, fill = NA, colour = "red") +
  theme_minimal() +
  gfplot::coord_sf_auto(dog_sf)

saveRDS(poly_grid_land_clipped, file.path(data_gen_dir, "dogfish-poly-grid.rds"))

