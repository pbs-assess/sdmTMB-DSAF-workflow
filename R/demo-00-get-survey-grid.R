library(dplyr)
library(ggplot2)

source(here::here("R", "00-settings.R"))

# 1. Load survey data (SYN WCVI dogfish sets; download-surveyjoin-data.R)
dog <- readRDS(file.path(data_gen_dir, "dogfish-syn-wcvi-sets.rds"))
utm_crs <- sdmTMB::get_crs(dog, ll_names = c("lon_start", "lat_start")) # get UTM crs for data region

# 2. Convert dataframe to sf object if not already; use equal distance projection (needed??)
dog_sf <-
  dog |>
  dplyr::select(year, lon_start, lat_start) |>
  dplyr::filter(!(is.na(lon_start) | is.na(lat_start))) |>
  sf::st_as_sf(coords = c("lon_start", "lat_start"), crs = 4326) |>
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
  dplyr::rename(geometry = x) |>
  # Cell id = position in the full grid, assigned before any clipping so it is
  # a stable key (one id per grid square, even if clipping splits a square
  # into several pieces) and is saved with the grid
  mutate(cell = row_number())


# 4. Overlay grid with sampling domain - 2 options:
  # a) intersect against sampled data points only
  # b) polygon covering full survey domain, clipped to the surveyed depth range

# 4a) intersect against sampled data points only  ------------------------------
intersect_points <- sf::st_intersects(full_grid, dog_sf)
point_intersect_id <- intersect_points |> purrr::map_dbl(length) > 0

point_grid <- full_grid[point_intersect_id, ]

ggplot(data = point_grid) +
  geom_sf() +
  geom_sf(data = dog_sf, alpha = 0.5) +
  theme_minimal()

# 4b) polygon covering full survey domain ------------------------------------
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

# 4c) clip the footprint to the surveyed GEBCO depth range ---------------------
# GEBCO (demo-00-get-GEBCO.R) is also the model's depth covariate
# (01-data-prep), so the grid and the covariate share one coastline: land and
# the nearshore drop out at the depth floor, and no cell is left with a depth
# <= 0 to be dropped later. The cap keeps prediction within the depths the
# survey sampled rather than extrapolating the depth effect down the slope.
gebco <- terra::rast(file.path(data_raw_dir, "gebco-depth.nc")) # m, positive down

# Limits from GEBCO at the sets (not recorded depth) so they're on the same
# scale as the clip. The deepest recorded sets (765-803 m) are 215-426 m in
# GEBCO: tow start positions on a steep slope.
set_depth_gebco <- terra::extract(
  gebco, terra::vect(sf::st_transform(dog_sf, terra::crs(gebco)))
)[, 2]
depth_limits <- c(floor(min(set_depth_gebco, na.rm = TRUE)), ceiling(max(set_depth_gebco, na.rm = TRUE)))
# depth_limits <- c(50, 1300) # or set by hand, e.g. the WCVI survey strata

# Work in GEBCO's lon/lat CRS (no raster resampling); crop first so polygonising is fast
gebco_crop <- terra::crop(gebco, terra::vect(sf::st_transform(dog_poly, terra::crs(gebco))), snap = "out")
depth_poly <- terra::ifel(gebco_crop >= depth_limits[1] & gebco_crop <= depth_limits[2], 1, NA) |>
  terra::as.polygons(dissolve = TRUE) |>
  sf::st_as_sf() |>
  sf::st_transform(utm_crs) |>
  sf::st_make_valid() |>
  sf::st_union()

# The raw clip has dozens of tiny pieces and holes, nearly all single GEBCO
# pixels (15 arc-sec, ~0.15 km2). Drop pieces and fill holes smaller than one
# grid cell.
min_area <- grid_spacing^2 # m^2
dog_poly_depth <- sf::st_intersection(dog_poly, depth_poly) |>
  sf::st_cast("POLYGON") |>
  (\(pieces) pieces[as.numeric(sf::st_area(pieces)) >= min_area])() |>
  purrr::map(\(p) {
    # ring 1 is the outer boundary; the rest are holes (planar area, UTM m)
    ring_area <- purrr::map_dbl(p, \(ring) sf::st_area(sf::st_polygon(list(ring))))
    sf::st_polygon(p[c(TRUE, ring_area[-1] >= min_area)])
  }) |>
  sf::st_sfc(crs = utm_crs) |>
  sf::st_union()

# st_intersection() cuts each cell's geometry down to the in-range area
poly_grid_clipped <- suppressWarnings(
  full_grid |>
    sf::st_intersection(dog_poly_depth) |>
    sf::st_collection_extract("POLYGON") # drop stray point/line slivers from edge-touching intersections
)

# Note that we now have grid cells with different areas
poly_grid_clipped |> mutate(area = sf::st_area(geometry)) |>
  arrange(area)

# official_grid <- gfdata::load_survey_blocks("synoptic") |> filter(survey_abbrev == "SYN WCVI")
# range(official_grid$depth_m, na.rm = TRUE)

ggplot() +
  geom_sf(data = pacea::bc_coast, fill = "grey85", colour = NA) +
  geom_sf(data = poly_grid_clipped) +
  geom_sf(data = dog_sf, size = 0.3) +
  geom_sf(data = dog_poly, fill = NA, colour = "black") +
  geom_sf(data = dog_poly_depth, fill = NA, colour = "red") +
  # geom_sf(data = official_grid, aes(fill = log10(depth_m)), alpha = 0.4) +
  labs(subtitle = sprintf("Footprint (black) clipped to GEBCO %g-%g m (red)", depth_limits[1], depth_limits[2])) +
  theme_minimal() +
  gfplot::coord_sf_auto(dog_sf)

# 5. Split cells at region boundaries ------------------------------------------
# Toy N/S split of WCVI (R/make-wcvi-region-polygons.R) standing in for real
# management units (e.g. NAFO divisions). Cells straddling a boundary get one
# row per region, each with its own recalculated area; pieces keep the parent
# `cell` id.
grid <- poly_grid_clipped

region_polys <- sf::st_read(
  file.path(data_raw_dir, "wcvi-region-polygons", "wcvi-region-polygons.shp"),
  quiet = TRUE
) |>
  sf::st_transform(utm_crs) |>
  select(region) # standardise the management-unit column name to "region"

grid_region <- suppressWarnings(
  sf::st_intersection(grid, region_polys) |>
    sf::st_collection_extract("POLYGON")
) |>
  mutate(area = as.numeric(sf::st_area(geometry)) / 1e6) |> # km^2
  filter(area > 1e-6) |> # drop GEOS sliver artifacts along the split line
  group_by(cell, region) |>
  summarise(area = sum(area), .groups = "drop")

# One X/Y per cell: the centroid of the kept area (the union of its pieces),
# shared by all pieces of a split cell. Covariates are extracted at this X/Y,
# so the halves of a split cell get identical covariate values and no
# artificial discontinuity at the region boundary.
cell_xy <- grid_region |>
  group_by(cell) |>
  summarise(.groups = "drop") |>
  sf::st_centroid() |>
  mutate(
    X = sf::st_coordinates(geometry)[, "X"] / 1000, # km
    Y = sf::st_coordinates(geometry)[, "Y"] / 1000
  ) |>
  sf::st_drop_geometry()
grid_region <- grid_region |>
  left_join(cell_xy, by = "cell") |>
  relocate(X, Y, .after = region)

# Area check: the split grid should cover the clipped grid (~0)
sum(as.numeric(sf::st_area(grid))) / 1e6 - sum(grid_region$area)
# Number of cells split across regions
grid_region |> sf::st_drop_geometry() |> count(cell) |> filter(n > 1) |> nrow()

ggplot() +
  geom_sf(data = grid_region, aes(fill = region), colour = NA) +
  geom_sf(data = grid_region |> group_by(region) |> summarise(), fill = NA, colour = "black") +
  theme_minimal() +
  gfplot::coord_sf_auto(grid_region)

# Only the region-split grid is saved; the unsplit grid is
# `grid_region |> group_by(cell) |> summarise(area = sum(area))`
saveRDS(grid_region, file.path(data_gen_dir, "dogfish-poly-grid.rds"))

