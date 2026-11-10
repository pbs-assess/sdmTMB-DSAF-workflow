# pak::pak("ncdf4")  # OPeNDAP reads (sec 3)
# Also called below via `pkg::`: here, cli, purrr, terra, tidyterra, pacea, gfplot
library(sf)
library(ggplot2)

source(here::here("R", "00-settings.R"))
dir.create(data_raw_dir, showWarnings = FALSE, recursive = TRUE)
output_filepath <- file.path(data_raw_dir, "gebco-depth.nc")

# This script demonstrates how to pull a regional subset of the GEBCO global
# bathymetry grid and save it as a depth covariate (`depth_m`, m, positive
# down) that attach_covariate() can use. GEBCO is static, so the file is a
# single layer with no time.
# No login is needed: GEBCO is served openly over OPeNDAP by CEDA, which lets
# us read only the rows/columns we need instead of the full ~7.5 GB global file.
# https://www.gebco.net/data-products/gridded-bathymetry-data (OPeNDAP section)
#
# (A packaged version with bbox-coverage checks lives in download_gebco() in
# NCAM/analysis/spatial_model/R/covariate_download_helpers.R, in the NCAM
# project.)

gebco_url <- "https://dap.ceda.ac.uk/thredds/dodsC/bodc/gebco/global/gebco_2026/ice_surface_elevation/netcdf/GEBCO_2026.nc"

# --- 1. Spatial extent ------------------------------------------------------
# Use the set locations: the survey grid (demo-00-get-survey-grid.R) is built
# from the sets and clipped with this GEBCO file, so run this first. The grid
# hull extends 1 km beyond the sets; the 0.25 deg buffer covers that and keeps
# edge cells and bilinear extraction inside the raster.
grid <- readRDS(file.path(data_gen_dir, "dogfish-syn-wcvi-sets.rds")) |>
  dplyr::filter(!is.na(lon_start), !is.na(lat_start)) |>
  st_as_sf(coords = c("lon_start", "lat_start"), crs = 4326)

region <- grid |>
  st_transform(4326) |> # the GEBCO grid is lon/lat, so the bbox must be too
  st_bbox() |>
  unclass()
region <- region + c(-0.25, -0.25, 0.25, 0.25)

# --- 2. Open the remote file and find the rows/columns in the region ---------
# Only the coordinate vectors are downloaded here (~86,400 lon + 43,200 lat).
nc <- ncdf4::nc_open(gebco_url)
nc # variables: lat, lon, elevation (m, positive up)

lon <- ncdf4::ncvar_get(nc, "lon") # cell centres, west to east
lat <- ncdf4::ncvar_get(nc, "lat") # cell centres, south to north

i_lon <- which(lon >= region["xmin"] & lon <= region["xmax"])
i_lat <- which(lat >= region["ymin"] & lat <= region["ymax"])
c(n_lon = length(i_lon), n_lat = length(i_lat)) # 240 cells per degree (15 arc-sec)

# --- 3. Download elevation in latitude bands --------------------------------
# The server rejects requests that are too large, so read ~1000 latitude rows
# at a time (~10 MB per band for a few degrees of longitude) and stitch them
# together. Note ncdf4 reverses the dimension order: arrays come back
# [lon, lat], and `start`/`count` are given as (lon, lat).
if (!file.exists(output_filepath)) {
  bands <- split(i_lat, ceiling(seq_along(i_lat) / 1000))

  elev <- bands |>
    purrr::map(\(rows) {
      ncdf4::ncvar_get(
        nc, "elevation",
        start = c(min(i_lon), min(rows)),
        count = c(length(i_lon), length(rows))
      )
    }, .progress = TRUE) |>
    purrr::reduce(cbind) # [lon, lat]; lat columns south to north

  # --- 4. Convert to a north-up raster of depth -------------------------------
  # Transpose to [lat, lon] and flip the rows so the first row is the
  # northernmost, which is the layout terra expects for a matrix.
  elev_matrix <- t(elev)[rev(seq_along(i_lat)), ]

  # GEBCO coordinates are cell centres; the extent needs cell edges, so pad by
  # half a cell (1/240 deg).
  half_cell <- 1 / 240 / 2
  gebco <- terra::rast(
    elev_matrix,
    crs = "EPSG:4326",
    extent = terra::ext(
      range(lon[i_lon]) + c(-half_cell, half_cell),
      range(lat[i_lat]) + c(-half_cell, half_cell)
    )
  )

  # GEBCO is elevation (positive up). Flip the sign so depth is positive down,
  # which log(depth_m) downstream needs. Land becomes negative depth, and
  # 01-data-prep gives those cells NA log depth.
  depth <- -gebco # This the standard way we will be saving and using the GEBCO data

  # NetCDF with a named variable, to match the covariate standard (see the
  # notes above attach_covariate() in functions.R)
  terra::writeCDF(
    depth, output_filepath,
    varname = "depth_m", longname = "Depth (GEBCO_2026, 15 arc-sec)", unit = "m",
    overwrite = TRUE
  )
  cli::cli_inform(c("v" = "Saved GEBCO depth to {.file {output_filepath}}."))
} else {
  cli::cli_inform(c(
    "i" = "{.file {output_filepath}} already exists; skipping download. Delete it and re-run to refresh (e.g. for a larger region)."
  ))
}

ncdf4::nc_close(nc)

# --- 5. Read back with terra and inspect -------------------------------------
depth <- terra::rast(output_filepath)
depth
terra::res(depth) * 3600 # 15 arc-sec
names(depth)             # "depth_m"

# Sign check: the deep offshore cells should be large positive numbers and
# land should be negative
terra::global(depth, range, na.rm = TRUE)

ggplot() +
  tidyterra::geom_spatraster(data = depth) +
  geom_sf(data = pacea::bc_coast) +
  geom_sf(data = st_transform(grid, 4326), fill = NA, linewidth = 0.1) +
  scale_fill_gradientn(
  colours = hcl.colors(6, "Blues 3", rev = TRUE),
  transform = "log10",
  limits = c(20, 1000), oob = scales::squish, na.value = NA,
  breaks = c(20, 50, 100, 200, 500, 1000),
  name = "Depth (m)"
) +
  gfplot::coord_sf_auto(grid)

# --- 6. Extraction method: "simple" vs "bilinear" -----------------------------
# This is not necessary for the GEBCO download, but I wanted to have this
# comparison somewhere for now.
# attach_covariate() uses terra::extract()
#   method = "simple": returned value is what is directly under the point
#   method = "bilinear": returned value is the weighted average of the 4 nearest cell centres
# Note - when a point falls exactly at the centre of a pixel, “bilinear” gives the same result as the “simple” method
# Here is a comparison of both methods relative to depths recorded at each set
# (`depth_m`, measured in situ).
# Takeaway - bilinear is probably fine.
# Caveat: the recorded depth describes the whole tow, but we extract at the
# tow start. #@question: True?
source(here::here("R", "functions.R")) # attach_covariate()

sets_depth <- readRDS(file.path(data_gen_dir, "dogfish-syn-wcvi-sets.rds")) |>
  st_as_sf(coords = c("lon_start", "lat_start"), crs = 4326) |>
  attach_covariate(output_filepath, "gebco_simple", method = "simple") |>
  attach_covariate(output_filepath, "gebco_bilinear", method = "bilinear") |>
  st_drop_geometry() |>
  tidyr::pivot_longer(
    c(gebco_simple, gebco_bilinear),
    names_to = "method", names_prefix = "gebco_", values_to = "gebco_depth"
  )

sets_depth |>
  dplyr::summarise(
    rmse = sqrt(mean((gebco_depth - depth_m)^2, na.rm = TRUE)),
    median_abs_error = median(abs(gebco_depth - depth_m), na.rm = TRUE),
    mean_bias = mean(gebco_depth - depth_m, na.rm = TRUE),
    r = cor(gebco_depth, depth_m, use = "complete.obs"),
    .by = method
  )
# SYN WCVI (1,616 sets), as of Nov 2026:
#   method    rmse  median_abs_error  mean_bias     r
#   simple    27.4  4.0               2.6       0.958
#   bilinear  26.0  4.1               2.6       0.962
# Very similar: a typical set is within ~4 m either way. Bilinear does a bit
# better on the larger misses (lower RMSE) on the steep shelf break, where a
# 15 arc-sec cell spans a big range of depths.

ggplot(sets_depth, aes(depth_m, gebco_depth)) +
  geom_abline(linetype = "dashed") +
  geom_point(alpha = 0.3, size = 0.8) +
  facet_wrap(~method) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Recorded set depth (m)", y = "GEBCO depth at tow start (m)")

# --- 7. Using it downstream ---------------------------------------------------
# 01-data-prep-v2.qmd reads this file as depth_file and extracts it with
# method = "bilinear" at both the sets (for the fit) and the grid cell centroids
# (for prediction), so the model sees one depth source throughout. The
# comparison above shows how far that is from the recorded set depth.

# --- 8. Cite --------------------------------------------------------------
# GEBCO Compilation Group (2026) GEBCO 2026 Grid. See gebco.net for the DOI.
