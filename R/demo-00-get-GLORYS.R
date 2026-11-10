# pak::pak("CopernicusMarine")       # CRAN release
# pak::pak("CFtime")                 # stars needs it to parse NetCDF time stamps; without it time comes back as raw numbers
# pak::pak(c('keyring', 'askpass'))  # secure password storage (sec 0b)
# pak::pak('ncmeta')                 # optional: not used here, but handy for inspecting .nc files (ncmeta::nc_vars())
# Also called below via `pkg::`: here, cli, gfplot, tidyterra, pacea, cmocean
library(CopernicusMarine)
library(stars)   # only needed for the download/write step below (sec 3): write_mdim() keeps the time dimension
library(sf)
library(ggplot2)
library(terra)   # attached so plot() dispatches to SpatRaster methods; other calls use terra::

source(here::here("R", "00-settings.R"))
output_filepath <- file.path(data_raw_dir, "glorys-bottom-temp-monthly.nc")

# This script demonstrates how you can query Copernicus and get bottom temperature
# The monthly file it writes is used as-is downstream: attach_covariate() picks
# the layer for each set's date (date_col =) and, on the prediction grid, the
# chosen month of each year (month =; several months are averaged).
# The extent and start year come from the survey sets data.

# ==========================================================================
# ONE-TIME SETUP -- run these lines once, interactively, in the console.
# Do NOT run automatically as part of the vignette/script. Comment back out
# once done
# ==========================================================================

# --- 0a. Store your username as a non-secret environment variable -------
# usethis::edit_r_environ() # or just open your .Renviron and edit by hand yourself
# # In the file that opens, add this line, save, then restart R:
# # COPERNICUS_USERNAME=your_username

# --- 0b. Store your password securely in the OS keychain ----------------
# A dialogue box will pop up asking for your password.
# keyring::key_set(service = "copernicus_marine", username = Sys.getenv("COPERNICUS_USERNAME"))
# NOTE - this is much safer than just hard coding your password somewhere.

# ==========================================================================
# EVERY SESSION -- after setting up credentials you can do the following
# ==========================================================================
# --- 1. Authenticate ------------------------------------------------------
username <- Sys.getenv("COPERNICUS_USERNAME")

cms_set_username(username)
cms_set_password(
  keyring::key_get(service = "copernicus_marine", username = username)
)
cms_login()  # errors if authentication fails; otherwise returns your account info

# --- 2. (Optional) confirm the variable name --------------
# GLORYS bottom temperature is in the multi-year reanalysis product
# (GLOBAL_MULTIYEAR_PHY_001_030), which has daily and monthly layers.
# Variable naming can shift between product revisions, so it's worth
# checking the product's variable list once before relying on "bottomT":
# vignette("product-info")  # or browse the product page directly


# --- 3. Download a bottom-temperature subset for a region -----------------
# Can load survey data to extract relevant parameters (time and spatial extent)
dog <- readRDS(file.path(data_gen_dir, "dogfish-syn-wcvi-sets.rds"))
dog_sf <- st_as_sf(dog, coords = c("lon_start", "lat_start"), crs = 4326) # lon/lat: the degree buffer below relies on this

# Spatial extent ------
# Load data domain or manually specify region in this order: c(xmin, ymin, xmax, ymax)
# Use your prediction grid's bbox (in lon/lat) if you have one; the set
# locations alone can miss outer grid cells. The 0.25 deg buffer (3 GLORYS
# cells) keeps edge cells and bilinear extraction inside the raster.
region <- st_bbox(dog_sf) |>
  unclass()
region <- region + c(-0.25, -0.25, 0.25, 0.25)

# Layer -----
# Monthly means rather than daily: daily over a large region and 30+ years is
# several GB, and monthly is enough for yearly or seasonal covariates.
# Daily:   "cmems_mod_glo_phy_my_0.083deg_P1D-m"
# Monthly: "cmems_mod_glo_phy_my_0.083deg_P1M-m"
copernicus_product <- "GLOBAL_MULTIYEAR_PHY_001_030"
copernicus_layer   <- "cmems_mod_glo_phy_my_0.083deg_P1M-m" # for cod - monthly seemed good enough

# Temporal extent -----
# Read the date range of the layer (not the whole product -- daily and monthly
# layers can end on different dates) from the STAC catalogue; no login needed.
# Catalogue IDs carry a version suffix (e.g. "_202311"), so match by prefix.
product_meta <- cms_product_metadata(copernicus_product)
layer_idx <- which(startsWith(product_meta$id, paste0(copernicus_layer, "_")))
if (length(layer_idx) != 1) {
  cli::cli_abort("Layer {.val {copernicus_layer}} not found (or not unique) in {.val {copernicus_product}} metadata.")
}
layer_extent <- unlist(product_meta$properties[[layer_idx]][["cube:dimensions"]]$time$extent)
layer_extent

layer_start <- as.Date(substr(layer_extent[[1]], 1, 10)) # ISO 8601 timestamps -- the date is always the first 10 characters
layer_end   <- as.Date(substr(layer_extent[[2]], 1, 10))

# Restrict to what the layer actually covers
# start date: January of the first survey year, but no earlier than the layer
#             start (1993)
# end date:   lags behind Sys.Date() by a processing delay.
# Downloading to the latest month (not just the last survey year) leaves the
# most recent conditions available as a starting point for forecasts.
start_date <- format(max(as.Date(sprintf("%d-01-01", min(dog$year))), layer_start), "%Y-%m-%d")
end_date   <- format(min(Sys.Date(), layer_end), "%Y-%m-%d")

# Download (via CopernicusMarine, which returns a stars object) and cache to
# disk as a proper multidimensional NetCDF (write_mdim() -- write_stars()
# flattens extra dimensions into bands and drops the time coordinate).
# Safer to write to mdim and read back in using terra rather than converting to
# SpatRaster in memory because the `terra::rast()` is picky about dimension
# regularity/naming.
# Allows for handling as terra object for rest of workflow.
# This downloads everything in one file and skips if it already exists; to
# add new months later, delete the file and re-run. (For one-file-per-year
# downloads that only fetch missing years, see download_glorys() in
# NCAM/analysis/spatial_model/R/covariate_download_helpers.R, in the NCAM
# project)
if (!file.exists(output_filepath)) {
  bottom_temp <- cms_download_subset(
    product       = copernicus_product,
    layer         = copernicus_layer,
    variable      = "bottomT",
    region        = region,
    timerange     = c(start_date, end_date),
    progress      = TRUE
  )
  # Expect GDAL warnings ("Dimension ... lacks a indexing variable",
  # "SetIndexingVariable() not implemented"): the netCDF writer can't formally
  # link dimensions to coordinate variables, but it still writes them under the
  # dimension names, which CF readers (terra, stars) treat as coordinates.
  write_mdim(bottom_temp["bottomT"], output_filepath)
  cli::cli_inform(c(
    "v" = "Saved monthly GLORYS bottom temperature ({start_date} to {end_date}) to {.file {output_filepath}}."
  ))
} else {
  cli::cli_inform(c(
    "i" = "{.file {output_filepath}} already exists; skipping download. Delete it and re-run to refresh."
  ))
}

# --- 4. Read back with terra and inspect -------
bt_monthly <- terra::rast(output_filepath)

# Fail fast if GDAL/terra didn't parse the CF time/CRS metadata cleanly
if (anyNA(terra::time(bt_monthly)) || !nzchar(terra::crs(bt_monthly))) {
  cli::cli_abort("{.file {output_filepath}} did not parse valid time/CRS metadata into terra.")
}

bt_monthly

plot(bt_monthly[[1]], col = hcl.colors(100, "RdYlBu", rev = TRUE), axes = TRUE) # base R option

# Aggregate over time (e.g. mean), then plot:
bt_mean <- terra::app(bt_monthly, mean, na.rm = TRUE)
plot(bt_mean, col = hcl.colors(100, "RdYlBu", rev = TRUE), axes = TRUE)

ggplot() +
  tidyterra::geom_spatraster(data = bt_monthly[[1]]) +
  scale_fill_gradientn(colours = hcl.colors(100, "RdYlBu", rev = TRUE), na.value = NA)

# --- 5. Plot June of each year -------------------------------------------
bt_june <- bt_monthly[[format(terra::time(bt_monthly), "%m") == "06"]]
names(bt_june) <- format(terra::time(bt_june), "%Y")

ggplot() +
  geom_sf(data = pacea::bc_coast) +
  tidyterra::geom_spatraster(data = bt_june) +
  facet_wrap(~ lyr, ncol = 4) +
  scale_fill_gradientn(colours = cmocean::cmocean("thermal")(100), na.value = NA) +
  gfplot::coord_sf_auto(dog_sf) +
  theme(axis.text = element_blank())

# --- 6. Cite the product (required under CMEMS terms of use) --------------
# Public catalogue lookup, so no login needed. Despite the name, `$doi` holds
# the full citation text (with the DOI at the end).
# cms_cite_product("GLOBAL_MULTIYEAR_PHY_001_030")$doi
# [1] "E.U. Copernicus Marine Service Information; Global Ocean Physics Reanalysis - GLOBAL_MULTIYEAR_PHY_001_030 (2012-08-26). DOI:10.48670/moi-00021"
