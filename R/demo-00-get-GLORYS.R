# pak::pak("CopernicusMarine")  # CRAN release
# pak::pak("CFtime")
# pak::pak(c('keyring', 'askpass'))
# pak::pak('ncmeta')
library(CopernicusMarine)
library(stars)   # only needed for the download/write step below (sec 3)
library(sf)
library(ggplot2)

library(terra)

source(here::here("R", "00-settings.R"))
output_filepath <- file.path(data_raw_dir, "glorys-bottom-temp.nc")
bt_june_yearly_nc_filepath <- file.path(data_gen_dir, "glorys-bottom-temp-june-yearly.nc")


# This script demonstrates how you can query Copernicus and get bottom temperature
# It DOES NOT need to be run to get "glorys-bottom-temp-june-yearly.rds" (which is on github)

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
cms_login()  # returns TRUE if authentication succeeded, FALSE otherwise

# --- 2. (Optional) confirm the variable name --------------
# GLORYS bottom temperature lives in the daily multi-year reanalysis product.
# Variable naming can shift between product revisions, so it's worth
# checking the product's variable list once before relying on "bottomT":
# vignette("product-info")  # or browse the product page directly


# --- 3. Download a bottom-temperature subset for a region -----------------
# Can load survey data to extract relevant parameters (time and spatial extent)
dog <- sdmTMB::dogfish
dog_sf <- gfplot::XY_to_sf(dog)

# Spatial extent ------
# Load data domain or manually specify region in this order: c(xmin, ymin, xmax, ymax)
region <- st_bbox(dog_sf)

# Temporal extent -----
# Get date range available for copernicus product
copernicus_product <- "GLOBAL_MULTIYEAR_PHY_001_030"
product_details <- cms_product_details(copernicus_product)
product_interval <- product_details$extent$temporal$interval[[1]]
product_interval

product_start <- as.Date(substr(product_interval[[1]], 1, 10)) # ISO 8601 timestamps -- the date is always the first 10 characters
product_end   <- as.Date(substr(product_interval[[2]], 1, 10))

# Restrict to what the product actually covers
# start date: 1993
# end date:   lags behind Sys.Date() by a processing delay.
start_date <- format(max(as.Date(sprintf("%d-01-01", min(dog$year))), product_start), "%Y-%m-%d")
end_date   <- format(min(Sys.Date(), product_end), "%Y-%m-%d")

# Download (via CopernicusMarine, which returns a stars object) and cache to
# disk as a proper multidimensional NetCDF (write_mdim() -- write_stars()
# only does 2D array).
# Safer to write to mdim and read back in using terra rather than converting to
# SpatRaster in memory because the `terra::rast()` is picky about dimension
# regularity/naming.
# Allows for handling as terra object for rest of workflow.
if (!file.exists(output_filepath)) {
  bottom_temp <- cms_download_subset(
    product       = "GLOBAL_MULTIYEAR_PHY_001_030",
    layer         = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
    variable      = "bottomT",
    region        = region,
    timerange     = c(start_date, end_date),
    progress      = TRUE
  )
  write_mdim(bottom_temp["bottomT"], output_filepath)
}

# --- 4. Save as terra and inspect -------
bt_daily <- terra::rast(output_filepath)

# Fail fast if GDAL/terra didn't parse the CF time/CRS metadata cleanly
if (anyNA(terra::time(bt_daily)) || !nzchar(terra::crs(bt_daily))) {
  cli::cli_abort("{.file {output_filepath}} did not parse valid time/CRS metadata into terra.")
}

bt_daily

plot(bt_daily[[1]], col = hcl.colors(100, "RdYlBu", rev = TRUE), axes = TRUE) # base R option

# Aggregate over time (e.g. mean), then plot:
bt_mean <- terra::app(bt_daily, mean, na.rm = TRUE)
plot(bt_mean, col = hcl.colors(100, "RdYlBu", rev = TRUE), axes = TRUE)

ggplot() +
  tidyterra::geom_spatraster(data = bt_daily[[1]]) +
  scale_fill_gradientn(colours = hcl.colors(100, "RdYlBu", rev = TRUE), na.value = NA)

# --- 5. Pull June values, aggregate to yearly means, and save -------------
june_idx <- which(format(terra::time(bt_daily), "%m") == "06")
bt_june <- bt_daily[[june_idx]]

years <- format(terra::time(bt_june), "%Y")
bt_june_yearly <- terra::tapp(bt_june, index = years, fun = mean, na.rm = TRUE)
names(bt_june_yearly) <- unique(years)
terra::time(bt_june_yearly) <- as.Date(paste0(unique(years), "-06-01"))

# Save derived/summarised data layer as a CF NetCDF (year as time step)
terra::writeCDF(
  bt_june_yearly, bt_june_yearly_nc_filepath,
  varname = "bottomT", longname = "Bottom temperature", unit = "degC",
  zname = "time", overwrite = TRUE
)

# --- 6. Inspect / plot the final product -----------------------------------
ggplot() +
  geom_sf(data = pacea::bc_coast) +
  tidyterra::geom_spatraster(data = bt_june_yearly) +
  facet_wrap(~ lyr, ncol = 4) +
  scale_fill_gradientn(colours = cmocean::cmocean("thermal")(100), na.value = NA) +
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  # scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  gfplot::coord_sf_auto(dog_sf) +
  theme(axis.text = element_blank())

# --- 7. Cite the product (required under CMEMS terms of use) --------------
# Only works if CopernicusMarine auth ran above
# cms_cite_product("GLOBAL_MULTIYEAR_PHY_001_030")$doi
# [1] "E.U. Copernicus Marine Service Information; Global Ocean Physics Reanalysis - GLOBAL_MULTIYEAR_PHY_001_030 (2012-08-26). DOI:10.48670/moi-00021"
