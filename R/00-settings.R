# Directory structure -------
# Define directories in once place to make it easier to change downstream pointers
data_raw_dir <- here::here("data-raw")
dir.create(data_raw_dir, showWarnings = FALSE, recursive = TRUE)

data_gen_dir <- here::here("data-generated")
dir.create(data_gen_dir, showWarnings = FALSE, recursive = TRUE)

model_input_dir <- file.path(data_gen_dir, "data-for-sdmTMB")
dir.create(model_input_dir, showWarnings = FALSE, recursive = TRUE)

# Note - in the end this isn't something we will cache for all species; but it
# is handy to save them while we are testing.
model_out_dir <- file.path(data_gen_dir, "model-outputs")
dir.create(model_out_dir, showWarnings = FALSE, recursive = TRUE)

fit_dir <- file.path(model_out_dir, "fits")
dir.create(fit_dir, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Other settings -------
# Set utm zone in one place so you are consistent across your data/meshes/models/predictions
# one way to choose this is:
# sdmTMB:::get_crs(<survey-set-dataframe>, ll_names = c("lon_start", "lat_start"))
utm_crs <- 32609

# -----------------------------------------------------------------------------
# plotting

theme_set(ggsidekick::theme_sleek())

# use bayesplot blues
blues <- c("lightest" = "#d1e1ec", "light" = "#b3cde0", "mid" = "#6497b1", "#005b96", "dark" = "#03396c", "darkest" = "#011f4b")