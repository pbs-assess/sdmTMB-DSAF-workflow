# Define directories in once place to make it easier to change downstream pointers
data_raw_dir <- here::here("data-raw")

data_gen_dir <- here::here("data-generated")
dir.create(data_gen_dir, showWarnings = FALSE, recursive = TRUE)

model_input_dir <- file.path(data_gen_dir, "data-for-sdmTMB")
dir.create(model_input_dir, showWarnings = FALSE, recursive = TRUE)

# Note - in the end this isn't something we will cache for all species; but it
# is handy to save them while we are testing.
model_out_dir <- file.path(data_gen_dir, "model-outputs")
dir.create(model_out_dir, showWarnings = FALSE, recursive = TRUE)

ggplot2::theme_set(gfplot::theme_pbs())

# @Check: this is the wording to use?
pvalue_labels <- tribble(
  ~p.value, ~p.significant,
  0.0005,   "Very strong",
  0.005,    "Strong",
  0.03,     "Moderate",
  0.5,      "Not significant"
)

# @TODO: PULL functions intow function script later
# @TODO: some of these might move to gfplot, so would rather do that than pull them into something local
meep <- function(user = "jilliandunic", ...) {
  current_user <- Sys.info()['user']

  if (current_user == user) {
    beepr::beep(...)
  }
}

# NEED TO PUT THIS IN gfdata
# BUT -- need to think more carefully if this might be used across the country,
# maybe the defaults don't make as much sense and we might need more defensive
# programming around it and our input grids.
# Also other notes:
# @TODO: add crs to documentation of `gfdata::load_survey_blocks()`
# @TODO: allow 'xy' OR 'XY' in `gfdata::load_survey_blocks()`
# @TODO: allow filtering of the different grids upon loading `gfdata::load_survey_blocks()`

#' Convert coordinates to sf object
#'
#' Converts sdmTMB XY coordinates(km) to sf object. Also works for general
#' conversions of point data to coordinates.
#'
#' @param data Data frame containing coordinate columns
#' @param coords Vector of coordinate column names (default: c("X", "Y"))
#' @param mult Multiplier for coordinates (default: 1000, converts km to m).
#'   Automatically set to 1 if crs_from = 4326.
#' @param crs_from Source coordinate reference system (default: 32609)
#' @param crs_to Target coordinate reference system (default: 4326)
#'
#' @return sf object, defaults to WGS84 (EPSG:4326)
#' @export
XY_to_sf <- function(data, coords = c("X", "Y"),
                     mult = 1000,
                     crs_from = 32609, crs_to = 4326) {
  if (!all(coords %in% names(data))) {
    missing <- coords[!coords %in% names(data)]
    stop("Coordinate column(s) not found: ", paste(missing, collapse = ", "))
  }

  if (crs_from == 4326) mult <- 1

  df <- data |>
    dplyr::mutate(
      x = .data[[coords[1]]] * mult,
      y = .data[[coords[2]]] * mult
    )

  df <- df |>
    sf::st_as_sf(coords = c("x", "y"), crs = crs_from) |>
    sf::st_transform(crs = crs_to)
  df
}