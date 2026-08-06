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

fit_dir <- file.path(model_out_dir, "fits")
dir.create(fit_dir, showWarnings = FALSE)


ggplot2::theme_set(gfplot::theme_pbs())

blues <- bayesplot::color_scheme_get("blue")

# @Check: this is the wording to use?
pvalue_labels <- dplyr::tribble(
  ~p.value, ~p.significant,
  0.0005,   "Very strong",
  0.005,    "Strong",
  0.03,     "Moderate",
  0.5,      "Not significant"
)

# @TODO: PULL functions intow function script later

# TODO: move to gfplot
# X/Y assumes sdmTMB::add_utm_columns() output.
# @TODO: add crs to documentation of `gfdata::load_survey_blocks()`
# @TODO: allow 'xy' OR 'XY' in `gfdata::load_survey_blocks()`
# @TODO: allow filtering of the different grids upon loading `gfdata::load_survey_blocks()`

#' Convert XY coordinates to sf object
#'
#' Converts `X`/`Y` coordinate columns, as produced by
#' `sdmTMB::add_utm_columns()` (typically in km), to an sf object. Also
#' works for general conversions of point data to sf.
#'
#' @param data Data frame containing coordinate columns.
#' @param coords Vector of coordinate column names (default: c("X", "Y")).
#' @param mult Multiplier applied to `coords` before constructing geometry
#'   (default: 1000, i.e. km -> m). Automatically set to 1 if `crs_from` is
#'   a geographic (lon/lat) CRS.
#' @param crs_from Source coordinate reference system. Defaults to
#'   EPSG:32609 (UTM zone 9N), matching Pacific region sdmTMB workflows -- pass
#'   explicitly for other regions.
#' @param crs_to Target coordinate reference system (default: 4326, WGS84).
#'
#' @return sf object with geometry in `crs_to`.
#' @export
XY_to_sf <- function(data, coords = c("X", "Y"),
                     mult = 1000,
                     crs_from = 32609, crs_to = 4326) {
  if (!all(coords %in% names(data))) {
    missing_cols <- coords[!coords %in% names(data)]
    stop("Coordinate column(s) not found: ", paste(missing_cols, collapse = ", "))
  }

  if (missing(crs_from)) {
    message(
      "Assuming crs_from = 32609 (UTM zone 9N). ",
      "Pass crs_from explicitly for other regions."
    )
  }

  if (sf::st_is_longlat(sf::st_crs(crs_from))) mult <- 1

  data |>
    dplyr::mutate(
      x = .data[[coords[1]]] * mult,
      y = .data[[coords[2]]] * mult
    ) |>
    sf::st_as_sf(coords = c("x", "y"), crs = crs_from) |>
    sf::st_transform(crs = crs_to)
}
