#' Append a covariate raster to survey points or a prediction grid
#'
#' Extracts one covariate at each row of `x`, choosing the raster layer that
#' matches the row's year or date.
#'
#' @section Covariate standard:
#' A covariate is a NetCDF file, or the `SpatRaster` from `terra::rast()` on
#' one, with:
#' - a defined CRS (any projection; `x` is reprojected to match)
#' - one variable
#' - time in `terra::time()`, not in layer names, which don't survive NetCDF
#' - one of three layouts:
#'   - **static**: a single layer with no time (e.g. depth)
#'   - **yearly**: one layer per calendar year
#'   - **monthly**: one layer per month, dated within that month (e.g. the 1st)
#'
#' A source that doesn't follow the standard can be fixed in memory (e.g.
#' setting `terra::time()`) and passed as a `SpatRaster`. Finer cadences
#' (daily, 5-day) aren't supported; aggregate them to monthly first.
#'
#' @param x An `sf` object with a CRS: survey points, or a polygon grid (values
#'   are extracted at each polygon's centroid; the polygons are kept). For a
#'   time-varying covariate, `x` must already have one row per location and
#'   year (expand a grid over years before calling, e.g. with
#'   `tidyr::expand_grid()` followed by `sf::st_as_sf()`).
#' @param covariate A NetCDF path or `SpatRaster` following the covariate
#'   standard.
#' @param name Name of the column to append.
#' @param method Passed to `terra::extract()`: `"bilinear"` or `"simple"`
#'   (value of the cell underneath).
#' @param year_col Column of `x` holding the year. Used for yearly covariates,
#'   and for monthly covariates averaged over `month`.
#' @param date_col Monthly covariates: column of `x` (`Date` or `POSIXct`); each
#'   row gets the layer for its own year and month (e.g. survey sets).
#' @param month Monthly covariates: month(s) to average within each year instead
#'   of `date_col` (e.g. a prediction grid with no dates). Each row then gets
#'   the mean for its `year_col`.
#'
#' @return `x` with column `name` appended. Rows, geometry and CRS are
#'   unchanged, so calls can be chained.
#' @export
attach_covariate <- function(x, covariate, name, method = "bilinear",
                             year_col = "year", date_col = NULL, month = NULL) {
  rlang::check_installed(c("sf", "terra"), reason = "to attach covariates.")
  if (!inherits(x, "sf") || is.na(sf::st_crs(x))) {
    cli::cli_abort("{.arg x} must be an {.cls sf} object with a CRS.")
  }
  r <- if (is.character(covariate)) terra::rast(covariate) else covariate
  if (!nzchar(terra::crs(r))) {
    cli::cli_abort("{.arg covariate} must have a CRS.")
  }

  time <- terra::time(r)
  static <- all(is.na(time))
  if (static && terra::nlyr(r) > 1) {
    cli::cli_abort("{.arg covariate} has {terra::nlyr(r)} layers but no {.fn terra::time}.")
  }
  if (!static && anyNA(time)) {
    cli::cli_abort("{.arg covariate} has layers with missing {.fn terra::time}.")
  }
  layer_year <- if (static) NA_integer_ else if (is.numeric(time)) as.integer(time) else as.integer(format(time, "%Y"))
  monthly <- !static && anyDuplicated(layer_year) > 0

  if (monthly) {
    if (is.numeric(time)) {
      cli::cli_abort("{.arg covariate} has repeated numeric {.fn terra::time}; monthly layers need dates.")
    }
    layer_month <- as.integer(format(time, "%m"))
    if (anyDuplicated(paste(layer_year, layer_month)) > 0) {
      cli::cli_abort("{.arg covariate} has more than one layer per month; aggregate it to monthly first.")
    }
  } else if (!is.null(month)) {
    cli::cli_warn("{.arg covariate} has no monthly layers; ignoring {.arg month}.")
  }

  # Monthly covariate averaged over `month`: collapse to one layer per year
  if (monthly && !is.null(month)) {
    keep <- layer_month %in% month
    if (!any(keep)) {
      cli::cli_abort("{.arg covariate} has no layers in month{?s} {.val {month}}.")
    }
    years <- sort(unique(layer_year[keep]))
    r <- terra::rast(purrr::map(years, \(y) terra::mean(r[[which(keep & layer_year == y)]], na.rm = TRUE)))
    layer_year <- years
    monthly <- FALSE
  }

  # Layer of `r` for each row of `x`
  if (static) {
    layer <- rep(1L, nrow(x))
  } else if (monthly) {
    if (is.null(date_col) || !date_col %in% names(x)) {
      cli::cli_abort(c(
        "{.arg covariate} has monthly layers; supply {.arg date_col} (a date column of {.arg x}) or {.arg month}.",
        i = "Use {.arg date_col} for survey points and {.arg month} for a prediction grid."
      ))
    }
    layer <- match(format(as.POSIXct(x[[date_col]], tz = "UTC"), "%Y-%m"), sprintf("%d-%02d", layer_year, layer_month))
  } else {
    if (!year_col %in% names(x)) {
      cli::cli_abort(c(
        "{.arg x} has no column {.field {year_col}}.",
        i = "For a prediction grid, expand it over years first, e.g. {.code tidyr::expand_grid(grid, year = years) |> sf::st_as_sf()}."
      ))
    }
    layer <- match(x[[year_col]], layer_year)
  }
  if (anyNA(layer)) {
    cli::cli_warn(
      "{sum(is.na(layer))} row{?s} of {.arg x} outside the years covered by {.arg covariate} ({min(layer_year)}-{max(layer_year)}); setting {.field {name}} to {.val {NA}}."
    )
  }

  pts <- sf::st_geometry(x)
  if (!all(sf::st_geometry_type(pts) == "POINT")) pts <- sf::st_centroid(pts)
  pts <- terra::vect(sf::st_transform(pts, terra::crs(r)))

  # One extract per layer used: the full points x layers matrix can be GBs
  rows_by_layer <- split(seq_len(nrow(x)), layer)
  vals <- rep(NA_real_, nrow(x))
  vals[unlist(rows_by_layer)] <- unlist(purrr::imap(rows_by_layer, \(rows, l) {
    terra::extract(r[[as.integer(l)]], pts[rows], method = method, ID = FALSE)[[1]]
  }))
  x[[name]] <- vals

  n_na <- sum(is.na(vals))
  if (n_na > 0) {
    cli::cli_warn("{.field {name}} is {.val {NA}} for {n_na} of {nrow(x)} row{?s}.")
  }
  x
}
