# Shift-indicator engine (make_sim_draws(), sim_stat(), summarise_trend(), ...)
# lives in the dshift package (./dshift). After editing it, reinstall with
# pak::local_install("dshift") or use devtools::load_all("dshift") interactively.
# @warning - temporarly location of library loading to capture dshift pkg since
#            we always call R/functions.R
library(dshift)

meep <- function(user = "jilliandunic", ...) {
  current_user <- Sys.info()['user']

  if (current_user == user) {
    beepr::beep(...)
  }
}

#' Shared internal borders between adjacent region polygons
#' #@todo - I am not sold this is even needed
#' Returns only the edges shared by two or more regions (the internal
#' dividers), so you can draw the boundaries between regions without the
#' surrounding outer box. Purely topological -- makes no assumption about
#' border orientation, so it finds perpendicular-to-coast splits and any
#' other adjacency equally.
#'
#' @param polys An `sf` of two or more (MULTI)POLYGON region features that
#'   meet edge-to-edge (adjacent regions must share exact boundaries, not
#'   overlap or leave gaps). CRS is preserved on the output.
#'
#' @return An `sf` of the shared borders as `LINESTRING`s, with an
#'   `n.overlaps` column (number of regions sharing each segment). Empty if
#'   no regions share an edge.
get_region_borders <- function(polys) {
  b <- sf::st_intersection(sf::st_boundary(polys))  # pairwise self-intersections
  b <- b[b$n.overlaps > 1, ]                        # keep only shared edges
  b <- sf::st_collection_extract(b, "LINESTRING") # not sure if we should suppress these yet until we test more polygons
  # b <- suppressWarnings(sf::st_collection_extract(b, "LINESTRING")) # option to suppress
  sf::st_sf(geometry = sf::st_geometry(b))          # drop misleading per-region attrs
}


#' Print a reminder of key model configuration
#'
#' #' #@todo - I am not sold this is even needed either - cod uses extensive caching instead
#' Surfaces the settings most likely to affect how shift indicators should
#' be interpreted -- model formula, spatial/spatiotemporal random field
#' structure, time-varying terms, and family -- as a `cli` message. Meant to be called
#' once near the top of a `.qmd` so this isn't easy to lose track of while
#' working through a long script.
#'
#' @param fit A fitted sdmTMB model object.
#'
#' @return `fit`, invisibly.
remind_model_info <- function(fit) {
  is_delta <- isTRUE(fit$family$delta)

  label_by_component <- function(x) {
    if (is_delta) paste0(x, " (", seq_along(x), ")", collapse = ", ") else x
  }

  time_varying <- if (is.null(fit$time_varying)) "none" else deparse(fit$time_varying)

  family_label <- if (!is.null(fit$family$clean_name)) {
    fit$family$clean_name
  } else {
    paste0(fit$family$family[1], "(link = '", fit$family$link[1], "')")
  }

  formula_label <- label_by_component(purrr::map_chr(fit$formula, deparse1))

  cli::cli_inform(c(
    "!" = "Model configuration reminder:",
    "*" = "Formula: {formula_label}",
    "*" = "Spatial: {label_by_component(fit$spatial)}",
    "*" = "Spatiotemporal: {label_by_component(fit$spatiotemporal)}",
    "*" = "Time-varying: {time_varying}",
    "*" = "Family: {family_label}"
  ))

  invisible(fit)
}

#' Check that saved inputs share one set of data sources
#'
#' `01-data-prep` attaches a `data_sources` attribute (survey, grid and
#' covariate files) to everything it saves, because a re-render with different
#' sources silently overwrites them. Call this after loading inputs downstream
#' so a fit and a grid built from different sources can't be combined, and so
#' the rendered document records the sources.
#'
#' @param ... Named objects carrying a `data_sources` attribute (e.g. data,
#'   prediction grid, fit).
#'
#' @return The shared `data_sources` list, invisibly.
check_data_sources <- function(...) {
  sources <- purrr::map(list(...), \(x) attr(x, "data_sources"))
  missing <- names(sources)[purrr::map_lgl(sources, is.null)]
  if (length(missing)) {
    cli::cli_abort("No {.field data_sources} attribute on {.val {missing}}; re-render {.file 01-data-prep.qmd}.")
  }
  # `created` differs between runs with identical sources
  same <- purrr::map_lgl(sources, \(x) identical(x[names(x) != "created"], sources[[1]][names(sources[[1]]) != "created"]))
  if (!all(same)) {
    cli::cli_abort(c(
      "Inputs were built from different data sources: {.val {names(sources)[!same]}} differ{?s/} from {.val {names(sources)[1]}}.",
      "i" = "Re-render {.file 01-data-prep.qmd} and refit, or load matching files."
    ))
  }
  src <- sources[[1]]
  cli::cli_inform(c(
    "i" = "Data sources (from 01-data-prep, {format(src$created, '%Y-%m-%d %H:%M')}):",
    purrr::set_names(
      paste0(names(src)[names(src) != "created"], ": ",
             purrr::map_chr(src[names(src) != "created"], \(x) if (is.null(x)) "NULL" else paste(x, collapse = ", "))),
      "*"
    )
  ))
  invisible(src)
}


# Covariate standard notes -------------------------------------------------------
# @note: IN PROGRESS - main testing is in NCAM - remember to port back to here
# A covariate is a NetCDF file (or the SpatRaster from terra::rast() on one):
# - CRS defined (any projection; `attach_covariate()` reprojects `x` to match)
# - one variable per raster
# - time in `terra::time()` (not layer names, which don't survive NetCDF):
#   the start of each layer's averaging window (e.g. the 1st of the month for
#   monthly means). A window ends where the next layer starts.
# - one layer per calendar year, or sub-annual layers (monthly, 5-day, ...);
#   a static covariate is a single layer with no time
#
# A source that isn't in the standard yet (e.g. no `terra::time()`) can be fixed
# in memory and passed as a SpatRaster; see `standardise_*()` in
# NCAM/analysis/spatial_model/R/covariate_helpers.R for examples.

#' Append a covariate raster to survey points or a prediction grid
#'
#' @param x An `sf` object with a CRS: survey points, or a polygon grid
#'   (values are extracted at each polygon's centroid; the polygons are kept).
#' @param covariate A NetCDF path or `SpatRaster` following the covariate
#'   standard above.
#' @param name Name of the column to append.
#' @param method Passed to `terra::extract()`: `"bilinear"` or `"simple"`
#'   (value of the cell underneath).
#' @param year_col Column of `x` holding the year, used to pick the matching
#'   layer of a yearly covariate. If `x` has no such column, `x` is repeated for
#'   every year in the covariate and `year_col` is added.
#' @param date_col Sub-annual covariates: column of `x` (`Date` or `POSIXct`)
#'   used to pick the layer whose window contains each date.
#' @param month Sub-annual covariates: month(s) to average within each year
#'   instead (e.g. for a prediction grid). The yearly means are then attached
#'   as for a yearly covariate, and a `month` column is added when `x` is
#'   expanded.
#'
#' @return `x` with `name` appended (and `year_col`, if `x` was expanded). The
#'   geometry and CRS of `x` are unchanged, so calls can be chained.
attach_covariate <- function(x, covariate, name, method = "bilinear", year_col = "year",
                             date_col = NULL, month = NULL) {
  if (!inherits(x, "sf")) {
    cli::cli_abort("{.arg x} must be an {.cls sf} object, not {.cls {class(x)}}.")
  }
  if (is.na(sf::st_crs(x))) {
    cli::cli_abort("{.arg x} must have a CRS.")
  }

  r <- if (is.character(covariate)) terra::rast(covariate) else covariate
  if (!inherits(r, "SpatRaster")) {
    cli::cli_abort("{.arg covariate} must be a {.cls SpatRaster} or a path to a NetCDF file.")
  }
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
  years <- if (static) NULL else if (is.numeric(time)) as.integer(time) else as.integer(format(time, "%Y"))
  sub_annual <- anyDuplicated(years) > 0
  if (!sub_annual && !is.null(month)) {
    cli::cli_warn("{.arg covariate} has no sub-annual layers; ignoring {.arg month}.")
  }

  if (sub_annual) {
    if (is.unsorted(time, strictly = TRUE)) {
      cli::cli_abort("{.arg covariate} layers must be in increasing {.fn terra::time}.")
    }
    # Window ends at the next layer's start; the last one, and any window before
    # a gap in the series, gets a typical step
    start <- as.numeric(as.POSIXct(time, tz = "UTC"))
    step <- diff(start)
    end <- start + c(pmin(step, 1.5 * stats::median(step)), stats::median(step))

    if (!is.null(month)) {
      # Average the chosen month(s) within each year; from here on it's yearly
      centre <- as.POSIXct((start + end) / 2, tz = "UTC")
      keep <- which(as.integer(format(centre, "%m")) %in% month)
      if (length(keep) == 0) {
        cli::cli_abort("{.arg covariate} has no layers in month{?s} {.val {month}}.")
      }
      yr <- format(centre[keep], "%Y")
      r <- terra::tapp(r[[keep]], index = yr, fun = mean, na.rm = TRUE)
      years <- as.integer(unique(yr))
      sub_annual <- FALSE
    } else if (is.null(date_col)) {
      cli::cli_abort("{.arg covariate} has sub-annual layers; supply {.arg date_col} (survey points) or {.arg month} (prediction grid).")
    } else if (!date_col %in% names(x)) {
      cli::cli_abort("{.arg x} has no column {.field {date_col}}.")
    }
  }

  pts <- sf::st_geometry(x)
  if (!all(sf::st_geometry_type(pts) == "POINT")) pts <- sf::st_centroid(pts)
  pts <- terra::vect(sf::st_transform(pts, sf::st_crs(terra::crs(r))))

  # Yearly covariate, no year column: repeat `x` for every year
  if (!static && !sub_annual && !year_col %in% names(x)) {
    gap_years <- setdiff(min(years):max(years), years)
    if (length(gap_years) > 0) {
      cli::cli_warn("{.arg covariate} has no layer for year{?s} {.val {gap_years}}.")
    }
    n <- nrow(x)
    x <- dplyr::slice(x, rep(seq_len(n), times = length(years)))
    x[[year_col]] <- rep(years, each = n)
    if (length(month) == 1) x$month <- month
    pts <- pts[rep(seq_len(n), times = length(years))]
  }

  # Layer of `r` for each row of `x`
  if (static) {
    layer <- rep(1L, nrow(x))
  } else if (sub_annual) {
    when <- as.numeric(as.POSIXct(x[[date_col]], tz = "UTC"))
    layer <- findInterval(when, start)
    layer[layer == 0] <- NA
    layer[when >= end[layer]] <- NA
    covers <- format(as.POSIXct(range(start, end), tz = "UTC"), "%Y-%m-%d")
  } else {
    layer <- match(x[[year_col]], years)
    covers <- as.character(range(years))
  }
  if (anyNA(layer)) {
    cli::cli_warn(
      "{sum(is.na(layer))} row{?s} of {.arg x} outside {.arg covariate} (which covers {.val {covers[1]}} to {.val {covers[2]}}); setting {.field {name}} to {.val {NA}}."
    )
  }

  # One extract per layer used: the full points x layers matrix can be GBs
  vals <- rep(NA_real_, nrow(x))
  for (rows in split(seq_len(nrow(x)), layer)) {
    vals[rows] <- terra::extract(r[[layer[rows[1]]]], pts[rows], method = method, ID = FALSE)[[1]]
  }
  x[[name]] <- vals

  n_na <- sum(is.na(x[[name]]))
  if (n_na > 0) {
    cli::cli_warn("{.field {name}} is {.val {NA}} for {n_na} of {nrow(x)} row{?s}.")
  }
  x
}
