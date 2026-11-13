# Shift-indicator engine (make_sim_draws(), sim_stat(), summarise_trend(), ...)
# and attach_covariate() (with the covariate standard in ?attach_covariate)
# live in the dshift package (./dshift). Loaded from source on every run so the
# pipeline always uses the dshift in this checkout (no reinstall after pulls or
# edits). export_all = FALSE mimics an installed package (exports only).
pkgload::load_all(here::here("dshift"), export_all = FALSE, quiet = TRUE)

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

