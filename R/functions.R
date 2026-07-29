meep <- function(user = "jilliandunic", ...) {
  current_user <- Sys.info()['user']

  if (current_user == user) {
    beepr::beep(...)
  }
}


#@TODO edit/review documentation

#' Identify which random field columns are active in a fitted model
#'
#' Spatial and spatiotemporal random fields can be switched on/off per
#' delta component. This checks `fit$spatial`/`fit$spatiotemporal` and
#' returns only the column name(s) for the component(s) that are actually
#' turned on -- the two delta components are never summed for plotting.
#'
#' @param fit A fitted sdmTMB model object.
#' @param prefix Column prefix identifying the field type: `"omega_s"` for
#'   the spatial field, or the spatiotemporal field's prefix (e.g.
#'   `"epsilon_st"`). For delta models these are suffixed with `1`/`2` to
#'   form the actual column name(s), e.g. `"omega_s1"`.
#'
#' @return A character vector of column name(s) to plot, or `NULL`
#'   (invisibly, with a `cli` info message) if the field is off in every
#'   component.
field_columns <- function(fit, prefix) {
  field_name <- if (prefix == "omega_s") "spatial" else "spatiotemporal"
  status <- fit[[field_name]]
  is_delta <- isTRUE(fit$family$delta)
  on <- if (is_delta) which(status != "off") else if (any(status != "off")) 1L else integer(0)
  if (length(on) == 0) {
    cli::cli_alert_info("{.strong {field_name}} random field not used in the model.")
    return(invisible(NULL))
  }
  if (!is_delta) return(prefix)
  cli::cli_alert_info("{.strong {field_name}} random field on for delta component{?s} {on}.")
  paste0(prefix, on)
}

#' Plot a spatial or spatiotemporal random field
#'
#' Maps the random field(s) identified by [field_columns()] over the
#' prediction grid, faceted by delta component (if any) and year.
#'
#' @param pred_poly An `sf` object of prediction grid polygons containing
#'   the field column(s) named by `prefix`/[field_columns()].
#' @param fit A fitted sdmTMB model object.
#' @param prefix Column prefix identifying the field type; see
#'   [field_columns()].
#' @param title Plot title.
#'
#' @return A ggplot object, or `NULL` (invisibly) if the field is not used
#'   in the model.
plot_field <- function(pred_poly, fit, prefix, title) {
  cols <- field_columns(fit, prefix)
  if (is.null(cols)) return(invisible(NULL))

  pred_poly |>
    tidyr::pivot_longer(cols = all_of(cols), names_to = "component", values_to = "value") |>
    ggplot() +
    geom_sf(aes(fill = value), colour = NA) +
    scale_fill_gradient2() +
    facet_grid(component ~ year) +
    ggtitle(title)
}


#@TODO: TEST ME with more data
#' Shared internal borders between adjacent region polygons
#'
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
