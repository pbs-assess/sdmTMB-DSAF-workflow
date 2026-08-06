meep <- function(user = "jilliandunic", ...) {
  current_user <- Sys.info()['user']

  if (current_user == user) {
    beepr::beep(...)
  }
}

# combine_delta <- function(.data, fit) {
#   if (isTRUE(fit$family$delta)) {
#     fit$family[[1]]$linkinv(.data[, 1, ]) * fit$family[[2]]$linkinv(.data[, 2, ])
#   } else {
#     fit$family$linkinv(.data[, 1, ])
#   }
# }

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

# @TODO - document
#' median/95% interval/sd summary across simulation draws
#'
#' One row per `by` group. With a single `cols` column, output columns are
#' plain `est`/`lwr`/`upr`/`se`. With multiple `cols`, each gets its own
#' suffixed set (`est_x`, `est_y`, ...) via `across()`.
summarise_sim_ci <- function(df, cols, by) {
  col_names <- tidyselect::eval_select(rlang::enquo(cols), df)
  name_pattern <- if (length(col_names) == 1) "{.fn}" else "{.fn}_{.col}"

  df |>
    summarise(
      across({{ cols }}, list(
        est = ~ stats::median(.x),
        lwr = ~ stats::quantile(.x, 0.025),
        upr = ~ stats::quantile(.x, 0.975),
        se  = ~ stats::sd(.x)
      ), .names = name_pattern),
      .by = {{ by }}
    )
}

#' lm() coefficient summary across simulation draws
#'
#' Fits `value ~ centred_year` separately within each simulation draw and
#' each combination of `by` columns, then summarises the resulting
#' distribution of *each* coefficient (intercept and slope) across sims
#' (median/se/95% interval and probability of direction). Both coefficients
#' come from a single `lm()` fit per sim/group -- filter to the one you want
#' (`term == "centred_year"` for the trend slope used by every shift
#' indicator; `term == "(Intercept)"` for the fitted level at
#' `last_historical_year`, e.g. to sanity-check that regional intercepts
#' line up as expected) rather than calling this twice.
#'
#' @param df A data frame with one row per simulation draw (and whatever
#'   `by` groups apply), containing `sim` and `centred_year`.
#' @param value Column or expression to regress on `centred_year`, tidy-eval,
#'   e.g. `log(x)` or `qlogis(prop)`.
#' @param by <tidy-select> Grouping column(s) identifying one time series
#'   per simulation draw, e.g. `c(time_period, region)`. Accepts a
#'   character vector via `all_of()`, so a shared set of grouping columns
#'   can be defined once (e.g. at the top of the script) and reused or
#'   extended per indicator.
#'
#' @return A data frame with `term`, `by` columns, `estimate`, `se`, `lwr`,
#'   `upr`, `prob_direction` -- one row per term (`"(Intercept)"` and
#'   `"centred_year"`) per `by` group.
summarise_sim_trend <- function(df, value, by) {
  coefs <- df |>
    mutate(.value = {{ value }}) |>
    reframe(
      {
        co <- coef(lm(.value ~ centred_year))
        tibble::tibble(term = names(co), coef_est = as.numeric(co))
      },
      .by = c(sim, {{ by }})
    )

  coefs |>
    summarise(
      estimate       = stats::median(coef_est),
      se             = stats::sd(coef_est),
      lwr            = stats::quantile(coef_est, 0.025),
      upr            = stats::quantile(coef_est, 0.975),
      prob_direction = 2 * min(mean(coef_est > 0), mean(coef_est < 0)),
      .by = c(term, {{ by }})
    )
}


#@TODO - clean up docs
#' position where cumulative abundance crosses each target quantile.
#'
#' @param df A data frame with one row per cell x sim (and `by` groups),
#'   containing `log_abundance` and the column referenced by `rank_var`.
#' @param rank_var Tidy-eval column to rank/interpolate along, e.g. `Y`, `X`
#' @param quantiles Named numeric vector of target cumulative-abundance
#'   quantiles, e.g. `c(trailing = 0.05, leading = 0.95)`.
#' @param by <tidy-select> Grouping column(s) identifying one ranking series
#'   per sim, e.g. `c(time_period, region, year, sim)`.
#'
#' @return A data frame with `by` columns, `edge`, and `{rank_var}_edge`
#'   (interpolated value of `rank_var` at each quantile).
get_range_edge_sim <- function(df, rank_var, quantiles, by) {
  edge_col <- paste0(rlang::as_label(rlang::enquo(rank_var)), "_edge")

  df |>
    mutate(.rank_var = {{ rank_var }}) |>
    group_by(pick({{ by }})) |>
    arrange(.rank_var, .by_group = TRUE) |>
    mutate(cumul_prop = cumsum(exp(log_abundance)) / sum(exp(log_abundance))) |> #L.201 - 203 in get-range-edge.R sdmTMB
    reframe(
      edge = names(quantiles),
      "{edge_col}" := approx(cumul_prop, .rank_var, xout = quantiles, rule = 2, ties = "ordered")$y
    )
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
  if (is.null(polys)) return(NULL)
  b <- sf::st_intersection(sf::st_boundary(polys))  # pairwise self-intersections
  b <- b[b$n.overlaps > 1, ]                        # keep only shared edges
  b <- sf::st_collection_extract(b, "LINESTRING") # not sure if we should suppress these yet until we test more polygons
  # b <- suppressWarnings(sf::st_collection_extract(b, "LINESTRING")) # option to suppress
  sf::st_sf(geometry = sf::st_geometry(b))          # drop misleading per-region attrs
}
