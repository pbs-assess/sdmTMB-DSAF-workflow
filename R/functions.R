meep <- function(user = "jilliandunic", ...) {
  current_user <- Sys.info()['user']

  if (current_user == user) {
    beepr::beep(...)
  }
}

#' Combine `project()` linear predictors into density or log density
#'
#' Applies each linear predictor's own `linkinv()` before combining, so this
#' works for single-linear-predictor models, standard delta models (e.g.
#' logit + log), and poisson-link delta models alike -- no assumption about
#' which links were used. For poisson-link families the log-density case
#' reduces to `est1 + est2` because both components use a log link; for
#' standard delta models it does not, so the two forms aren't interchangeable.
#'
#' Always sums in log space first, then exponentiates only if `log = FALSE`.
#' For poisson-link delta models both components are unbounded `exp()`
#' outputs, so their product can overflow to `Inf` even when each component
#' is individually representable -- computing the product before logging
#' would lose an otherwise-finite log-density to that overflow.
#'
#' @param projections Output of `project()`: an `est1`/`est2` list for delta
#'   models, or an object with `est1` for single-linear-predictor models.
#' @param fit A fitted sdmTMB model object.
#' @param log If `TRUE` (default), return log density. If `FALSE`, return
#'   density on the natural scale (`exp()` of the log density) -- this can
#'   legitimately be `Inf` for poisson-link models if the log density itself
#'   overflows double precision once exponentiated.
#'
#' @return A matrix of log-density (or density) draws, one row per `newdata`
#'   row and one column per `nsim` draw.
return_density <- function(projections, fit, log = TRUE) {
  if (isTRUE(fit$family$delta)) {
    log_density <- log(fit$family[[1]]$linkinv(projections$est1)) +
      log(fit$family[[2]]$linkinv(projections$est2))
  } else {
    log_density <- log(fit$family$linkinv(projections$est1))
  }
  if (log) log_density else exp(log_density)
}

#' Print a reminder of key model configuration
#'
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

#' Auto-derive a model name from its configuration
#'
#' Builds a filesystem-safe slug encoding the modeling axes this project
#' actually varies (family, spatiotemporal RF, time-varying), so a model's
#' filename records what it actually is instead of a hand-picked mnemonic
#' that can drift out of sync with the fit (e.g. implying AR1 for a model
#' that's actually RW). Spatiotemporal and time-varying get separate
#' `st-`/`tv-` tags since both can independently be RW or AR1 -- and the
#' latter's type isn't even a top-level field on the fit object (it's buried
#' in `fit$tmb_data$ar1_time`; even `print.sdmTMB()` doesn't label it) --
#' conflating the two is the exact ambiguity this replaces.
#'
#' This name is meant to be the single record of model config -- used as the
#' filename, and carried into `sims_re`/`sims_fe` as the `model_name` column
#' -- rather than also keeping a separate descriptive label in sync.
#'
#' Scoped deliberately to axes this project currently varies. Add a tag here
#' if/when a genuinely new axis gets used (e.g. threshold effects) -- don't
#' pre-build support for structure nothing here uses yet.
#'
#' @param fit A fitted sdmTMB model object.
#' @param label A short human label for the species/dataset -- the one piece
#'   this can't recover from `fit` itself, e.g. `"dogfish"`.
#'
#' @return A length-1 character string safe for use as a filename (no
#'   extension).
build_model_name <- function(fit, label) {
  family <- if (!is.null(fit$family$clean_name)) {
    sub("\\(.*", "", fit$family$clean_name)
  } else {
    fit$family$family[1]
  }

  st_tag <- paste(unique(unlist(fit$spatiotemporal)), collapse = "-")

  tv_tag <- if (is.null(fit$time_varying)) {
    NULL
  } else if (isTRUE(fit$tmb_data$ar1_time == 1)) "tv-ar1" else "tv-rw"

  paste(c(label, family, paste0("st-", st_tag), tv_tag), collapse = "-")
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
  rank_name <- rlang::as_name(rlang::ensym(rank_var))
  by_names <- names(dplyr::select(df, {{ by }}))
  edge_col <- paste0(rank_name, "_edge")

  dt <- data.table::as.data.table(
    df[, c(by_names, rank_name, "log_abundance"), drop = FALSE]
  )

  # Sorting once globally is substantially cheaper than arranging within each
  # dplyr group. The cumulative abundance calculation still occurs per group.
  data.table::setorderv(dt, c(by_names, rank_name))
  dt[, cumul_prop := cumsum(exp(log_abundance)) / sum(exp(log_abundance)),
     by = by_names]

  out <- dt[, {
    edge_values <- approx(
      cumul_prop,
      get(rank_name),
      xout = quantiles,
      rule = 2,
      ties = "ordered"
    )$y

    data.table::data.table(
      edge = names(quantiles),
      edge_value = edge_values
    )
  }, by = by_names]

  data.table::setnames(out, "edge_value", edge_col)
  tibble::as_tibble(out)
}

#@TODO review docs
#@TODO use data.table because these matrices can be kind of big -- see:
#  https://github.com/pbs-assess/dogfish-trends/blob/9e472c7e5202cd232d2035038b41f28df0ddfba0/analysis/016-calculate-weighted-mean-depth.R#L55
#@optimize for larger datasets working with lots of draws this could be slow

#' Tidy a `project()` draws matrix into a long data frame
#'
#' Takes a matrix of log-density draws (e.g. `projections_re` or
#' `projections_fe` after `return_density()`), binds it onto `prediction_grid`,
#' and pivots to one row per cell x year x draw -- the shared shape every
#' shift-indicator function downstream (`summarise_sim_trend()`,
#' `get_range_edge_sim()`, etc.) expects. Lets different `project()` scenarios
#' be pushed through the same indicator pipeline and compared.
#'
#' @param prediction_grid The `newdata` grid passed to `project()` -- one row
#'   per cell x region x year, must contain `year`, `area`, and `region`.
#' @param projections Matrix of log-density draws (e.g. `return_density()`
#'   output) with rows matching `prediction_grid` and one column per draw
#'   (`project()`'s `nsim`).
#' @param last_historical_year Last year with survey data; splits
#'   `time_period` into `"historical"`/`"forecast"`.
#'
#' @note `projections` should come from expected-density draws
#'   (`observation_error = FALSE`, the `project()` default) -- no
#'   observation-level noise is added, so this is not the same object
#'   `simulate(..., observation_error = TRUE)` would give you.
#'
#'   What varies across `draw` depends entirely on the `project()` call that
#'   produced `projections`, not on this function:
#'   * `uncertainty = "random"` -- theta (fixed effects) and any `sim_re`-on
#'     random effects are sampled from their approximate joint distribution
#'     (`mle-mvn`); these are genuine posterior-predictive-style draws.
#'   * `uncertainty = "none"` -- theta is fixed at the MLE in every draw.
#'     Historical rows are then bit-for-bit identical across all `draw`
#'     values; only forecast-year rows vary, and only if `sim_re` turns on
#'     forward simulation of the AR1 process there. Don't be surprised when a
#'     CI computed on the historical portion of such a run comes out zero
#'     width -- that is correct, not a bug.
#'
#' @return A data frame shaped like `prediction_grid`, with added `draw`,
#'   `time_period`, `log_density`, `log_abundance` (density x cell area, in
#'   logspace like `get_index_sims()`) columns. Doubled with a
#'   `region == "all"` copy alongside the named regions (or renamed to
#'   `"all"` outright if `prediction_grid` only has one region).
make_sims_df <- function(prediction_grid, projections, last_historical_year) {
  sims <- prediction_grid |>
    mutate(time_period = ifelse(year <= last_historical_year, "historical", "forecast")) |>
    bind_cols(
      as_tibble(projections, .name_repair = ~ paste0("sim", seq_along(.x)))
    ) |>
    tidyr::pivot_longer(
      cols = starts_with("sim"),
      names_to = "sim",
      names_pattern = "^sim(\\d+)",
      names_transform = list(sim = as.integer),
      values_to = "log_density" # assuming offset = 1km2 (log(offset) = 0 used in prediction)
    ) |>
    mutate(
      time_period = factor(time_period, levels = c("historical", "forecast")),
      log_abundance = log_density + log(area) # add area in logspace like get_index_sims()
    )

  if (dplyr::n_distinct(sims$region) > 1) {
    bind_rows(sims, sims |> mutate(region = "all")) # full domain and regions processed in single step
  } else {
    sims |> mutate(region = "all") # if only a single region provided
  }
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
