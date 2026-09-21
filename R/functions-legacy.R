# Long-format (tidy, one row per cell x sim) shift-indicator functions that
# preceded the matrix-native engine in functions.R. Kept only so
# verify-shift-indicators-matrix.R can check the two paths agree; not used
# by the pipeline.

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
