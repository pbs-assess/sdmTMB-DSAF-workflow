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


# Open questions for make_sim_draws() / the sim_draws object:
# - model-agnostic #@todo - write test to show model agnostic
# - should last_historical year be an argument? optional?
# - are optional labels messy?
# - should we use a class? does `sim_draws` make sense for a name?

#' Build list of projection output with associated metadata
#'
#' Takes a matrix of projected density values and combines it with metadata
#' from the prediction grid into a single list, so downstream functions can
#' do fast matrix operations on `density`/`abundance` and look up metadata
#' (grouping columns, area) via `grid` when needed.
#'
#' @param prediction_grid The `newdata` grid used to generate `log_density`
#'   (row order must match `log_density`'s rows) -- one row per cell x
#'   region x year, must contain `year` and `area`.
#' @param log_density Log-scale density matrix, one row per `prediction_grid`
#'   row and one column per simulation draw. For sdmTMB, this is
#'   `return_density(projections, fit, log = TRUE)` (the default) -- kept as
#'   the canonical scale here since it's what `return_density()` computes
#'   natively, and since it's cheaper to keep only one scale on a large grid
#'   than to carry a duplicate natural-scale matrix nobody may need;
#'   downstream functions `exp()` what they need on the fly.
#' @param last_historical_year Last year with survey data; splits
#'   `time_period` into `"historical"`/`"forecast"`.
#' @param type,model_name Optional labels carried through for bookkeeping
#'   (e.g. `"re"`/`"fe"`, a model name), not used internally.
#'
#' @return A `sim_draws` object: a list with `grid` (`prediction_grid` plus
#'   `time_period`), `log_density`, `log_abundance` (`log_density +
#'   log(grid$area)`), `type`, `model_name`.
#'
#' @export
make_sim_draws <- function(prediction_grid, log_density, last_historical_year,
                            type = NULL, model_name = NULL) {
  if (nrow(log_density) != nrow(prediction_grid)) {
    cli::cli_abort("{.arg log_density} must have one row per row of {.arg prediction_grid}.")
  }
  grid <- prediction_grid
  grid$time_period <- factor(
    ifelse(grid$year <= last_historical_year, "historical", "forecast"),
    levels = c("historical", "forecast")
  )
  structure(
    list(
      grid = grid,
      log_density = log_density,
      log_abundance = log_density + log(grid$area),
      type = type,
      model_name = model_name
    ),
    class = "sim_draws"
  )
}

#' Keyed-matrix container for simulation-draw summaries
#'
#' `list(key, mats)` -- row `i` of `key` describes row `i` of every matrix
#' in `mats`. Lets [summarise_estimates()] and [summarise_trend()] operate
#' directly on `groups x nsim` matrices.
#'
#' @param key A data frame, one row per group.
#' @param mats Named list of `nrow(key) x nsim` matrices.
#'
#' @return A `sim_stat` object: `list(key, mats)`.
#'
#' @keywords internal
#'
sim_stat <- function(key, mats) {
  key <- as.data.frame(key)
  rownames(key) <- NULL
  nr <- vapply(mats, nrow, integer(1))
  if (!all(nr == nrow(key))) {
    cli::cli_abort("Every matrix in {.arg mats} must have {nrow(key)} rows to match {.arg key}.")
  }
  nc <- vapply(mats, ncol, integer(1))
  if (length(unique(nc)) > 1) {
    cli::cli_abort("Every matrix in {.arg mats} must share the same number of columns (simulation draws).")
  }
  mats <- lapply(mats, function(m) { dimnames(m) <- NULL; m })
  structure(list(key = tibble::as_tibble(key), mats = mats), class = "sim_stat")
}

#' Melt a `sim_stat` to a long tibble (one row per group x simulation draw)
#'
#' Converts key-matrix `sim_stat` to a long tibble. Useful for inspection,
#' debugging, and the verify script's per-sim comparisons; the summary
#' functions below do not use it.
#'
#' @param x A `sim_stat` object.
#'
#' @return A tibble: `x$key` repeated once per simulation draw, plus `sim`
#'   and one column per entry of `x$mats`.
#'
#' @export
sim_stat_tibble <- function(x) {
  G <- nrow(x$key); nsim <- ncol(x$mats[[1]])
  out <- x$key[rep(seq_len(G), times = nsim), , drop = FALSE]
  out$sim <- rep(seq_len(nsim), each = G)
  for (nm in names(x$mats)) out[[nm]] <- as.vector(x$mats[[nm]])
  tibble::as_tibble(out)
}

#' Pull individual simulation draws (plus aggregates) back onto the grid
#'
#' For inspecting a raw `nrow(grid) x nsim` matrix spatially -- e.g. is a
#' `project()` `epsilon_st` "hot spot" a real feature of the field, or an
#' artifact of one noisy draw? Deliberately takes a plain matrix rather than
#' a `sim_draws`/`sim_stat` object: those bake in density semantics
#' (`log_abundance = log_density + log(area)`) that don't apply to an
#' arbitrary field like a raw spatiotemporal random effect, so this works
#' directly off whatever `project()` (or `sim_draws$log_density`, etc.)
#' hands back, no conversion required.
#'
#' @param grid A data frame, one row per row of `M`.
#' @param M A `nrow(grid) x nsim` matrix, or a `nrow(grid) x k x nsim` array
#'   (e.g. `project(..., sims_var = "proj_epsilon_st_A_vec")` on a delta
#'   model, `k` = linear predictor) -- use `component` to pick a slice of the
#'   latter.
#' @param sims Integer vector of individual draw indices to pull; `NULL` to
#'   skip.
#' @param agg Named list of row-wise aggregate functions to also compute
#'   across all draws (e.g. `rowMeans`, `matrixStats::rowMedians`); `NULL`
#'   to skip.
#' @param component For a 3-D `M`, which index along its 2nd dimension to use
#'   (e.g. `1`/`2` for a delta model's two linear predictors). Required when
#'   `M` is 3-D; ignored otherwise.
#' @param value_name Name for the output value column.
#'
#' @return `grid`, repeated once per requested draw/aggregate, plus `draw`
#'   (`"sim_<i>"` or the `agg` name) and `value_name`.
#'
#' @export
sim_draw_slices <- function(grid, M, sims = 1,
                             agg = list(mean = rowMeans, median = matrixStats::rowMedians),
                             component = NULL, value_name = "value") {
  if (length(dim(M)) == 3) {
    if (is.null(component)) {
      cli::cli_abort("{.arg M} is 3-D (e.g. a multi-linear-predictor {.arg sims_var} array) -- specify {.arg component} to pick a slice.")
    }
    M <- M[, component, ]
  }
  if (nrow(grid) != nrow(M)) {
    cli::cli_abort("{.arg grid} and {.arg M} must have the same number of rows.")
  }
  pull <- function(v, label) { out <- grid; out$draw <- label; out[[value_name]] <- v; out }
  draws <- purrr::map(sims, \(s) pull(M[, s], paste0("sim_", s)))
  aggs  <- purrr::imap(agg, \(f, nm) pull(f(M), nm))
  dplyr::bind_rows(draws, aggs)
}

#' Stack `sim_stat` objects group-wise
#'
#' The row-wise analogue of `dplyr::bind_rows()` for `sim_stat` objects --
#' e.g. combining a by-region call with a full-domain (`by_region = FALSE`)
#' call into one result, as described in [abundance_by_group()]'s docs.
#' [group_key()] already labels the full-domain result's `region` as
#' `"all"`, so no relabelling is needed before binding.
#'
#' @param ... `sim_stat` objects sharing the same `mats` names and `ncol`.
#'
#' @return A single `sim_stat` with `key`/`mats` row-bound across inputs.
#'
#' @keywords internal
#'
sim_stat_bind <- function(...) {
  parts <- list(...)
  key <- dplyr::bind_rows(lapply(parts, `[[`, "key"))
  mats <- stats::setNames(
    lapply(names(parts[[1]]$mats), function(nm) {
      do.call(rbind, lapply(parts, function(p) p$mats[[nm]]))
    }),
    names(parts[[1]]$mats)
  )
  sim_stat(key, mats)
}

#@TODO: do we need to think more carefully about what "ALL" is, if we are including season?
# - Otherwise are we going to be accidentally double counting? --> I don't think so,
#   as long as region is meant to be the polygon boundaries and NOT included in the model
#   I don't think this is a problem.
# - This is an open question should we be choosing a season? should be we averaging the year?
#   - Does it even matter if we only included season as an intercept in the model?
#   - Would we have enough data to consider this more concretely? - see barn door skate presentation discussion notes
#@TODO: make grouping variable an argument - so that we can provide more options like season
#   - making this an argument would propagate into `abundance_by_group`

#' Grouping key for `sim_draws` matrix aggregation (internal)
#'
#' `rowsum()` requires single grouping vector; this builds one from
#' `time_period`/`region`/`year` (collapsing `region` to a constant `"all"`
#' level for the domain-wide total) plus a lookup table to re-attach the key
#' columns to `rowsum()`'s output by matching on row name.
#'
#' @param grid A `sim_draws$grid` data frame.
#' @param by_region If `FALSE`, `region` is collapsed to a single `"all"`
#'   level -- used to compute the full-domain total without duplicating grid
#'   rows; grouping on a constant column doesn't change the grouping, so the
#'   result is a `sim_stat` whose `key$region` is always `"all"`, ready to
#'   `sim_stat_bind()` alongside a `by_region = TRUE` result.
#'
#' @return A list with `key` (grouping vector, one entry per grid row),
#'   `lookup` (unique key-column combinations with a `.key` character label
#'   matching `rowsum()`'s row names), `key_cols`, `grid` (`grid` with
#'   `region` collapsed to `"all"` when `by_region = FALSE` -- callers that
#'   need row-level, non-key columns alongside the key should subset this,
#'   not the original `grid`, so they see the same collapsed `region`).
#'
#' @keywords internal
group_key <- function(grid, by_region = TRUE) {
  if (!by_region) grid$region <- "all"
  key_cols <- c("time_period", "region", "year")
  key <- interaction(grid[key_cols], drop = TRUE)
  lookup <- unique(grid[key_cols])
  lookup$.key <- as.character(interaction(lookup[key_cols], drop = TRUE))
  list(key = key, lookup = lookup, key_cols = key_cols, grid = grid)
}

#@TODO - this function name isn't very intuitive for all that it does -- make
#   indicator specific wrapper functions, get_COG, get_weighted_mean: e.g.,:
# get_total_abundance <- function(sims, by_region = TRUE) {
#   abundance_by_group(sims, weighted_cols = NULL, by_region = by_region)
# }

# get_centre_of_gravity <- function(sims, by_region = TRUE) {
#   abundance_by_group(sims, weighted_cols = c("X", "Y"), by_region = by_region)
# }

# get_mean_depth <- function(sims, by_region = TRUE) {
#   abundance_by_group(sims, weighted_cols = "depth_m", by_region = by_region)
# }

#' Calculate abundance totals and abundance-weighted means across simulation draws
#'
#' Shared engine for abundance, centre-of-gravity, and mean-depth indicators
#' -- all are a `rowsum()` of `abundance` (optionally weighted by another
#' grid column) over the same grouping key, computed for every simulation
#' draw in one matrix call.
#'
#' @param sims A `sim_draws` object (see [make_sim_draws()]).
#' @param weighted_cols Character vector of `sims$grid` column names to take
#'   the abundance-weighted mean of (e.g. `c("X", "Y")` for centre of
#'   gravity, `"depth_m"` for mean depth). `NULL` (default) returns only
#'   `total_abundance`.
#' @param by_region If `FALSE`, aggregates over the full domain (`region`
#'   collapsed to `"all"` in the grouping key, see [group_key()]).
#'   Call once with `TRUE` and once with `FALSE` and combine with
#'   [sim_stat_bind()] to get both by-region and full-domain totals in one
#'   `sim_stat`.
#'
#' @return A `sim_stat` object: `key` has the grouping columns
#'   (`time_period`, `region` [if `by_region`], `year`), one row per group;
#'   `mats` has `total_abundance` and one entry per `weighted_cols` column,
#'   each a `groups x nsim` matrix.
#'
#' @export
abundance_by_group <- function(sims, weighted_cols = NULL, by_region = TRUE) {
  if (!inherits(sims, "sim_draws")) {
    cli::cli_abort("{.arg sims} must be a {.cls sim_draws} object (see {.fn make_sim_draws}).")
  }
  gk <- group_key(sims$grid, by_region)
  total <- rowsum(exp(sims$log_abundance), gk$key) # groups x nsim; shared denominator

  weighted <- stats::setNames(
    lapply(weighted_cols, \(v) rowsum(exp(sims$log_abundance) * sims$grid[[v]], gk$key) / total),
    weighted_cols
  )

  mats <- c(list(total_abundance = total), weighted)
  key <- gk$lookup[match(rownames(total), gk$lookup$.key), gk$key_cols, drop = FALSE]
  sim_stat(key, mats)
}

#' Mean log density across simulation draws (matrix-native)
#'
#' Unweighted per-cell mean of `sims$log_density`, grouped the same way as
#' [abundance_by_group()]. Unlike summed abundance, a mean of per-cell
#' log-density values isn't dominated by a handful of extreme cells in a
#' single draw -- useful as a cross-check when summed abundance shows
#' inflation that per-cell diagnostics suggest is driven by outlier cells
#' rather than broad-based spatial spread.
#'
#' Cells split across regions (see `01-data-prep.qmd`'s N/S WCVI split)
#' contribute one row per region but keep the parent `cell` id and identical
#' `log_density` in every half. That's correct for `by_region = TRUE` (each
#' half legitimately represents that region's share of the cell), but with
#' `by_region = FALSE` both halves fall into the same domain-wide group and
#' would otherwise be double-counted in an unweighted per-cell mean -- so the
#' full-domain case dedupes to one row per physical (`cell`, `year`) first.
#'
#' @param sims A `sim_draws` object (see [make_sim_draws()]).
#' @param by_region If `FALSE`, aggregates over the full domain (see
#'   [group_key()]), deduplicating split cells first (see above).
#'
#' @return A `sim_stat` object: `key` has the grouping columns; `mats` has a
#'   single `mean_log_density` entry, a `groups x nsim` matrix.
#'
#' @export
log_density_by_group <- function(sims, by_region = TRUE) {
  if (!inherits(sims, "sim_draws")) {
    cli::cli_abort("{.arg sims} must be a {.cls sim_draws} object (see {.fn make_sim_draws}).")
  }
  gk <- group_key(sims$grid, by_region)
  log_density <- sims$log_density
  if (!by_region) {
    keep <- !duplicated(sims$grid[c("cell", "year")])
    gk$key <- gk$key[keep]
    log_density <- log_density[keep, , drop = FALSE]
  }
  total <- rowsum(log_density, gk$key)
  n_g <- as.vector(rowsum(rep(1, nrow(log_density)), gk$key))
  key <- gk$lookup[match(rownames(total), gk$lookup$.key), gk$key_cols, drop = FALSE]
  sim_stat(key, list(mean_log_density = total / n_g))
}

#' Range-edge position across simulation draws (matrix-native)
#'
#' Density-independent sort key (`rank_var`, e.g. `X`/`Y`) doesn't vary by
#' simulation draw, so a single global sort serves all `nsim` columns; the
#' grouped cumulative sum then runs across every sim column at once via
#' `data.table`'s `.SDcols`. Only uses `sims$log_abundance` and grid columns,
#' so this is model-agnostic.
#'
#' @param sims A `sim_draws` object.
#' @param rank_var Tidy-eval column in `sims$grid` to rank/interpolate along
#'   (e.g. `Y`, `X`), unquoted.
#' @param quantiles Named numeric vector of target cumulative-abundance
#'   quantiles, e.g. `c(trailing = 0.05, leading = 0.95)`.
#' @param by_region If `TRUE`, computes edges separately per region;
#'   default `FALSE` matches historical usage (range edge computed over the
#'   full domain only).
#'
#' @return A `sim_stat` object: `key` has `time_period`, `region` (if
#'   `by_region`), `year`, `edge`, one row per group x quantile; `mats` has
#'   a single `{rank_var}_edge` entry, a `nrow(key) x nsim` matrix.
get_range_edge_sims <- function(sims, rank_var, quantiles, by_region = FALSE) {
  if (!inherits(sims, "sim_draws")) {
    cli::cli_abort("{.arg sims} must be a {.cls sim_draws} object (see {.fn make_sim_draws}).")
  }
  rank_name <- rlang::as_name(rlang::ensym(rank_var))
  edge_col <- paste0(rank_name, "_edge")
  gk <- group_key(sims$grid, by_region)
  nsim <- ncol(sims$log_abundance)
  n_q <- length(quantiles)
  simcols <- paste0("sim", seq_len(nsim))

  dt <- data.table::as.data.table(gk$grid[c(gk$key_cols, rank_name)])
  dt[, (simcols) := data.table::as.data.table(exp(sims$log_abundance))]
  data.table::setorderv(dt, c(gk$key_cols, rank_name))
  dt[, (simcols) := lapply(.SD, function(x) { cx <- cumsum(x); cx / cx[.N] }),
     by = eval(gk$key_cols), .SDcols = simcols]

  out <- dt[, {
    edge_mat <- vapply(simcols, function(cn) { # n_q x nsim, colnames = simcols
      stats::approx(get(cn), get(rank_name), xout = quantiles, rule = 2, ties = "ordered")$y
    }, numeric(n_q))
    cbind(data.table::data.table(edge = names(quantiles)), data.table::as.data.table(edge_mat))
  }, by = eval(gk$key_cols)]

  key <- out[, c(gk$key_cols, "edge"), with = FALSE]
  M <- as.matrix(out[, simcols, with = FALSE])
  sim_stat(key, stats::setNames(list(M), edge_col))
}

#' Area-occupied / quantile-EAO across simulation draws (matrix-native)
#'
#' Density (the rank variable here) *does* vary per simulation draw, unlike
#' [get_range_edge_sims()]'s `rank_var`, so this needs a per-sim-column sort:
#' one vectorized `data.table` grouped sort + cumulative sum over the whole
#' grid per draw, looped over `nsim` columns -- not a per-(group, sim) scan.
#'
#' @param sims A `sim_draws` object.
#' @param thresholds Numeric vector of target cumulative-abundance
#'   percentages (0-100), e.g. `c(50, 75, 90, 95)`.
#' @param by_region If `FALSE`, aggregates over the full domain -- see
#'   [abundance_by_group()] for the [sim_stat_bind()] pattern to combine both.
#'
#' @return A `sim_stat` object: `key` has `time_period`, `region` (if
#'   `by_region`), `year`, `threshold`, one row per group x threshold;
#'   `mats` has `area_threshold`, `percent_area_used`, `area_efficiency`,
#'   each a `nrow(key) x nsim` matrix.
get_area_occupied_sims <- function(sims, thresholds, by_region = TRUE) {
  if (!inherits(sims, "sim_draws")) {
    cli::cli_abort("{.arg sims} must be a {.cls sim_draws} object (see {.fn make_sim_draws}).")
  }
  gk <- group_key(sims$grid, by_region)
  nsim <- ncol(sims$log_density)

  per_sim <- purrr::map(seq_len(nsim), function(s) {
    dt <- data.table::as.data.table(gk$grid[gk$key_cols])
    dt[, `:=`(abund = exp(sims$log_abundance[, s]), area = sims$grid$area,
              dens = exp(sims$log_density[, s]))]
    data.table::setorderv(dt, c(gk$key_cols, "dens"),
                           order = c(rep(1L, length(gk$key_cols)), -1L))
    dt[, `:=`(cum_abund = cumsum(abund), cum_area = cumsum(area)), by = eval(gk$key_cols)]

    dt[, {
      total_abund <- cum_abund[.N]; total_area <- sum(area)
      pct <- cum_abund / total_abund * 100
      at <- stats::approx(pct, cum_area, xout = thresholds, rule = 2, ties = "ordered")$y
      data.table::data.table(threshold = thresholds, area_threshold = at,
                              percent_area_used = at / total_area, area_efficiency = at / total_abund)
    }, by = eval(gk$key_cols)]
  })

  key_cols_thr <- c(gk$key_cols, "threshold")
  key <- per_sim[[1]][, key_cols_thr, with = FALSE]
  if (nsim > 1) {
    same_key <- vapply(per_sim[-1], function(dt) {
      identical(as.data.frame(dt[, key_cols_thr, with = FALSE]), as.data.frame(key))
    }, logical(1))
    if (!all(same_key)) {
      cli::cli_abort("Group/threshold ordering differs across simulation draws; cannot combine into a matrix.")
    }
  }

  metrics <- c("area_threshold", "percent_area_used", "area_efficiency")
  mats <- stats::setNames(
    lapply(metrics, function(m) do.call(cbind, lapply(per_sim, `[[`, m))),
    metrics
  )
  sim_stat(key, mats)
}

#' Row-wise central estimate with fast paths for median/mean (internal)
#' `est_function` needs to stay a plain function to match
#' `sdmTMB::get_index_sims()`'s API, but `apply()` defeats the point of a
#' matrix-native engine. This special-cases the two estimators that are
#' actually used and falls back to `apply()` for anything else.
#'
#' @param M A `groups x nsim` matrix.
#' @param f A function taking a numeric vector and returning a scalar.
#'
#' @return A numeric vector, one value per row of `M`.
row_est <- function(M, f) {
  if (identical(f, stats::median)) matrixStats::rowMedians(M)
  else if (identical(f, mean)) rowMeans(M)
  else apply(M, 1L, f)
}

#' Row-wise estimate/95% interval/SE across simulation draws (internal)
#'
#' Shared by [summarise_estimates()] and [summarise_trend()] -- both end in
#' the same operation: collapse a `groups x nsim` matrix to per-row
#' est/lwr/upr/se. `lwr`/`upr` are always the 2.5%/97.5% sim quantiles,
#' independent of `est_function` -- see [summarise_estimates()] for why that's
#' intentional.
#'
#' @param M A `groups x nsim` matrix.
#' @param est_function Function summarising each row -- `stats::median`
#'   (default) or `mean` take the vectorised fast path in [row_est()].
#' @param probs Length-2 vector of quantile probabilities for `lwr`/`upr`
#'   (default `c(0.025, 0.975)`).
#'
#' @return A list with `est`, `lwr`, `upr`, `se`, each a numeric vector of
#'   length `nrow(M)`.
row_stats <- function(M, est_function = stats::median, probs = c(0.025, 0.975)) {
  q <- matrixStats::rowQuantiles(M, probs = probs) # type 7, matches stats::quantile()
  list(
    est = row_est(M, est_function),
    lwr = q[, 1],
    upr = q[, 2],
    se = matrixStats::rowSds(M)
  )
}

#' Estimate/95% interval/SE summary across simulation draws (matrix-native)
#'
#' Collapses a `sim_stat`'s `groups x nsim` matrices directly to per-group
#' estimate/CI summaries.
#'
#' @param x A `sim_stat` object (see [abundance_by_group()] etc.).
#' @param est_function Function summarising each group's central estimate
#'   across sims -- `stats::median` (default) or
#'   `mean` for `E[X]` (e.g. total abundance, following
#'   `sdmTMB::get_index_sims()`'s `est_function` argument and the same
#'   `bias_correct = TRUE` motivation: `project()` returns bare matrices
#'   with no `sd_report`, so TMB's epsilon bias correction isn't available
#'   for forecast years, making the draw mean the only route to `E[X]`
#'   there). `lwr`/`upr` stay the 2.5%/97.5% draw quantiles regardless of
#'   `est_function` -- exactly what `get_index()` does too (`lwr`/`upr` are
#'   built from the uncorrected estimate even when `bias_correct = TRUE`).
#'   Total abundance is a linear functional, so its draw mean is `E[X]`
#'   cleanly; centre of gravity, mean depth, `percent_area_used` and
#'   `area_efficiency` are ratios computed per draw, so `mean` there gives
#'   `E[ratio]`, not the ratio of expectations.
#' @param probs Length-2 vector of quantile probabilities for `lwr`/`upr`
#'   (default `c(0.025, 0.975)`) -- passed to [row_stats()].
#'
#' @return A tibble: `x$key` plus, for each entry of `x$mats`, `est`/`lwr`/
#'   `upr`/`se` (single-metric `x`) or `est_<metric>`/`lwr_<metric>`/...
#'   (multi-metric `x`) -- this naming lets downstream column references
#'   (`est_x`/`est_y`, the EAO `names_pattern = "(est|lwr|upr|se)_(.*)"`
#'   pivot) work directly.
summarise_estimates <- function(x, est_function = stats::median, probs = c(0.025, 0.975)) {
  if (!inherits(x, "sim_stat")) {
    cli::cli_abort("{.arg x} must be a {.cls sim_stat} object (see {.fn abundance_by_group} etc.).")
  }
  multi <- length(x$mats) > 1
  out <- x$key
  for (nm in names(x$mats)) {
    s <- row_stats(x$mats[[nm]], est_function, probs)
    for (f in names(s)) out[[if (multi) paste0(f, "_", nm) else f]] <- s[[f]]
  }
  tibble::as_tibble(out)
}

#' Closed-form grouped-OLS trend summary across simulation draws (matrix-native)
#'
#' `centred_year` is identical across sims within a group, so every group's
#' per-draw OLS regression shares one design vector, and the whole set of
#' regressions reduces to `rowsum()` calls (closed-form simple-OLS algebra).
#'
#' @param x A `sim_stat` object whose `key` contains a `year` column plus
#'   the columns named in `by`. Each entry of `x$mats` is a metric to
#'   regress on `centred_year` -- apply any response transform (`log()`,
#'   `qlogis()`) to the matrix before calling, e.g.
#'   `x$mats$total_abundance <- log(x$mats$total_abundance)`.
#' @param by Character vector of `x$key` columns identifying one time
#'   series per simulation draw (e.g. `c("time_period", "region")`) --
#'   `year` is collapsed into the regression and dropped from the output.
#' @param centre_year Year subtracted from `x$key$year` before fitting
#'   (typically `last_historical_year`).
#' @param est_function Function summarising each coefficient's distribution
#'   across sims -- see [row_stats()].
#'
#' @return A tibble with `term` (`"(Intercept)"`/`"centred_year"`), `metric`
#'   (if `x$mats` has more than one entry), the `by` columns, `estimate`,
#'   `se`, `lwr`, `upr`, `prob_direction` -- one row per term per group (per
#'   metric).
summarise_trend <- function(x, by, centre_year, est_function = stats::median) {
  if (!inherits(x, "sim_stat")) {
    cli::cli_abort("{.arg x} must be a {.cls sim_stat} object (see {.fn abundance_by_group} etc.).")
  }
  k <- interaction(x$key[by], drop = TRUE)
  cy <- x$key$year - centre_year

  n_g_mat <- rowsum(rep(1, length(cy)), k) # G x 1; rownames = sorted group labels
  group_labels <- rownames(n_g_mat)
  n_g   <- as.vector(n_g_mat)
  Sx_g  <- as.vector(rowsum(cy, k))
  Sxx_g <- as.vector(rowsum(cy^2, k))
  denom <- Sxx_g - Sx_g^2 / n_g
  degenerate <- denom == 0 # single-year group: slope unidentifiable, matches lm()'s NA

  lookup <- unique(x$key[by])
  lookup$.key <- as.character(interaction(lookup[by], drop = TRUE))
  key_out <- lookup[match(group_labels, lookup$.key), by, drop = FALSE]

  multi <- length(x$mats) > 1
  build_term <- function(B, term) {
    s <- row_stats(B, est_function)
    tibble::as_tibble(key_out) |>
      dplyr::mutate(
        term = term,
        estimate = s$est, se = s$se, lwr = s$lwr, upr = s$upr,
        prob_direction = 2 * pmin(rowMeans(B > 0), rowMeans(B < 0))
      )
  }

  purrr::map(names(x$mats), function(nm) {
    M <- x$mats[[nm]]
    Sv  <- rowsum(M, k); dimnames(Sv) <- NULL # rowsum() sets rownames = group labels regardless
    Sxv <- rowsum(M * cy, k); dimnames(Sxv) <- NULL # of M's own dimnames -- strip so est/se/etc. below stay unnamed

    slope <- (Sxv - Sx_g * Sv / n_g) / denom
    if (any(degenerate)) slope[degenerate, ] <- NA
    slope_for_intercept <- slope
    slope_for_intercept[degenerate, ] <- 0 # single-year group: intercept is just the group mean
    intercept <- Sv / n_g - slope_for_intercept * (Sx_g / n_g)

    out <- dplyr::bind_rows(
      build_term(intercept, "(Intercept)"),
      build_term(slope, "centred_year")
    )
    if (multi) out$metric <- nm
    out
  }) |> purrr::list_rbind()
}
