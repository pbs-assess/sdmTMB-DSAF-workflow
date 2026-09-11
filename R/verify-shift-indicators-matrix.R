# Verify the matrix-native shift-indicator engine (R/functions2.R) against
# the long-format pipeline currently in R/03-get-shift-indicators.qmd
# (R/functions.R). See ~/.claude/plans/we-have-implemented-the-parallel-orbit.md.
#
# Both paths start from the same `project()` draws and are compared with
# waldo::compare() -- correctness, not speed, is the point of this script
# (nsim kept small on purpose; bump it separately for a timing check).
#
# Covers both the per-sim indicator draws (functions2.R's sim_stat engine)
# AND the uncertainty summaries built on top of them (summarise_estimates(),
# summarise_trend()) -- the two things the matrix-native refactor was
# meant to preserve while making fast.

library(dplyr)
library(sdmTMB)

source(here::here("R", "00-settings.R"))
source(here::here("R", "functions.R"))
# functions2.R no longer reuses these names, so aliasing is only for
# readability (matches this script's `_old`/`_new` naming throughout).
get_range_edge_sim_old <- get_range_edge_sim
summarise_sim_ci_old <- summarise_sim_ci
summarise_sim_trend_old <- summarise_sim_trend
source(here::here("R", "functions2.R"))

model_name <- "dogfish-dgg-rw"
fit_path <- file.path(fit_dir, paste0(model_name, ".rds"))
fit <- readRDS(fit_path)

prediction_grid <- readRDS(file.path(model_input_dir, "prediction-grid.rds")) |>
  filter(year <= 2035)

last_historical_year <- max(fit$data$year)

nsim <- 50
set.seed(2026)
projections <- project(
  fit,
  newdata = prediction_grid,
  nsim = nsim,
  sample_fe = FALSE,
  sample_historical_re = TRUE,
  sample_future_re = TRUE,
  offset = rep(0, nrow(prediction_grid))
)

# ---- OLD path (functions.R: make_sims_df + long-format group_by/summarise) ----

time_old <- system.time({
log_density <- return_density(projections, fit, log = TRUE)
sims_old <- make_sims_df(prediction_grid, log_density, last_historical_year)

a_by_sim_old <- sims_old |>
  group_by(time_period, region, year, sim) |>
  summarise(total_abundance = sum(exp(log_abundance)), .groups = "drop")

totals_by_sim_old <- a_by_sim_old |>
  filter(region == "all") |>
  select(time_period, year, sim, total_abundance_all = total_abundance)

p_abund_by_sim_old <- a_by_sim_old |>
  filter(region != "all") |>
  left_join(totals_by_sim_old, by = c("time_period", "year", "sim")) |>
  mutate(prop = total_abundance / total_abundance_all)

cog_by_sim_old <- sims_old |>
  group_by(time_period, region, year, sim) |>
  summarise(
    x = weighted.mean(X, w = exp(log_abundance)),
    y = weighted.mean(Y, w = exp(log_abundance)),
    .groups = "drop"
  )

d_by_sim_old <- sims_old |>
  group_by(time_period, region, year, sim) |>
  summarise(mean_depth = weighted.mean(depth_m, w = exp(log_abundance)), .groups = "drop")

thresholds <- c(50, 75, 90, 95)

ranked_accumulation_by_sims_old <- sims_old |>
  group_by(time_period, region, year, sim) |>
  arrange(desc(log_density), .by_group = TRUE) |>
  mutate(
    total_abundance = sum(exp(log_abundance)),
    total_area = sum(area),
    cumulative_abundance = cumsum(exp(log_abundance)),
    percent_abundance = cumulative_abundance / total_abundance * 100,
    cumulative_area = cumsum(area)
  )

quantile_eao_by_sim_old <- ranked_accumulation_by_sims_old |>
  reframe(
    threshold = thresholds,
    area_threshold = approx(x = percent_abundance, y = cumulative_area,
                             xout = thresholds, rule = 2, ties = "ordered")$y,
    percent_area_used = area_threshold / dplyr::first(total_area),
    area_efficiency = area_threshold / dplyr::first(total_abundance)
  ) |>
  select(time_period, region, year, sim, threshold, area_threshold, percent_area_used, area_efficiency)

range_edge_quantiles <- c(trailing = 0.05, leading = 0.95)

range_edge_NS_by_sim_old <- sims_old |>
  filter(region == "all") |>
  get_range_edge_sim_old(
    rank_var = Y,
    quantiles = range_edge_quantiles,
    by = c(time_period, region, year, sim)
  )
range_edge_EW_by_sim_old <- sims_old |>
  filter(region == "all") |>
  get_range_edge_sim_old(
    rank_var = X,
    quantiles = range_edge_quantiles,
    by = c(time_period, region, year, sim)
  )
}) # time_old

# ---- NEW path (functions2.R: make_sim_draws + matrix-native engine) ----

time_new <- system.time({
log_density <- return_density(projections, fit, log = TRUE)
sims_new <- make_sim_draws(prediction_grid, log_density, last_historical_year)

abund_region_new <- abundance_by_group(sims_new)
abund_all_raw_new <- abundance_by_group(sims_new, by_region = FALSE)
a_stat_new <- sim_stat_bind(abund_region_new, abund_all_raw_new)

idx_all <- match(
  interaction(abund_region_new$key$time_period, abund_region_new$key$year, drop = TRUE),
  interaction(abund_all_raw_new$key$time_period, abund_all_raw_new$key$year, drop = TRUE)
)
prop_stat_new <- sim_stat(
  abund_region_new$key,
  list(prop = abund_region_new$mats$total_abundance / abund_all_raw_new$mats$total_abundance[idx_all, , drop = FALSE])
)

cog_stat_new <- sim_stat_bind(
  abundance_by_group(sims_new, weighted_cols = c("X", "Y")),
  abundance_by_group(sims_new, weighted_cols = c("X", "Y"), by_region = FALSE)
)
names(cog_stat_new$mats) <- c("total_abundance", "x", "y")
cog_stat_new$mats$total_abundance <- NULL

d_stat_new <- sim_stat_bind(
  abundance_by_group(sims_new, weighted_cols = "depth_m"),
  abundance_by_group(sims_new, weighted_cols = "depth_m", by_region = FALSE)
)
names(d_stat_new$mats) <- c("total_abundance", "mean_depth")
d_stat_new$mats$total_abundance <- NULL

eao_stat_new <- sim_stat_bind(
  get_area_occupied_sims(sims_new, thresholds),
  get_area_occupied_sims(sims_new, thresholds, by_region = FALSE)
)

range_ns_stat_new <- get_range_edge_sims(sims_new, rank_var = Y, quantiles = range_edge_quantiles)
range_ew_stat_new <- get_range_edge_sims(sims_new, rank_var = X, quantiles = range_edge_quantiles)

a_by_sim_new <- sim_stat_tibble(a_stat_new)
cog_by_sim_new <- sim_stat_tibble(cog_stat_new)
d_by_sim_new <- sim_stat_tibble(d_stat_new)
quantile_eao_by_sim_new <- sim_stat_tibble(eao_stat_new) |>
  select(time_period, region, year, sim, threshold, area_threshold, percent_area_used, area_efficiency)
range_edge_NS_by_sim_new <- sim_stat_tibble(range_ns_stat_new)
range_edge_EW_by_sim_new <- sim_stat_tibble(range_ew_stat_new)
}) # time_new

cli::cli_h2("Timing (indicator computation only, project() excluded)")
cli::cli_alert_info("nsim = {nsim}")
cli::cli_alert_info("old (long-format pivot): {round(time_old[['elapsed']], 2)}s")
cli::cli_alert_info("new (matrix-native): {round(time_new[['elapsed']], 2)}s")

# ---- Compare (per-sim draws) ----

key_cols_group <- c("time_period", "region", "year", "sim")
key_cols_edge <- c("time_period", "region", "year", "sim", "edge")

compare_one <- function(name, old, new, key_cols) {
  # stats::quantile() names its output ("2.5%"/"97.5%"); matrixStats::rowQuantiles()
  # doesn't -- unname() so that incidental attribute doesn't register as a mismatch.
  unname_cols <- function(df) dplyr::mutate(df, across(everything(), unname))
  cmp <- waldo::compare(
    old |> arrange(across(all_of(key_cols))) |> select(all_of(names(new))) |> unname_cols(),
    new |> arrange(across(all_of(key_cols))) |> unname_cols(),
    tolerance = 1e-8
  )
  cli::cli_h2(name)
  if (length(cmp) == 0) {
    cli::cli_alert_success("identical")
  } else {
    print(cmp)
  }
  invisible(cmp)
}

compare_one("abundance (a_by_sim)", a_by_sim_old, a_by_sim_new, key_cols_group)
compare_one("centre of gravity (cog_by_sim)", cog_by_sim_old, cog_by_sim_new, key_cols_group)
compare_one("mean depth (d_by_sim)", d_by_sim_old, d_by_sim_new, key_cols_group)
compare_one("area occupied / EAO (quantile_eao_by_sim)", quantile_eao_by_sim_old, quantile_eao_by_sim_new,
            c(key_cols_group, "threshold"))
compare_one("range edge N/S", range_edge_NS_by_sim_old, range_edge_NS_by_sim_new, key_cols_edge)
compare_one("range edge E/W", range_edge_EW_by_sim_old, range_edge_EW_by_sim_new, key_cols_edge)

# ---- Uncertainty summaries: median CI (default est_function) ----

key_cols_ci <- c("time_period", "region", "year")

time_old_stats <- system.time({
abund_ci_old <- a_by_sim_old |> summarise_sim_ci_old(cols = total_abundance, by = c(time_period, region, year))
p_abund_ci_old <- p_abund_by_sim_old |> summarise_sim_ci_old(cols = prop, by = c(time_period, region, year))
cog_ci_old <- cog_by_sim_old |> summarise_sim_ci_old(c(x, y), c(time_period, region, year))
depth_ci_old <- d_by_sim_old |> summarise_sim_ci_old(mean_depth, c(time_period, region, year))
eao_ci_old <- quantile_eao_by_sim_old |> summarise_sim_ci_old(
  cols = c(area_threshold, percent_area_used, area_efficiency),
  by = c(time_period, region, year, threshold)
)
range_ns_ci_old <- range_edge_NS_by_sim_old |> summarise_sim_ci_old(Y_edge, c(time_period, region, year, edge))
range_ew_ci_old <- range_edge_EW_by_sim_old |> summarise_sim_ci_old(X_edge, c(time_period, region, year, edge))

abund_trend_old <- a_by_sim_old |> mutate(centred_year = year - last_historical_year) |>
  summarise_sim_trend_old(value = log(total_abundance), by = c(time_period, region))
p_abund_trend_old <- p_abund_by_sim_old |> mutate(centred_year = year - last_historical_year) |>
  summarise_sim_trend_old(value = qlogis(prop), by = c(time_period, region))
cog_trend_old <- cog_by_sim_old |> mutate(centred_year = year - last_historical_year) |>
  tidyr::pivot_longer(c(x, y), names_to = "coord", values_to = "cog") |>
  summarise_sim_trend_old(value = cog, by = c(time_period, region, coord))
depth_trend_old <- d_by_sim_old |> mutate(centred_year = year - last_historical_year) |>
  summarise_sim_trend_old(value = mean_depth, by = c(time_period, region))
eao_trend_old <- quantile_eao_by_sim_old |> mutate(
    centred_year = year - last_historical_year,
    percent_area_used = qlogis(percent_area_used)
  ) |>
  tidyr::pivot_longer(c(area_threshold, percent_area_used, area_efficiency),
                       names_to = "metric", values_to = "value") |>
  summarise_sim_trend_old(value = value, by = c(time_period, region, metric, threshold))
range_ns_trend_old <- range_edge_NS_by_sim_old |> mutate(centred_year = year - last_historical_year) |>
  summarise_sim_trend_old(value = Y_edge, by = c(time_period, region, edge))
range_ew_trend_old <- range_edge_EW_by_sim_old |> mutate(centred_year = year - last_historical_year) |>
  summarise_sim_trend_old(value = X_edge, by = c(time_period, region, edge))
}) # time_old_stats

time_new_stats <- system.time({
abund_ci_new <- summarise_estimates(a_stat_new)
p_abund_ci_new <- summarise_estimates(prop_stat_new)
cog_ci_new <- summarise_estimates(cog_stat_new)
depth_ci_new <- summarise_estimates(d_stat_new)
eao_ci_new <- summarise_estimates(eao_stat_new)
range_ns_ci_new <- summarise_estimates(range_ns_stat_new)
range_ew_ci_new <- summarise_estimates(range_ew_stat_new)

abund_trend_stat <- a_stat_new
abund_trend_stat$mats$total_abundance <- log(abund_trend_stat$mats$total_abundance)
abund_trend_new <- summarise_trend(abund_trend_stat, by = c("time_period", "region"),
                                    centre_year = last_historical_year)

prop_trend_stat <- prop_stat_new
prop_trend_stat$mats$prop <- qlogis(prop_trend_stat$mats$prop)
p_abund_trend_new <- summarise_trend(prop_trend_stat, by = c("time_period", "region"),
                                      centre_year = last_historical_year)

cog_trend_new <- summarise_trend(cog_stat_new, by = c("time_period", "region"),
                                  centre_year = last_historical_year) |>
  dplyr::rename(coord = metric)

depth_trend_new <- summarise_trend(d_stat_new, by = c("time_period", "region"),
                                    centre_year = last_historical_year)

eao_trend_stat <- eao_stat_new
eao_trend_stat$mats$percent_area_used <- qlogis(eao_trend_stat$mats$percent_area_used)
eao_trend_new <- summarise_trend(eao_trend_stat, by = c("time_period", "region", "threshold"),
                                  centre_year = last_historical_year)

range_ns_trend_new <- summarise_trend(range_ns_stat_new, by = c("time_period", "region", "edge"),
                                       centre_year = last_historical_year)
range_ew_trend_new <- summarise_trend(range_ew_stat_new, by = c("time_period", "region", "edge"),
                                       centre_year = last_historical_year)
}) # time_new_stats

cli::cli_h2("Timing (uncertainty summaries: summarise_estimates() + summarise_trend())")
cli::cli_alert_info("nsim = {nsim}")
cli::cli_alert_info("old (per-sim lm() loop): {round(time_old_stats[['elapsed']], 2)}s")
cli::cli_alert_info("new (closed-form rowsum()): {round(time_new_stats[['elapsed']], 2)}s")

cli::cli_h2("Uncertainty summaries: median CI (default est_function)")
compare_one("abundance CI (median)", abund_ci_old, abund_ci_new, key_cols_ci)
compare_one("proportional abundance CI (median)", p_abund_ci_old, p_abund_ci_new, key_cols_ci)
compare_one("centre of gravity CI (median)", cog_ci_old, cog_ci_new, key_cols_ci)
compare_one("mean depth CI (median)", depth_ci_old, depth_ci_new, key_cols_ci)
compare_one("area occupied / EAO CI (median)", eao_ci_old, eao_ci_new, c(key_cols_ci, "threshold"))
compare_one("range edge N/S CI (median)", range_ns_ci_old, range_ns_ci_new, c(key_cols_ci, "edge"))
compare_one("range edge E/W CI (median)", range_ew_ci_old, range_ew_ci_new, c(key_cols_ci, "edge"))

# ---- Uncertainty summaries: mean (E[X]) branch, checked against an inline reference ----
# (the old function is hardcoded to the median, so there's no old `_mean` to shadow --
# build the reference directly from the old long-format tibbles instead)

cli::cli_h2("Uncertainty summaries: mean (E[X]) branch")

abund_ci_mean_old <- a_by_sim_old |>
  dplyr::summarise(
    est = mean(total_abundance),
    lwr = stats::quantile(total_abundance, 0.025),
    upr = stats::quantile(total_abundance, 0.975),
    se  = stats::sd(total_abundance),
    .by = c(time_period, region, year)
  )
compare_one("abundance CI (mean)", abund_ci_mean_old, summarise_estimates(a_stat_new, est_function = mean), key_cols_ci)

cog_ci_mean_old <- cog_by_sim_old |>
  dplyr::summarise(
    est_x = mean(x), lwr_x = stats::quantile(x, 0.025), upr_x = stats::quantile(x, 0.975), se_x = stats::sd(x),
    est_y = mean(y), lwr_y = stats::quantile(y, 0.025), upr_y = stats::quantile(y, 0.975), se_y = stats::sd(y),
    .by = c(time_period, region, year)
  )
compare_one("centre of gravity CI (mean)", cog_ci_mean_old, summarise_estimates(cog_stat_new, est_function = mean), key_cols_ci)

cli::cli_h2("row_est() apply() fallback vs mean() fast path")
fallback_cmp <- waldo::compare(
  summarise_estimates(a_stat_new, est_function = mean),
  summarise_estimates(a_stat_new, est_function = function(z) mean(z)), # not identical() to mean -- takes the apply() fallback
  tolerance = 1e-8
)
if (length(fallback_cmp) == 0) {
  cli::cli_alert_success("identical")
} else {
  print(fallback_cmp)
}

# ---- Uncertainty summaries: trend (closed-form OLS vs lm() loop) ----

cli::cli_h2("Uncertainty summaries: trend")
compare_one("abundance trend", abund_trend_old, abund_trend_new, c("term", "time_period", "region"))
compare_one("proportional abundance trend", p_abund_trend_old, p_abund_trend_new, c("term", "time_period", "region"))
compare_one("centre of gravity trend", cog_trend_old, cog_trend_new, c("term", "time_period", "region", "coord"))
compare_one("mean depth trend", depth_trend_old, depth_trend_new, c("term", "time_period", "region"))
compare_one("area occupied / EAO trend", eao_trend_old, eao_trend_new,
            c("term", "time_period", "region", "metric", "threshold"))
compare_one("range edge N/S trend", range_ns_trend_old, range_ns_trend_new, c("term", "time_period", "region", "edge"))
compare_one("range edge E/W trend", range_ew_trend_old, range_ew_trend_new, c("term", "time_period", "region", "edge"))
