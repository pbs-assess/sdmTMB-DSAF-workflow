# Presentation copy of R/03-get-shift-indicators-matrix.R: computes the six DSAF
# shift indicators from the cached sim_draws and saves slide-sized figures
# (historical-only and historical + projection) to presentations/figs/.
# Run R/03-get-shift-indicators-matrix.R first to build the sims cache.

library(dplyr)
library(ggplot2)

source(here::here("R", "00-settings.R"))
source(here::here("R", "functions.R"))

# ---- load-sim-draws ----

model_name <- "dogfish-dgg-rw"
# model_name <- "dogfish-dgg-ar1"

sims_re_path <- file.path(model_out_dir, paste0("sims_re-", model_name, "-matrix.qs2"))
if (!file.exists(sims_re_path)) {
  cli::cli_abort(c(
    "Sim-draws cache {.file {sims_re_path}} not found.",
    "i" = "Run {.file R/03-get-shift-indicators-matrix.R} to create it."
  ))
}
sims <- qs2::qs_read(sims_re_path)

last_historical_year <- max(sims$grid$year[sims$grid$time_period == "historical"])

interval_probs <- c(0.1, 0.9) # 80% simulation interval, as in the other slide figures

# ---- indicators ----
# Each indicator ends up as est/lwr/upr by year with `series` (colour) and
# `panel` (facet) columns so one plotting function handles all six.

# Abundance: mean log density rather than summed abundance -- the sum is
# dominated by a few extreme cells in each draw (see jensens-inequality.R).
abundance <- log_density_by_group(sims, by_region = FALSE) |>
  summarise_estimates(est_function = mean, probs = interval_probs) |>
  transmute(year, est, lwr, upr, series = "all", panel = "Mean log density")

abund_region <- abundance_by_group(sims, by_region = TRUE)
abund_all <- abundance_by_group(sims, by_region = FALSE)
idx_all <- match(
  interaction(abund_region$key$time_period, abund_region$key$year, drop = TRUE),
  interaction(abund_all$key$time_period, abund_all$key$year, drop = TRUE)
)
prop_abundance <- sim_stat(
  abund_region$key,
  list(prop = abund_region$mats$total_abundance / abund_all$mats$total_abundance[idx_all, , drop = FALSE])
) |>
  summarise_estimates(probs = interval_probs) |>
  transmute(year, est, lwr, upr, series = region, panel = "Proportion of total abundance")

cog <- abundance_by_group(sims, weighted_cols = c("X", "Y"), by_region = FALSE)
names(cog$mats) <- c("total_abundance", "x", "y")
cog$mats$total_abundance <- NULL
centre_of_gravity <- summarise_estimates(cog, probs = interval_probs) |>
  tidyr::pivot_longer(
    cols = -c(time_period, region, year),
    names_to = c(".value", "coord"),
    names_pattern = "(est|lwr|upr|se)_(.*)"
  ) |>
  transmute(
    year, est, lwr, upr, series = "all",
    panel = if_else(coord == "x", "Easting (km)", "Northing (km)")
  )

area_occupied <- get_area_occupied_sims(sims, thresholds = 90, by_region = FALSE) |>
  summarise_estimates(probs = interval_probs) |>
  transmute(
    year, est = est_area_threshold, lwr = lwr_area_threshold, upr = upr_area_threshold,
    series = "all", panel = "Area containing 90% of abundance (km²)"
  )

depth <- abundance_by_group(sims, weighted_cols = "depth_m", by_region = FALSE)
names(depth$mats) <- c("total_abundance", "mean_depth")
depth$mats$total_abundance <- NULL
mean_depth <- summarise_estimates(depth, probs = interval_probs) |>
  transmute(year, est, lwr, upr, series = "all", panel = "Abundance-weighted mean depth (m)")

# Y increases north and X east, so the 5% quantile is the south/west edge
range_edge_quantiles <- c(`5% (south / west)` = 0.05, `95% (north / east)` = 0.95)
range_edge <- list(
  "Northing (km)" = get_range_edge_sims(sims, rank_var = Y, quantiles = range_edge_quantiles),
  "Easting (km)" = get_range_edge_sims(sims, rank_var = X, quantiles = range_edge_quantiles)
) |>
  purrr::imap(\(stat, panel) {
    summarise_estimates(stat, probs = interval_probs) |>
      transmute(year, est, lwr, upr, series = factor(edge, levels = names(range_edge_quantiles)), panel = panel)
  }) |>
  bind_rows()

# ---- plots ----

# 16:9 slide is 13.33 x 7.5 in; each figure is sized for ~half a slide
theme_slide <- theme(
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 14),
  strip.text = element_text(size = 12),
  legend.text = element_text(size = 12),
  legend.title = element_blank(),
  legend.position = "top",
  plot.margin = margin(4, 16, 4, 4)
)

plot_indicator <- function(df, y_lab = NULL, reverse_y = FALSE, boundary_year = last_historical_year) {
  p <- ggplot(df) +
    aes(x = year, y = est, ymin = lwr, ymax = upr, colour = series, fill = series) +
    geom_ribbon(alpha = 0.3, colour = NA) +
    geom_line(linewidth = 1) +
    coord_cartesian(expand = FALSE) +
    labs(x = "Year", y = y_lab) +
    theme_slide

  if (max(df$year) > boundary_year) {
    p <- p + geom_vline(xintercept = boundary_year, linetype = "dashed", colour = "grey30")
  }
  if (reverse_y) p <- p + scale_y_reverse()
  if (n_distinct(df$series) == 1) {
    p <- p +
      scale_colour_manual(values = blues[[5]]) +
      scale_fill_manual(values = blues[[3]]) +
      guides(colour = "none", fill = "none")
  } else {
    p <- p + scale_colour_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2")
  }
  if (n_distinct(df$panel) > 1) {
    p <- p + facet_wrap(vars(panel), ncol = 1, scales = "free_y")
  }
  p
}

indicators <- list(
  abundance = list(df = abundance, y_lab = "Mean log density"),
  `centre-of-gravity` = list(df = centre_of_gravity),
  `area-occupied` = list(df = area_occupied, y_lab = "Area occupied, 90% of abundance (km²)"),
  depth = list(df = mean_depth, y_lab = "Abundance-weighted mean depth (m)", reverse_y = TRUE),
  `range-edge` = list(df = range_edge),
  `proportional-abundance` = list(df = prop_abundance, y_lab = "Proportion of total abundance")
)

periods <- c(historical = last_historical_year, projection = Inf)

plots <- purrr::map(indicators, \(ind) {
  purrr::map(periods, \(max_year) {
    plot_indicator(
      filter(ind$df, year <= max_year),
      y_lab = ind$y_lab %||% NULL,
      reverse_y = ind$reverse_y %||% FALSE
    )
  })
})

# ---- save ----

fig_dir <- here::here("presentations", "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

purrr::iwalk(plots, \(by_period, indicator) {
  purrr::iwalk(by_period, \(p, period) {
    ggsave(
      file.path(fig_dir, paste0("shift-", indicator, "-", period, ".png")),
      p, width = 6.5, height = 4.5, dpi = 300, bg = "white"
    )
  })
})
