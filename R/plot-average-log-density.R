library(dplyr)
library(ggplot2)

source(here::here("R", "00-settings.R"))
source(here::here("R", "functions.R"))

fig_dir <- here::here('presentations', 'figs')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

sims <- sims_re

# sims <- readRDS(file.path(model_out_dir, "sims.rds"))
last_historical_year <- max(sims$year[sims$time_period == "historical"])

average_log_density_by_sim <- sims |>
  filter(region == "all") |>
  group_by(year, sim) |>
  summarise(
    average_log_density = mean(log_density, na.rm = TRUE),
    .groups = "drop"
  )

average_log_density_summary <- average_log_density_by_sim |>
  group_by(year) |>
  summarise(
    mean = mean(average_log_density),
    q10 = quantile(average_log_density, 0.10),
    q90 = quantile(average_log_density, 0.90),
    .groups = "drop"
  )

set.seed(2026)
example_sims <- sample(
  unique(average_log_density_by_sim$sim),
  size = min(40L, n_distinct(average_log_density_by_sim$sim))
)

example_draws <- average_log_density_by_sim |>
  filter(sim %in% example_sims)

p_average_log_density <- ggplot(
  average_log_density_summary,
  aes(x = year, y = mean)
) +
  geom_ribbon(
    aes(ymin = q10, ymax = q90),
    fill = blues[[2]],
    alpha = 0.35
  ) +
  geom_line(
    data = example_draws,
    aes(y = average_log_density, group = sim),
    colour = blues[[3]],
    alpha = 0.25,
    linewidth = 0.35
  ) +
  geom_line(colour = blues[[5]], linewidth = 1) +
  geom_vline(
    xintercept = last_historical_year,
    linetype = "dashed",
    colour = "grey40"
  ) +
  labs(
    x = "Year",
    y = "Average log density",
    subtitle = NULL
  ) +
  coord_cartesian(expand = FALSE, xlim = c(NA, 2035.3))

p_average_log_density

ggsave(file.path(fig_dir, "dogfish-dgg-rw-projection.png"),
  width = 4.55, height = 3.25
)

# Historical-only version, no example draws: reuse the built plot rather than
# rebuilding it from scratch.
p_average_log_density_historical <- p_average_log_density
p_average_log_density_historical$layers[[2]] <- NULL # drop example draws layer
p_average_log_density_historical <- p_average_log_density_historical +
  filter(average_log_density_summary, year <= last_historical_year)
p_average_log_density_historical$layers[[3]] <- NULL # drop vline (no future cutoff to mark)
p_average_log_density_historical <- p_average_log_density_historical +
  coord_cartesian(expand = FALSE, xlim = c(NA, last_historical_year))

p_average_log_density_historical

ggsave(
  file.path(fig_dir, "dogfish-dgg-rw-projection_historical-single-line.png"),
  plot = p_average_log_density_historical,
  width = 4.55, height = 3.25
)

# By region
average_log_density_by_sim_by_region <- sims |>
  filter(year <= last_historical_year) |>
  filter(region != "all") |>
  group_by(year, region, sim) |>
  summarise(
    average_log_density = mean(log_density, na.rm = TRUE),
    .groups = "drop"
  )

average_log_density_summary_by_region <- average_log_density_by_sim_by_region |>
  group_by(year, region) |>
  summarise(
    mean = mean(average_log_density),
    q10 = quantile(average_log_density, 0.10),
    q90 = quantile(average_log_density, 0.90),
    .groups = "drop"
  )

set.seed(2026)
p_average_log_density_by_region <- ggplot(average_log_density_summary_by_region) +
  aes(x = year, y = mean, colour = region, fill = region, group = region) +
  geom_ribbon(
    aes(ymin = q10, ymax = q90),
    alpha = 0.35,
    colour = NA
  ) +
  geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    x = "Year",
    y = "Average log density",
    subtitle = NULL
  ) +
  coord_cartesian(expand = FALSE, xlim = c(NA, last_historical_year)) +
  guides(colour = "none", fill = "none")

p_average_log_density_by_region

ggsave(file.path(fig_dir, "dogfish-dgg-rw-projection_by_region.png"), width = 4.55, height = 3.25)

# --------

range_edge_quantiles <- c(trailing = 0.20, leading = 0.80)

range_edges_by_sim <- sims |>
  filter(region == "all") |>
  get_range_edge_sim(
    rank_var = Y,
    quantiles = range_edge_quantiles,
    by = c(time_period, region, year, sim)
  ) |>
  transmute(time_period, region, year, sim, edge, edge_position = Y_edge)

range_edge_summary <- range_edges_by_sim |>
  group_by(time_period, region, year, edge) |>
  summarise(
    mean = mean(edge_position, na.rm = TRUE),
    q20 = quantile(edge_position, 0.20, na.rm = TRUE),
    q80 = quantile(edge_position, 0.80, na.rm = TRUE),
    .groups = "drop"
  )

p_range_edges <- ggplot(
  range_edge_summary,
  aes(x = year, y = mean, colour = edge, fill = edge)
) +
  geom_ribbon(aes(ymin = q20, ymax = q80), alpha = 0.25, colour = NA) +
  geom_line() +
  geom_vline(
    xintercept = last_historical_year,
    linetype = "dashed",
    colour = "grey40"
  ) +
  labs(
    x = "Year",
    y = "North–south range edge position (km)",
    colour = "Range edge",
    fill = "Range edge",
    subtitle = "Mean and 80% simulation interval; edges at 20% and 80% abundance"
  ) +
  coord_cartesian(expand = FALSE) +
  scale_colour_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2")

# p_range_edges +
#   theme(legend.position = "top")
# ggsave(file.path(fit_dir, "dogfish-dgg-range-edge.png"), width = )


mean_depth_by_sim <- sims_fe |>
  filter(region == "all") |>
  group_by(year, sim) |>
  summarise(
    mean_depth = weighted.mean(
      depth_m,
      w = exp(log_abundance),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

mean_depth_summary <- mean_depth_by_sim |>
  group_by(year) |>
  summarise(
    mean = mean(mean_depth, na.rm = TRUE),
    q10 = quantile(mean_depth, 0.10, na.rm = TRUE),
    q90 = quantile(mean_depth, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

mean_depth_examples <- mean_depth_by_sim |>
  filter(sim %in% example_sims)

p_mean_depth <- ggplot(
  mean_depth_summary,
  aes(x = year, y = mean)
) +
  geom_ribbon(
    aes(ymin = q10, ymax = q90),
    fill = blues[[2]],
    alpha = 0.35
  ) +
  geom_line(
    data = mean_depth_examples,
    aes(y = mean_depth, group = sim),
    colour = blues[[3]],
    alpha = 0.25,
    linewidth = 0.35
  ) +
  geom_line(colour = blues[[5]], linewidth = 1) +
  geom_vline(
    xintercept = last_historical_year,
    linetype = "dashed",
    colour = "grey40"
  ) +
  scale_y_reverse() +
  labs(
    x = "Year",
    y = "Density-weighted mean depth (m)",
    subtitle = "Mean, 80% simulation interval, and 40 example draws"
  ) +
  coord_cartesian(expand = FALSE)

p_mean_depth

ggsave(file.path(fig_dir, "range-edges-dog-eg.png"), width = 6, height = 4)
