library(dplyr)
library(ggplot2)

source(here::here("R", "00-settings.R"))
source(here::here("R", "functions.R"))

fig_dir <- here::here('presentations', 'figs')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

prediction_grid <- readRDS(file.path(model_input_dir, "prediction-grid.rds")) |>
  filter(year <= 2050)

# title = "Dogfish"
title = "Pacific Cod"

if (title == "Pacific Cod") {
  dat <- readRDS(file.path(model_input_dir, "pacific-cod.rds")) |>
    filter(!is.na(log_depth)) |>
    mutate(offset = log(area_km2))
}
if (title == "Dogfish") {
  dat <- readRDS(file.path(model_input_dir, "dogfish.rds"))
}

NSIM <- 400

dat |>
  mutate(presence = ifelse(catch_weight == 0, 0, 1)) |>
  group_by(year) |>
  summarise(total_sets = n(), positive_catch = sum(presence)) |>
  mutate(mean_sets = mean(positive_catch / total_sets))

# A tibble: 11 × 3
#     year total_sets positive_catch
#    <int>      <int>          <dbl>
#  1  2004         87             53
#  2  2006        164             88
#  3  2008        159             64
#  4  2010        124             92
#  5  2012        150             94
#  6  2014        146            110
#  7  2016        140             99
#  8  2018        190             91
#  9  2021        168             86
# 10  2022        126             79
# 11  2024        148            106


mesh <- make_mesh(dat, c("X", "Y"), cutoff = 15)
missing_years <- min(dat$year):max(dat$year)
last_historical_year <- max(dat$year)

# hist(dat$catch_weight)

fit_dgg_rw <-
  sdmTMB(
    formula = catch_weight ~ 1 + poly(log_depth, 2) + poly(bottom_temp, 2),
    family = delta_gengamma(type = "poisson-link"),
    data = dat,
    mesh = mesh,
    spatial = "on",
    spatiotemporal = "rw",
    time = "year",
    extra_time = missing_years,
    offset = "offset",
    control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
    silent = FALSE
  )
sanity(fit_dgg_rw)

fit_dgg_ar1 <- update(fit_dgg_rw, spatiotemporal = "ar1")
sanity(fit_dgg_ar1)

fit_rw <- fit_dgg_rw
fit_ar1 <- fit_dgg_ar1

get_proj_draws <- function(fit, nsim = 400,
  sample_fe = FALSE, sample_historical_re = TRUE, sample_future_re = TRUE,
  future_re = "include", model_name, seed = 107) {
  last_historical_year = max(fit$data$year)
  set.seed(seed)
  proj_mat <- project(
      fit,
      newdata = prediction_grid,
      nsim = nsim,
      sample_fe = sample_fe,
      sample_historical_re = sample_historical_re,
      sample_future_re = sample_future_re,
      future_re = future_re,
      offset = rep(0, nrow(prediction_grid))
    ) |>
    return_density(fit, log = TRUE)

  prediction_grid |>
    mutate(sample_fe = sample_fe,
        sample_historical_re = sample_historical_re,
        sample_future_re = sample_future_re) |>
    make_sim_draws(proj_mat, last_historical_year, model_name)
}

sims_rw <- get_proj_draws(fit_rw, nsim = NSIM, model_name = "rw",
  sample_fe = FALSE,
  sample_historical_re = TRUE,
  sample_future_re = TRUE
)

sims_ar1 <- get_proj_draws(fit_ar1, nsim = NSIM, model_name = "ar1",
  sample_fe = FALSE,
  sample_historical_re = TRUE,
  sample_future_re = TRUE)

sims_rw_fe <- get_proj_draws(fit_rw, nsim = NSIM, model_name = "rw",
  sample_fe = FALSE,
  sample_historical_re = TRUE,
  sample_future_re = FALSE,
  future_re = "fix"
)

sims_ar1_fe <- get_proj_draws(fit_ar1, nsim = NSIM, model_name = "ar1",
  sample_fe = FALSE,
  sample_historical_re = TRUE,
  sample_future_re = FALSE,
  future_re = "fix"
)
sims_rw_fe_zero <- get_proj_draws(fit_rw, nsim = NSIM, model_name = "rw",
  sample_fe = FALSE,
  sample_historical_re = TRUE,
  sample_future_re = FALSE,
  future_re = "zero"
)

sims_ar1_fe_zero <- get_proj_draws(fit_ar1, nsim = NSIM, model_name = "ar1",
  sample_fe = FALSE,
  sample_historical_re = TRUE,
  sample_future_re = FALSE,
  future_re = "zero"
)


ld_all_rw <- log_density_by_group(sims_rw, by_region = FALSE) |>
  summarise_estimates(est_function = mean, probs = c(0.1, 0.9)) |>
  mutate(type = "RW", future_re = "include")
ld_all_ar1 <- log_density_by_group(sims_ar1, by_region = FALSE) |>
  summarise_estimates(est_function = mean, probs = c(0.1, 0.9)) |>
  mutate(type = "AR1", future_re = "include")
ld_all_rw_fe <- log_density_by_group(sims_rw_fe, by_region = FALSE) |>
  summarise_estimates(est_function = mean, probs = c(0.1, 0.9)) |>
  mutate(type = "RW", future_re = "fixed")
ld_all_ar1_fe <- log_density_by_group(sims_ar1_fe, by_region = FALSE) |>
  summarise_estimates(est_function = mean, probs = c(0.1, 0.9)) |>
  mutate(type = "AR1", future_re = "fixed")
ld_all_rw_zero <- log_density_by_group(sims_rw_fe_zero, by_region = FALSE) |>
  summarise_estimates(est_function = mean, probs = c(0.1, 0.9)) |>
  mutate(type = "RW", future_re = "zero")
ld_all_ar1_zero <- log_density_by_group(sims_ar1_fe_zero, by_region = FALSE) |>
  summarise_estimates(est_function = mean, probs = c(0.1, 0.9)) |>
  mutate(type = "AR1", future_re = "zero")

source('R/scale-indices.R')
test1 <- bind_rows(ld_all_rw, ld_all_ar1) |>
  scale_indices(.groups = "type", scale_by = "geomean", space = "log", min_year = last_historical_year, max_year = last_historical_year) |>
  mutate(type = factor(type, levels = c("RW", "AR1")))
test2 <- bind_rows(ld_all_rw_fe, ld_all_ar1_fe) |>
  scale_indices(.groups = "type", scale_by = "geomean", space = "log", min_year = last_historical_year, max_year = last_historical_year) |>
  mutate(type = factor(type, levels = c("RW", "AR1")))
test3 <- bind_rows(ld_all_rw, ld_all_rw_fe)#, ld_all_rw_zero)
test4 <- bind_rows(ld_all_ar1, ld_all_ar1_fe)#, ld_all_ar1_zero)

pal_fill <- c("#E8A87C", "#4E79A7", "black")
pal_line <- c("#C1703B", "#2C4A6E", "black")

plot_mean_ld <- function(df, title, colour_var, title_var) {
  ggplot(df) +
  aes(x = year, y = est, ymin = lwr, ymax = upr,
    colour = .data[[colour_var]], fill = .data[[colour_var]]
  ) +
  # geom_line() +
  geom_line(linewidth = 1) +
  geom_ribbon(colour = NA, alpha = 0.3) +
  geom_vline(xintercept = last_historical_year, linetype = "dashed", colour = "grey30") +
  scale_fill_manual(values = pal_fill) +
  scale_color_manual(values = pal_line) +
  coord_cartesian(expand = FALSE) +
  labs(
    x = "Year",
    y = "Scaled log density",
    subtitle = paste0(title, " - ", unique(df[[title_var]]))
  )
}

p1 <- plot_mean_ld(test1, title = title,colour_var = "type", title_var = "future_re")#; p1
p2 <- plot_mean_ld(test2, title = title, colour_var = "type", title_var = "future_re")#; p2
p3 <- plot_mean_ld(test3, title = title, colour_var = "future_re", title_var = "type")#; p3
p4 <- plot_mean_ld(test4, title = title, colour_var = "future_re", title_var = "type")#; p4

# (p1 + p2) / (p3 + p4)

# Right half of a 16:9 slide (13.33 x 7.5 in): ~6 x 4.5 in, with slide-legible text
p1_slide <- p1 +
  labs(y = "Scaled log density", subtitle = NULL) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    # legend.position = "inside",
    # legend.position.inside = c(0.02, 0.02),
    # legend.justification = c(0, 0),
    legend.key.spacing.y = unit(0.5, "lines"),
    plot.margin = margin(0, 0, 0, 0)
  )

ggsave(
  file.path(fig_dir, "pacific-cod-log-density-rw-ar1.png"),
  plot = p1_slide,
  width = 6.5, height = 4.5, dpi = 300
)

# Mean Bottom temperature trend
prediction_grid |>
  group_by(year) |>
  summarise(mean_bt = mean(bottom_temp)) |>
ggplot(data = _) +
  aes(x = year, y = mean_bt) +
  stat_smooth() +
  geom_point() +
  labs(y = "Mean bottom temperature (˚C)", x = "Year")


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
