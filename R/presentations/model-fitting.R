library(dplyr)
library(ggplot2)
library(sf)
library(ggsidekick)

source('R/00-settings.R')

fit <- readRDS(file.path(fit_dir, paste0("dogfish-dgg", ".rds")))

prediction_grid0 <- readRDS(file.path(model_input_dir, "prediction-grid.rds"))
cell_polygons0 <- readRDS(file.path(model_input_dir, "wcvi-cell-polygons.rds"))

missing_cells <- anti_join(cell_polygons0, prediction_grid0)
cell_polygons <- cell_polygons0 |> filter(!(cell %in% missing_cells$cell))

pgrid <- prediction_grid0 |> filter(year %in% c(2005, 2022, 2042))
pgrid_2022 <- prediction_grid0 |> filter(year == 2022) |> tidyr::drop_na(bottom_temp)
pgrid_sf <- left_join(cell_polygons, pgrid)
pgrid_outline <- pgrid_sf |> sf::st_union() |> sf::st_sf()

dog <- sdmTMB::dogfish
dog_sf <- XY_to_sf(dog)

coast <- pacea::bc_coast

# Colours

# Base map
bmap <- ggplot() +
  geom_sf(data = coast, fill = "grey85", colour = "grey65", linewidth = 0.2) +
    theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  plot.title = element_text(size = 20, face = "bold"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12),
  legend.key = element_blank()
) +
  guides(colour = "none", size = "none", fill = "none")

# ------------------------------------------------------------------------------
# non-linear depth
nd <- fit$data |>
  distinct(X, Y, year, .keep_all = TRUE) |>
  mutate(bottom_temp = mean(pgrid_2022$bottom_temp))

pred_2022 <- predict(fit, newdata = nd, type = "link", nsim = 0, model = NA) |>
  mutate(est = fit$family[[1]]$linkinv(est1) * fit$family[[2]]$linkinv(est2))

ggplot(data = pred_2022 |> filter(exp(log_depth) < 600)) +
  aes(x = exp(log_depth), y = log(est), alpha = 0.5) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_minimal()

# ------------------------------------------------------------------------------
# prediction map
pred_2022 <- predict(fit, newdata = pgrid_2022, type = "link", nsim = 0, model = NA) |>
  mutate(est = fit$family[[1]]$linkinv(est1) * fit$family[[2]]$linkinv(est2))
pred_2022_sf <- left_join(cell_polygons, pred_2022) |>
  tidyr::drop_na(est)

dgg_pred_2022 <- bmap +
  geom_sf(data = pred_2022_sf, aes(fill = est)) +
  scale_fill_viridis_c(trans = "fourth_root_power") +
  gfplot::coord_sf_auto(pred_2022_sf, buffer = -100000)
dgg_pred_2022

# spatial field
fit_tw <- readRDS(file.path(fit_dir, "dogfish-tweedie.rds")) # use tweedie because dgg spatial RF collapsed
pred_tw <- predict(fit_tw, newdata = pgrid_2022, type = "link", nsim = 0, model = NA)
pred_tw_sf <- left_join(cell_polygons, pred_tw)

tw_omega_s <- bmap +
  geom_sf(data = pgrid_outline, linewidth = 0.8) +
  geom_sf(data = pred_tw_sf, aes(fill = omega_s, colour = omega_s)) +
  scale_colour_gradient2() +
  scale_fill_gradient2() +
  gfplot::coord_sf_auto(pred_2022_sf, buffer = -100000)

tw_epsilon_st <- bmap +
  geom_sf(data = pgrid_outline, linewidth = 0.8) +
  geom_sf(data = pred_tw_sf, aes(fill = epsilon_st, colour = epsilon_st)) +
  scale_colour_gradient2() +
  scale_fill_gradient2() +
  gfplot::coord_sf_auto(pred_2022_sf, buffer = -100000)


# ------------------------------------------------------------------------------
# Get plots for index using project
prediction_grid <- prediction_grid0 |> filter(year <= 2035)

fit_rw <-

projections_fe <- project(
  fit,
  newdata = prediction_grid,
  nsim = 200, # change to >100 for demo
  sample_parameters = FALSE, # if FALSE: theta at mode, b at mode → lp = lpb exactly, replicated. Zero stochasticity = single line.
  sample_historical_re = TRUE, # TRUE:
  sample_future_re = TRUE,
  sims_var = "eta_i",
  simulate_re = c(
      spatial = FALSE,
      spatiotemporal = TRUE,
      spatial_varying = FALSE,
      iid = FALSE,
      time_varying = TRUE,
      smoothers = FALSE
    ),
  model = NA
)
meep()

sims_fe <- make_sims_df(prediction_grid, projections_fe, last_historical_year) |>
  mutate(type = "fe")

a_by_sim <- sims |>
  group_by(time_period, region, year, sim) |>
  summarise(total_abundance = sum(exp(log_abundance)), .groups = "drop")

a_by_region <- a_by_sim |>
  summarise_sim_ci(cols = total_abundance, by = c(time_period, region, year))

a_by_region |>
  filter(region != "all") |>
ggplot() +
  aes(x = year, y = est, ymin = lwr, ymax = upr, colour = region, fill = region) +
  geom_line() +
  geom_ribbon(colour = NA, alpha = 0.3) +
  geom_vline(xintercept = last_historical_year, linetype = "dashed", colour = "grey30") +
  coord_cartesian(expand = FALSE)

a_by_sim |>
  filter(region == "all") |>
ggplot() +
  aes(x = year, y = log(total_abundance), colour = region, fill = region, group = interaction(sim, region)) +
geom_ribbon(data = a_by_region |> filter(region == "all"),
  aes(y = log(est), ymin = log(lwr), ymax = log(upr), group = NULL),
  fill = blues[[1]], colour = NA, alpha = 0.7) +
  geom_line(colour = blues[[4]], alpha = 0.3) +
  geom_line(data = a_by_region |> filter(region == "all"),
    aes(y = log(est), group = NULL),
    colour = blues[[5]], linewidth = 1) +
  geom_vline(xintercept = last_historical_year, linetype = "dashed", colour = "grey30") +
  coord_cartesian(expand = FALSE, clip = "on") +
  facet_wrap(~ region, ncol = 1) +
  ggtitle(unique(sims$type))


# ------------------------------------------------------------------------------
# save plots
fig_dir <- here::here('presentations', 'figs')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

panel_width <- 4.2   # in; 3 panels x 4.2in + gaps/margins fits a 13.33in-wide slide
panel_height <- 3.75 # in; keeps ~1.12 aspect ratio, leaves room for title/caption text boxes

save_panel <- function(plot, name) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot,
         width = panel_width, height = panel_height, dpi = 300, bg = "white")
}

save_panel(dgg_pred_2022, "model-pred-map-2022")
save_panel(tw_omega_s, "model-omega_s-map-2022")
save_panel(tw_epsilon_st, "model-epsilon_st-map-2022")


