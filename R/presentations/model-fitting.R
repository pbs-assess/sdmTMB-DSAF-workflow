library(dplyr)
library(ggplot2)
library(sf)
library(ggsidekick)
library(sdmTMB)

source('R/00-settings.R')
source(here::here("R", "functions.R"))

fit <- readRDS(file.path(fit_dir, paste0("dogfish-dgg", ".rds")))

prediction_grid0 <- readRDS(file.path(model_input_dir, "prediction-grid.rds"))
cell_polygons0 <- readRDS(file.path(model_input_dir, "wcvi-cell-polygons.rds"))

missing_cells <- anti_join(cell_polygons0, prediction_grid0)
cell_polygons <- cell_polygons0 |> filter(!(cell %in% missing_cells$cell))

prediction_grid <- prediction_grid0 |> filter(year <= 2040)

pgrid <- prediction_grid0 |> filter(year %in% c(2005, 2022, 2042))
pgrid_2022 <- prediction_grid0 |> filter(year == 2022) |> tidyr::drop_na(bottom_temp)
pgrid_sf <- left_join(cell_polygons, pgrid)
pgrid_outline <- pgrid_sf |> sf::st_union() |> sf::st_sf()

region_polys_path <- file.path(data_raw_dir, "wcvi-region-polygons", "wcvi-region-polygons.shp")
region_polys <- sf::st_read(region_polys_path, quiet = TRUE)

dog <- sdmTMB::dogfish
dog_sf <- gfplot::XY_to_sf(dog)

coast <- pacea::bc_coast

last_historical_year <- max(dog$year)

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
# sim_draws objects cached by R/03-get-shift-indicators-matrix.R
model_name <- "dogfish-dgg-rw"
# model_name <- "dogfish-dgg-ar1"
fit_path <- file.path(fit_dir, paste0(model_name, ".rds"))
fit <- readRDS(fit_path)
remind_model_info(fit)

last_historical_year <- max(fit$data$year) # use for centring lm trend intercept

set.seed(2026) # @Question - should we have set seed in project.sdmTMB?
proj <- project(
  fit,
  newdata = prediction_grid,
  nsim = 400,
  sample_fe = FALSE,
  sample_historical_re = TRUE,
  sample_future_re = TRUE,
  offset = rep(0, nrow(prediction_grid))
)

proj_mat <- return_density(proj, fit, log = TRUE)
sims_re <- make_sim_draws(prediction_grid, proj_mat, last_historical_year,
  type = "re", model_name = model_name)

# Map slices from sims_re: the per-cell median across draws for 2022 (closest
# analogue to a point prediction; the median of exp() equals exp() of the
# median) and a single draw for 2035. distinct() because cells split across
# regions carry identical log_density in each half.
map_years <- c(2022, 2035)
rows_map <- which(sims_re$grid$year %in% map_years)
slices <- sim_draw_slices(
  sims_re$grid[rows_map, ],
  sims_re$log_density[rows_map, , drop = FALSE],
  sims = 1, agg = list(median = matrixStats::rowMedians), value_name = "log_density"
) |>
  distinct(cell, year, draw, log_density)

pred_2022_sf <- cell_polygons |>
  left_join(filter(slices, year == 2022, draw == "median"), by = "cell") |>
  tidyr::drop_na(log_density)
sim1_2035_sf <- cell_polygons |>
  left_join(filter(slices, year == 2035, draw == "sim_1"), by = "cell") |>
  tidyr::drop_na(log_density)

# prediction map
dgg_pred_2022 <- bmap +
  geom_sf(data = pred_2022_sf, aes(fill = exp(log_density))) +
  scale_fill_viridis_c(trans = "fourth_root_power") +
  gfplot::coord_sf_auto(pred_2022_sf, buffer = -100000)
dgg_pred_2022

dgg_pred_2022 +
  geom_sf(data = region_polys, fill = NA, aes(colour = region), linewidth = 2) +
  scale_colour_brewer(palette = "Dark2") +
  gfplot::coord_sf_auto(pred_2022_sf)

#
dgg_proj_2035 <- bmap +
  geom_sf(data = sim1_2035_sf, aes(fill = exp(log_density))) +
  scale_fill_viridis_c(trans = "fourth_root_power") +
  gfplot::coord_sf_auto(pred_2022_sf, buffer = -100000)
dgg_proj_2035

# spatial field
fit_tw <- readRDS(file.path(fit_dir, "dogfish-tweedie.rds")) # use tweedie because dgg spatial RF collapsed
pred_tw <- predict(fit_tw, newdata = pgrid_2022, type = "link", nsim = 0, model = NA)
pred_tw_sf <- left_join(cell_polygons, pred_tw)

tw_omega_s <- bmap +
  geom_sf(data = pgrid_outline, linewidth = 0.8, colour = "grey80") +
  geom_sf(data = pred_tw_sf, aes(fill = omega_s, colour = omega_s)) +
  scale_colour_gradient2() +
  scale_fill_gradient2() +
  gfplot::coord_sf_auto(pred_2022_sf, buffer = -100000)

tw_epsilon_st <- bmap +
  geom_sf(data = pgrid_outline, linewidth = 0.8, colour = "grey80") +
  geom_sf(data = pred_tw_sf, aes(fill = epsilon_st, colour = epsilon_st)) +
  scale_colour_gradient2() +
  scale_fill_gradient2() +
  gfplot::coord_sf_auto(pred_2022_sf, buffer = -100000)

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
save_panel(dgg_proj_2035, "model-proj-map-2035")


