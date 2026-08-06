library(dplyr)
library(ggplot2)
library(sf)

source('R/00-settings.R')

prediction_grid0 <- readRDS(file.path(model_input_dir, "prediction-grid.rds"))
cell_polygons <- readRDS(file.path(model_input_dir, "wcvi-cell-polygons.rds"))

dog <- sdmTMB::dogfish
dog_sf <- XY_to_sf(dog)

bt0 <- pacea::bccm_bottom_temperature()
bt <-  bt0 |> select(matches("_6$")) |>
  tidyr::pivot_longer(cols = -geometry, names_to = c("year", "month"), names_sep = "_",
values_to = "bottom_temp", names_transform = as.numeric) |>
  filter(year >= 2004)

pgrid <- prediction_grid0 |> filter(year %in% c(2005, 2022, 2042, 2060))
pgrid_sf <- left_join(cell_polygons, pgrid)

pgrid_2022 <- pgrid_sf |> filter(year == 2022)
pgrid_outline <- pgrid_2022 |> sf::st_union() |> sf::st_sf()

coast <- pacea::bc_coast
grid26_depth <- pacea::grid26_depth |>
  dplyr::mutate(depth_m = pmax(-mean_depth, 0)) |>
  gfplot::rotate_sf()


# Colours

depth_pal <- cmocean::cmocean("deep")(256)
ecdf_d <- ecdf(grid26_depth$depth_m)
grid26_depth$depth_rank <- ecdf_d(grid26_depth$depth_m) # show shelf contours better

brk_m <- c(10, 50, 200, 500, 1000, 3000)   # fewer breaks to avoid legend label crowding
brk_rank <- ecdf_d(brk_m)

temp_pal <- cmocean::cmocean("thermal")(256)

# ------------------------------------------------------------------------------

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

# Catch locations
catch_map <- bmap +
  geom_sf(
    data = dog_sf |> filter(year == 2022),
    aes(size = catch_weight),
    shape = 21, fill = "#043878", colour = "white",
    stroke = 0.3, alpha = 0.8
  ) +
  scale_size_area(trans = "sqrt", max_size = 7, breaks = scales::pretty_breaks(n = 3)) +
  gfplot::coord_sf_auto(pgrid_sf, buffer = -10000)
catch_map

catch_map_grid <- catch_map +
  geom_sf(data = pgrid_2022, fill = NA, colour = "grey30", linewidth = 0.08) +
  geom_sf(data = pgrid_outline, fill = NA, colour = "grey30", linewidth = 0.2) +
  gfplot::coord_sf_auto(pgrid_sf, buffer = -10000)
catch_map_grid

depth_map <- bmap +
  geom_sf(data = grid26_depth, aes(fill = depth_rank, colour = depth_rank)) +
  scale_fill_gradientn(colours = depth_pal, breaks = brk_rank, labels = brk_m, name = "Depth (m)") +
  scale_colour_gradientn(colours = depth_pal, breaks = brk_rank, labels = brk_m, name = "Depth (m)") +
  gfplot::coord_sf_auto(pgrid_sf, buffer = -10000)
depth_map

depth_map_grid <- depth_map +
  geom_sf(data = pgrid_2022, fill = NA, colour = "grey80", linewidth = 0.08) +
  geom_sf(data = pgrid_outline, fill = NA, colour = "white", linewidth = 0.2) +
  gfplot::coord_sf_auto(pgrid_sf, buffer = -10000)
depth_map_grid

temp_map <- bmap +
  geom_sf(data = bt |> filter(year == 2019), aes(fill = bottom_temp, colour = bottom_temp)) +
  scale_fill_gradientn(colours = temp_pal, limits = c(1.9, 11), oob = scales::squish,
                        breaks = scales::pretty_breaks(6), name = "Bottom temp (°C)") +
  scale_colour_gradientn(colours = temp_pal, limits = c(1.9, 11), oob = scales::squish,
                          breaks = scales::pretty_breaks(6), name = "Bottom temp (°C)") +
  gfplot::coord_sf_auto(pgrid_sf, buffer = -10000)
temp_map

temp_map_grid <- temp_map +
  geom_sf(data = pgrid_2022, fill = NA, colour = "grey80", linewidth = 0.08) +
  geom_sf(data = pgrid_outline, fill = NA, colour = "white", linewidth = 0.2) +
  gfplot::coord_sf_auto(pgrid_sf, buffer = -10000)
temp_map_grid

# ------------------------------------------------------------------------------
# Save for PowerPoint (widescreen 16:9, three panels in one row per slide).
# Panel aspect ratio (~1.12 w:h) matches the WCVI grid bbox. No legends, axes,
# or titles baked in -- add those in PPT. Base + *_grid pairs share identical
# framing so they can be layered as a duplicate-slide/animation reveal of the
# prediction grid.
fig_dir <- here::here('presentations', 'figs')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

panel_width <- 4.2   # in; 3 panels x 4.2in + gaps/margins fits a 13.33in-wide slide
panel_height <- 3.75 # in; keeps ~1.12 aspect ratio, leaves room for title/caption text boxes

save_panel <- function(plot, name) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot,
         width = panel_width, height = panel_height, dpi = 300, bg = "white")
}

save_panel(catch_map, "catch-map")
save_panel(catch_map_grid, "catch-map-grid")
save_panel(depth_map, "depth-map")
save_panel(depth_map_grid, "depth-map-grid")
save_panel(temp_map, "temp-map")
save_panel(temp_map_grid, "temp-map-grid")

# ------------------------------------------------------------------------------
# Thumbnails for regression-equation graphic: cropped (not squished) from
# the same full-extent panels above, trimming off mostly-empty ocean/land on
# each side, to use as small icons next to each term (depth, temperature,
# spatial) in a slide showing the model formula.
# coord_sf() locks x:y to equal map units, so to keep the same in-frame
# scale as the full panels (i.e. an actual crop, not a zoom), the physical
# output width/height must shrink by the same fraction as the retained map
# extent -- both are derived from map_scale (in/map-unit), calibrated once
# from the full bbox and panel_width.
# Note: gfplot::coord_sf_auto()'s buffer step is gated on `buffer > 0`, so the
# panels above (buffer = -10000) never actually buffer -- bbox is taken as-is.
full_bbox <- sf::st_bbox(pgrid_sf)
map_scale <- panel_width / unname(full_bbox["xmax"] - full_bbox["xmin"]) # in per map-unit

get_thumbnail_bbox <- function(sf_obj, height_frac = 1 / 3, left_trim_frac = 0.15, right_trim_frac = 0.15) {
  bbox <- sf::st_bbox(sf_obj)
  y_mid <- mean(bbox[c("ymin", "ymax")])
  y_half <- unname(bbox["ymax"] - bbox["ymin"]) * height_frac / 2
  x_range <- unname(bbox["xmax"] - bbox["xmin"])
  x_min <- unname(bbox["xmin"]) + x_range * left_trim_frac
  x_max <- unname(bbox["xmax"]) - x_range * right_trim_frac
  list(
    xlim = c(x_min, x_max),
    ylim = c(y_mid - y_half, y_mid + y_half)
  )
}

thumb_bbox <- get_thumbnail_bbox(pgrid_sf)
thumb_xlim <- thumb_bbox$xlim
thumb_ylim <- thumb_bbox$ylim
thumb_width <- map_scale * diff(thumb_xlim)
thumb_height <- map_scale * diff(thumb_ylim)

catch_map_thumb <- catch_map +
  gfplot::coord_sf_auto(pgrid_sf, buffer = -10000, xlim = thumb_xlim, ylim = thumb_ylim)
depth_map_thumb <- depth_map +
  gfplot::coord_sf_auto(pgrid_sf, buffer = -10000, xlim = thumb_xlim, ylim = thumb_ylim)
temp_map_thumb <- temp_map +
  gfplot::coord_sf_auto(pgrid_sf, buffer = -10000, xlim = thumb_xlim, ylim = thumb_ylim)

save_thumb <- function(plot, name) {
  ggsave(file.path(fig_dir, paste0(name, "-thumb.png")), plot,
         width = thumb_width, height = thumb_height, dpi = 300, bg = "white")
}

save_thumb(catch_map_thumb, "catch-map")
save_thumb(depth_map_thumb, "depth-map")
save_thumb(temp_map_thumb, "temp-map")
