source(here::here('R/presentations/enact-privacy-filter.R'))
trawl_grid <- gfdata::load_survey_blocks(dataset = "syn_hbll") |>
  filter(grepl("SYN", survey_abbrev))

hbll_grid <- gfdata::load_survey_blocks(dataset = "syn_hbll") |>
  filter(grepl("HBLL", survey_abbrev))

# comm_dat_pcod0 <- readRDS(here::here("../../gfsynopsis-2025-data/report/data-cache-2026-07/pacific-cod.rds"))$catch_spatial
# comm_dat_dog0 <- readRDS(here::here("../../gfsynopsis-2025-data/report/data-cache-2026-07/north-pacific-spiny-dogfish.rds"))$catch_spatial
# comm_dat_petrale0 <- readRDS(here::here("../../gfsynopsis-2025-data/report/data-cache-2026-07/petrale-sole.rds"))$catch_spatial
# comm_dat_shortspine0 <- readRDS(here::here("../../gfsynopsis-2025-data/report/data-cache-2026-07/shortspine-thornyhead.rds"))$catch_spatial
# comm_dat_ratfish0 <- readRDS(here::here("../../gfsynopsis-2025-data/report/data-cache-2026-07/spotted-ratfish.rds"))$catch_spatial
# comm_dat_arrowtooth0 <- readRDS(here::here("../../gfsynopsis-2025-data/report/data-cache-2026-07/arrowtooth-flounder.rds"))$catch_spatial
# comm_dat <- bind_rows(comm_dat_dog0, comm_dat_pcod0, comm_dat_petrale0,
#   comm_dat_shortspine0, comm_dat_ratfish0, comm_dat_arrowtooth0) |>
#   filter(year >= 2013) |>
#   sf::st_as_sf(coords = c("lon", "lat"), crs = 'WGS84') |>
#   sf::st_transform(crs = 32609)


data_dir <- here::here("../../gfsynopsis-2025-data/report/data-cache-2026-07")

comm_dat <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE) |>
  purrr::set_names(~ sub("\\.rds$", "", basename(.x))) |>
  purrr::map_dfr(~ readRDS(.x)$catch_spatial) |>
  filter(year >= 2013) |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 'WGS84') |>
  sf::st_transform(crs = 32609)

comm_list <- get_public_grid(dat = comm_dat, cell_size = 2000, .groups = "grid_id")
percent_omitted <- comm_list$omitted_df |>
  mutate(percent = percent * 100) |>
  filter(omitted == 1) |>
  pull(percent)

public_grid <- comm_list$public_grid_df

idx <- !duplicated(sf::st_drop_geometry(public_grid)["grid_id"])
public_grid_distinct <- public_grid[idx, ]

survey_cells <- trawl_grid
comm_cells <- public_grid_distinct
coast <- pacea::bc_coast

# Commercial cells are 2 km; classify by whether the cell centre falls in a survey block
comm_cells <- comm_cells |>
  mutate(covered = lengths(sf::st_intersects(sf::st_centroid(sf::st_geometry(comm_cells)), survey_cells)) > 0)

pct_cells_covered <- round(100 * mean(comm_cells$covered))
pct_events_covered <- round(100 * sum(comm_cells$count[comm_cells$covered]) / sum(comm_cells$count))

fill_labels <- c(
  survey = "Survey grid",
  covered = glue::glue(
    "Commercial, inside grid\n{pct_cells_covered}% of cells · {pct_events_covered}% of events"
  ),
  uncovered = glue::glue(
    "Commercial, outside grid\n{100 - pct_cells_covered}% of cells · {100 - pct_events_covered}% of events"
  )
)

legend_note <- glue::glue(
  "Darker = more events per cell\n{round(percent_omitted, 1)}% of events omitted for privacy"
)

p_overlap <- ggplot() +
  geom_sf(aes(fill = "survey"), data = survey_cells, colour = NA) +
  geom_sf(aes(fill = "covered", alpha = count), data = filter(comm_cells, covered), colour = NA) +
  geom_sf(aes(fill = "uncovered", alpha = count), data = filter(comm_cells, !covered), colour = NA) +
  geom_sf(data = coast, fill = "grey90", colour = "grey55", linewidth = 0.15) +
  scale_fill_manual(
    values = c(survey = "#C9D3DD", covered = "#2F4B7C", uncovered = "#E8761E"),
    labels = fill_labels,
    breaks = names(fill_labels),
    name = legend_note
  ) +
  scale_alpha_continuous(range = c(0.35, 1), transform = "log10", guide = "none") +
  guides(fill = guide_legend(override.aes = list(colour = NA, alpha = 1))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  gfplot::coord_sf_auto(survey_cells) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.02),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = alpha("white", 0.9), colour = NA),
    legend.margin = margin(6, 8, 6, 6),
    legend.title.position = "bottom",
    legend.title = element_text(size = 7, face = "plain", colour = "grey35", hjust = 0, lineheight = 1.05, margin = margin(t = 6)),
    legend.text = element_text(size = 8, lineheight = 1.05, colour = "grey15"),
    # uniform swatches, gap between entries, labels vertically centred on their swatch
    legend.key.width = unit(0.45, "cm"),
    legend.key.height = unit(0.85, "cm"),
    legend.key.spacing.y = unit(3, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    axis.title = element_blank(),
    axis.text = element_text(size = 8)
  )

p_overlap

# 4.5 x 5.2 in fits roughly 1/3-1/2 of a 16:9 Google Slide (10 x 5.625 in)
ggsave(
  here::here("presentations", "figs", "bc-commercial-overlap.png"),
  p_overlap, width = 4.5, height = 5.2, dpi = 300
)
beepr::beep()
