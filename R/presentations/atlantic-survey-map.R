library(marmap)
library(rnaturalearth)

theme_set(gfplot::theme_pbs())

d0 <- read.csv(here::here("data-raw", "ATLANTIC_HALIBUT_Joined_05.27.2026.csv"))

glimpse(d0)

test <- d0 |> filter((latitude < 42 & longitude > -60 & longitude < -50))
glimpse(test)
test |>  distinct(year)

# Make a field for nafo zone groupings, for indexing later
d <- d0 |>
  filter(#year >= 1993,
         !NAFO_div %in% c("3Pn", "4R", "6B", "6C", "4T", "6F",
                          "0B",  "1C",  "1D",  "1B",  "1A",  "0A",
                          NA, "2G" )) |> #Dropping 3Pn(too much mixing with GSL) and 4R (GSL) , 6BC are too south
  mutate(
    NAFO_ZONES = case_when(
      NAFO_div %in% c("2J","2H") ~ "2HJ", #Labrador
      NAFO_div %in% c("3K", "3L", "3M") ~ "3KLM", # Newfoundland NorthEast
      NAFO_div %in% c("3Ps", "3N", "3O") ~ "3NOPs", #Newfoundland South
      NAFO_div %in% c("4Vn","4W","4Vs") ~ "4VW",  #Western Scotian Shelf
      NAFO_div == "4X" ~ "4X", #Eastern Scotian Shelf
      NAFO_div %in% c("5Y","5Ze", "5Zw", "6A") ~ "5YZ6A" ) #GUlf of Maine
  ) |>
  rename(fe_id = X) |> # Not sure what the X column is
  sdmTMB::add_utm_columns(ll_names = c("longitude", "latitude"))

# drop the isolated offshore outlier cluster (~39-41N, -58 to -50W)
weird <- d |>
  filter((latitude < 42 & longitude > -60 & longitude < -50))
d <- d0 |>
  filter(!(latitude < 42 & longitude > -60 & longitude < -50))


modal_utm_zone <- sdmTMB::get_crs(d)

d_sf <- st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326)
weird_sf <- st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326)

# bounding box with a small buffer, from the data itself
buf <- 1
lon_range <- range(d$longitude, na.rm = TRUE) + c(-buf, buf)
lat_range <- range(d$latitude, na.rm = TRUE) + c(-buf, buf)

# coastline (high-res land polygons)
coast <- rnaturalearth::ne_countries(scale = "large", continent = "north america", returnclass = "sf")

# bathymetry (downloads NOAA ETOPO once, caches locally as a .csv)
bathy <- getNOAA.bathy(
  lon1 = lon_range[1], lon2 = lon_range[2],
  lat1 = lat_range[1], lat2 = lat_range[2],
  resolution = 4, keep = TRUE
)
bathy_df <- fortify.bathy(bathy)

# depth values at the d_sf tow locations
d_coords <- st_coordinates(d_sf)
d_depths <- get.depth(bathy, x = d_coords, locator = FALSE)
d_sf$depth <- d_depths$depth

d_sf$depth |> range()
d_sf$year |> range()

ggplot(d_sf) +
  geom_histogram(aes(x = -depth)) +
  scale_x_continuous(
    trans = "log10",
    breaks = c(10, 50, 100, 500, 1000, 2000),
    labels = scales::label_comma(),
    name = "Depth (m)"
  ) +
  gfplot::theme_pbs()

survey_footprint <- d_sf |>
  filter(year >= 2000) |>
  sf::st_transform(modal_utm_zone) |>
  group_by(survey) |>
  sf::st_union() |>
  sf::st_concave_hull(ratio = 0.1, allow_holes = FALSE) |>
  st_buffer(10)

survey_footprint <- st_difference(survey_footprint, coast |> st_transform(modal_utm_zone))

# setdet <- st_intersection(all_utm_coords, survey_footprint)

ggplot() +
  geom_raster(data = bathy_df, aes(x = x, y = y, fill = z)) +
  scale_fill_gradientn(colors = c("grey30", "grey90"), guide = "none") +
  geom_contour(data = bathy_df, aes(x = x, y = y, z = z),
               breaks = c(-100, -200, -500, -1000, -2000),
               color = "grey60", linewidth = 0.2) +
  geom_sf(data = coast, fill = "grey40", color = NA) +
  geom_sf(data = d_sf, size = 0.8, alpha = 0.3, aes(colour = survey)) +
  # geom_sf(data = survey_footprint, aes(colour = survey), fill = NA, linewidth = 0.6) +
  scale_colour_brewer(palette = "Dark2") +
  # geom_sf(data = d_sf, size = 0.5, alpha = 0.4, color = "firebrick") +
  coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
  theme_minimal() #+ facet_wrap(~ year)
# I think the spatiotemporal overlap looks a bit misleading if not faceted by year

# Note: for presentation I just put this ^ plot on my vertical monitor full screen and took a screenshot

glimpse(d_sf)

distinct(d, survey)

# ggplot() +
#   # geom_raster(data = bathy_df, aes(x = x, y = y, fill = z)) +
#   scale_fill_gradientn(colors = c("grey30", "grey90"), guide = "none") +
#   geom_contour(data = bathy_df, aes(x = x, y = y, z = z),
#                breaks = c(-100, -200, -500, -1000, -2000),
#                color = "grey60", linewidth = 0.2) +
#   geom_sf(data = coast, fill = "grey40", color = NA) +
#   geom_sf(data = d_sf, size = 0.5, alpha = 0.4, aes(colour = depth)) +
#   coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
#   facet_wrap(~ year)

# number of sets per year, by survey
sets_per_year <- d_sf |>
  st_drop_geometry() |>
  distinct(survey, year, set = trawl_id) |>
  dplyr::count(survey, year, name = "n_sets")

ggplot(sets_per_year, aes(x = year, y = n_sets, colour = survey)) +
  geom_line() +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Year", y = "Number of sets", colour = "Survey")

sets_per_month <- d_sf %>%
  st_drop_geometry() %>%
  mutate(month = lubridate::month(date, label = TRUE)) %>%
  distinct(survey, month, set = trawl_id) %>%
  dplyr::count(survey, month, name = "n_sets") %>%
  tidyr::complete(survey, month, fill = list(n_sets = 0))

ggplot(sets_per_month, aes(x = survey, y = month, fill = n_sets)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_sets), size = 3) +
  scale_y_discrete(limits = rev) +
  scale_fill_gradientn(colors = c("grey95", "steelblue4"), guide = "none") +
  labs(x = "Survey", y = "Month", title = "Number of sets by month and survey") +
  gfplot::theme_pbs()

library(patchwork)

sets_per_year_month <- d_sf %>%
  st_drop_geometry() %>%
  mutate(month = lubridate::month(date, label = TRUE)) %>%
  distinct(survey, year, month, set = trawl_id) %>%
  dplyr::count(survey, year, month, name = "n_sets")

# one panel per survey: each keeps only the months it actually samples (free y)
# and gets its own colour scale (free fill), since raw magnitudes differ a lot
# by survey — stacked with patchwork so the shared x axis (year) still lines up
survey_year_month_panel <- function(dat, surv) {
  dat_surv <- dat %>%
    filter(survey == surv) %>%
    tidyr::complete(year = min(year):max(year), month, fill = list(n_sets = 0))

  active_months <- dat_surv %>%
    filter(n_sets > 0) %>%
    distinct(month) %>%
    pull(month)

  dat_surv %>%
    filter(month %in% active_months) %>%
    ggplot(aes(x = year, y = month, fill = n_sets)) +
    geom_tile(color = "white") +
    scale_y_discrete(limits = rev, drop = TRUE) +
    scale_fill_gradientn(colors = c("grey95", "steelblue4"), guide = "none") +
    labs(x = NULL, y = NULL, subtitle = surv) +
    gfplot::theme_pbs()
}

panels <- lapply(sort(unique(sets_per_year_month$survey)), \(s) survey_year_month_panel(sets_per_year_month, s))

wrap_plots(panels, ncol = 1) +
  plot_annotation(title = "Number of sets by year and month, by survey") &
  scale_x_continuous(limits = range(sets_per_year_month$year))

# ---- spatial overlap between surveys, within year ----
# co-occurring in the same year isn't enough for a survey factor to be
# identifiable separately from the spatiotemporal field -- the surveys also
# need to actually sample overlapping locations that year. Grid the study
# area (equal-area projection, cell size on the same order as the sdmTMB
# mesh cutoff used in Halibut_Model.R) and, for each year, measure what
# fraction of the smaller survey's sampled cells are shared with the other
# survey.
overlap_cell_km <- 25
aea_crs <- "+proj=aea +lat_1=41 +lat_2=52 +lat_0=46 +lon_0=-66 +datum=NAD83 +units=km +no_defs"

d_aea <- st_transform(d_sf, aea_crs)

grid_sfc <- st_make_grid(d_aea, cellsize = overlap_cell_km)
grid <- st_sf(cell_id = seq_along(grid_sfc), geometry = grid_sfc)

set_cells <- st_join(d_aea, grid, join = st_intersects) %>%
  st_drop_geometry() %>%
  distinct(survey, year, cell_id) %>%
  filter(!is.na(cell_id))

year_survey_cells <- set_cells %>%
  group_by(year, survey) %>%
  summarise(cells = list(unique(cell_id)), .groups = "drop")

# for one year's worth of surveys, overlap_coef = |shared cells| / |smaller
# survey's cells| -- asymmetric on purpose, since a small survey fully
# nested inside a larger one's footprint is still identifiable
survey_pair_overlap <- function(dat_year) {
  if (nrow(dat_year) < 2) {
    return(tibble(survey_a = character(), survey_b = character(),
                   n_shared = integer(), overlap_coef = double()))
  }
  combn(seq_len(nrow(dat_year)), 2, simplify = FALSE) %>%
    purrr::map_dfr(function(idx) {
      a <- dat_year[idx[1], ]
      b <- dat_year[idx[2], ]
      shared <- length(intersect(a$cells[[1]], b$cells[[1]]))
      tibble(
        survey_a = a$survey, survey_b = b$survey,
        n_shared = shared,
        overlap_coef = shared / min(length(a$cells[[1]]), length(b$cells[[1]]))
      )
    })
}

overlap_by_year <- year_survey_cells %>%
  tidyr::nest(data = -year) %>%
  mutate(overlap = purrr::map(data, survey_pair_overlap)) %>%
  select(year, overlap) %>%
  tidyr::unnest(overlap) %>%
  mutate(pair = paste(survey_a, survey_b, sep = " × "))

all_pairs <- combn(sort(unique(d_sf$survey)), 2, simplify = FALSE) %>%
  purrr::map_chr(~ paste(.x, collapse = " × "))

# every year x pair combo, so years/pairs that never co-occur show up blank
# rather than just being absent from the plot
overlap_grid <- tidyr::expand_grid(
  year = sort(unique(d_sf$year)),
  pair = all_pairs
) %>%
  left_join(overlap_by_year %>% select(year, pair, overlap_coef), by = c("year", "pair"))

ggplot(overlap_grid, aes(x = year, y = pair, fill = overlap_coef)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colors = c("grey90", "steelblue4"),
    na.value = "white",
    limits = c(0, 1),
    name = "Overlap"
  ) +
  labs(x = "Year", y = NULL,
       title = "Spatial overlap between survey pairs, by year",
       subtitle = paste0("Cell size ", overlap_cell_km, " km; white = surveys not co-occurring that year, ",
                          "grey = co-occurring but zero spatial overlap")) +
  gfplot::theme_pbs()

# ---- checking for genuine co-located calibration tows across surveys ----
# distinct from the coarse grid-cell overlap above: this looks tow-by-tow for
# pairs from different surveys that are suspiciously close in BOTH space and
# time (i.e. a formal calibration study), which would be a much stronger
# identifiability signal than "these surveys operated the same year".
#
# NOTE: one DFO_NF row (trawl_id 1999-39-084-040) has date = 3199-10-20,
# a typo'd year (the `year` column correctly says 1999) -- likely an error
# upstream of this script, not introduced here. Dropped for this check via
# the year(date) == year filter below; worth flagging/fixing at the source.
tows <- d_aea %>%
  group_by(survey, year, trawl_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(survey, year, trawl_id, date) %>%
  mutate(date = as.Date(date)) %>%
  filter(lubridate::year(date) == year)

nearest_pair_distances <- function(year_val, survey_a, survey_b, tows) {
  pts_a <- tows %>% filter(year == year_val, survey == survey_a)
  pts_b <- tows %>% filter(year == year_val, survey == survey_b)
  if (nrow(pts_a) == 0 || nrow(pts_b) == 0) return(NULL)
  nn <- st_nearest_feature(pts_a, pts_b)
  dist_km <- as.numeric(st_distance(pts_a, pts_b[nn, ], by_element = TRUE))
  day_gap <- abs(as.numeric(difftime(pts_a$date, pts_b$date[nn], units = "days")))
  tibble(
    year = year_val, survey_a = survey_a, survey_b = survey_b,
    trawl_id_a = pts_a$trawl_id, trawl_id_b = pts_b$trawl_id[nn],
    dist_km = dist_km, day_gap = day_gap
  )
}

survey_pairs <- combn(sort(unique(tows$survey)), 2, simplify = FALSE)

nearest_pairs <- purrr::map_dfr(survey_pairs, function(pr) {
  purrr::map_dfr(sort(unique(tows$year)), function(y) nearest_pair_distances(y, pr[1], pr[2], tows))
})

ggplot(nearest_pairs, aes(x = day_gap, y = dist_km)) +
  geom_point(alpha = 0.15, size = 0.6) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 5, 10, 25, 50, 100, 250)) +
  facet_wrap(~ paste(survey_a, survey_b, sep = " x ")) +
  labs(x = "Day gap to nearest different-survey tow",
       y = "Distance to nearest different-survey tow (km, log1p)",
       title = "Are there candidate calibration tow pairs?",
       subtitle = "A cluster near (0, 0) would indicate genuine paired/calibration tows -- none found") +
  gfplot::theme_pbs()
