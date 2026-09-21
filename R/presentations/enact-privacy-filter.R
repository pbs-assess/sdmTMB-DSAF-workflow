# Enact privacy filter
get_public_grid <- function(dat = st_clean_df, cell_size = 5000,
                            .groups = c('grid_id')) {
  group_syms <- syms(.groups)  # can choose to ignore 'bin' group

  if ('cpue' %in% names(dat)) {
    metric <- 'cpue'
  }
  if ('catch' %in% names(dat)) {
    metric <- 'catch'
  }

  # Make grid based on cell size
  pgrid <- dat |>
  # pgrid <- p |>
  st_make_grid(cellsize = c(cell_size, cell_size), square = TRUE) |>
  sf::st_as_sf() |>
  dplyr::rename(geometry = x) |>
  mutate(grid_id = row_number(),
         cell_size = cell_size)

  intersect_df <- st_intersection(dat, pgrid)

  # Count the number of points from each boat within each grid cell
  # counts <- intersect_df |>
  #   group_by(grid_id, cell_size) |>
  #   summarise(count = n(), .groups = "drop") |>
  #   mutate(binning = "all time")
  privacy_by_group <- function(data, group_syms) {
    data |>
      group_by(!!!group_syms) |>
      summarise(count = n(), .groups = "drop") |>
      mutate(privacy_bins := paste(.groups, collapse = "-"))
  }
  counts <- privacy_by_group(intersect_df, group_syms)
  counts_lu <- counts |>
    st_drop_geometry()

  # Get grid cell ids that can be public
  # public_cells <- filter(counts, count >= 3) |>
  #   pull('grid_id')
  public_df <- intersect_df |>
    left_join(counts_lu) |>
    filter(count >= 3) |>
    group_by(!!!group_syms) |>
    #mutate(sum_cpue = sum(cpue)) |>
    mutate(!!metric := ifelse(metric == 'cpue', exp(mean(log(cpue), na.rm = FALSE)),
      sum(catch / 1000, na.rm = FALSE))) |>
    #mutate(cpue = exp(mean(log(cpue), na.rm = FALSE))) |>
    #mutate(catch = sum(catch / 1000, na.rm = FALSE)) |>
    ungroup()

  omitted_df <- intersect_df |>
    left_join(counts_lu) |>
    mutate(omitted = ifelse(count >= 3, 0, 1)) |>
    select(grid_id, count, omitted) |>
    st_drop_geometry() |>
    as_tibble()
  if (nrow(omitted_df) == 0L) {
    omitted_df <- NA
  } else {
    omitted_df <- omitted_df |>
      janitor::tabyl(omitted) |>
      janitor::adorn_totals("row")
      # janitor::adorn_percentages("row") |>
      # janitor::adorn_pct_formatting() |>
  }

  public_grid_df <- public_df |>
    st_drop_geometry() |>
    left_join(pgrid) |>
    st_sf(geometry = _)

  out <- list(pgrid = pgrid, public_grid_df = public_grid_df, omitted_df = omitted_df)
}