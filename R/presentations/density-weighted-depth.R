mean_depth_by_sim <- sims_re |>
  filter(region != "all") |>
  group_by(year, region, sim) |>
  summarise(
    mean_depth = weighted.mean(
      depth_m,
      w = exp(log_abundance),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

mean_depth_summary <- mean_depth_by_sim |>
  group_by(year, region) |>
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
  aes(x = year, y = mean, colour = region, fill = region, group = region)
) +
  geom_ribbon(
    aes(ymin = q10, ymax = q90),
    # fill = blues[[2]],
    colour = NA,
    alpha = 0.35
  ) +
  geom_line(
    data = mean_depth_examples,
    aes(y = mean_depth, group = interaction(sim, region)),
    # colour = blues[[3]],
    alpha = 0.25,
    linewidth = 0.35
  ) +
  geom_line(linewidth = 1) +
  geom_vline(
    xintercept = last_historical_year,
    linetype = "dashed",
    colour = "grey40"
  ) +
  scale_y_reverse() +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    x = "Year",
    y = "Density-weighted mean depth (m)",
    # subtitle = "Mean, 80% simulation interval, and 40 example draws",
    colour = "Region",
    fill = "Region"
  ) +
  coord_cartesian(expand = FALSE, clip = "off") +
  theme(
    legend.position = "top",
    plot.margin = margin(5.5, 16, 5.5, 5.5)
  )

p_mean_depth

ggsave(file.path(fig_dir, "density-weighted-depth-dogfish-eg.png"),
  width = 5.6, height = 3.75
)
