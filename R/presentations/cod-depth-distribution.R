library(dplyr)
library(ggplot2)

theme_set(gfplot::theme_pbs())

setdet0 <- readRDS(here::here("../NCAM", "analysis", "spatial_model", "derived", "atlantic-cod.rds"))

cod_depth <- setdet0 |>
  filter(number > 0 | weight > 0) |>
  transmute(log_depth = log10(depth))

central <- quantile(cod_depth$log_depth, probs = c(0.05, 0.95))

# equal-width log10 bins with edges landing exactly on the 10% and 90% quantiles,
# so the highlighted bins match the 80% range
n_inside <- 12
bin_width <- diff(central) / n_inside
n_below <- ceiling((central[[1]] - min(cod_depth$log_depth)) / bin_width)
n_above <- ceiling((max(cod_depth$log_depth) - central[[2]]) / bin_width)
bin_breaks <- central[[1]] + bin_width * (-n_below:(n_inside + n_above))

# strict/inclusive edges match stat_bin()'s default right-closed bins
cod_depth <- cod_depth |>
  mutate(in_central = log_depth > central[[1]] & log_depth <= central[[2]])

tallest_bin <- max(table(findInterval(cod_depth$log_depth, bin_breaks, left.open = TRUE)))
y_arrow <- tallest_bin * 1.12

range_label <- paste0(
  diff(as.numeric(sub("%", "", names(central)))), "% of cod occurrence\n",
  round(10^central[[1]]), " \u2013 ", round(10^central[[2]]), " m"
)

p_cod_depth <- ggplot(cod_depth) +
  aes(x = 10^log_depth, fill = in_central) +
  stat_bin(breaks = 10^bin_breaks) +
  geom_vline(xintercept = 10^central, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
  annotate(
    "segment",
    x = 10^central[[1]], xend = 10^central[[2]], y = y_arrow, yend = y_arrow,
    colour = "grey50", linewidth = 0.4,
    arrow = arrow(ends = "both", length = unit(0.06, "in"), type = "closed")
  ) +
  annotate(
    "text",
    x = 10^mean(central), y = y_arrow * 1.03, label = range_label,
    vjust = 0, lineheight = 0.9, size = 3.6, colour = "grey40"
  ) +
  scale_fill_manual(values = c(`TRUE` = "#2C7FB8", `FALSE` = "grey80"), guide = "none") +
  scale_x_continuous(
    trans = "log10", breaks = c(30, 100, 300, 1000),
    expand = expansion(mult = 0.01)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, y_arrow * 1.4)) +
  labs(x = "Depth (m)", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "grey40", linewidth = 0.3)
  )

p_cod_depth

# ~1/3 of a 7.5 in tall 16:9 slide
ggsave(
  here::here("presentations", "figs", "cod-depth-distribution.png"),
  p_cod_depth,
  width = 8, height = 2.5, dpi = 300
)
