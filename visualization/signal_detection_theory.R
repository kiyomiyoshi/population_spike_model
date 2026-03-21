library(tidyverse)
library(cowplot)
library(patchwork)

x_vals <- seq(-6, 6, length.out = 1000)
df <- tibble(
  x = x_vals,
  density_1 = dnorm(x_vals, mean = -0.75, sd = 1),
  density_2 = dnorm(x_vals, mean =  0.75, sd = 1)
) %>%
  pivot_longer(cols = starts_with("density"),
               names_to = "distribution", values_to = "density")

labels <- c("―CW―", "―CCW―") 
x_pos  <- c(-0.60, 1.38)

# color_1 <- "#E41A1C"
# color_2 <- "#377EB8"

p1 <- ggplot(df, aes(x = x, y = density, color = distribution)) +
  geom_segment(x =  0,   xend =  0,   y = 0, yend = 0.42, linetype = "solid",  color = "black",  size = 0.6) +
  geom_segment(x = -1, xend = -1, y = 0, yend = 0.42, linetype = "dashed", color = "grey70", size = 0.5) +
  geom_segment(x = -2, xend = -2, y = 0, yend = 0.42, linetype = "dashed", color = "grey70", size = 0.5) +
  geom_segment(x =  1, xend =  1, y = 0, yend = 0.42, linetype = "dashed", color = "grey70", size = 0.5) +
  geom_segment(x =  2, xend =  2, y = 0, yend = 0.42, linetype = "dashed", color = "grey70", size = 0.5) +
  geom_line(size = 0.7, alpha = 1) +
  annotate(
    "text",
    fontface = "bold",
    x = x_pos,
    y = 1,
    label = labels) +
# scale_color_manual(values = c("density_1" = color_1, "density_2" = color_2)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2), limits = c(-6, 6)) +
  scale_y_continuous(limits = c(0, 0.42)) +
  labs(x = "Signal strength", y = "Density") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

g1 <- list(p1)
g1 <- lapply(g1, function(p) {
  p + theme(
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )
})

sdt <- wrap_plots(g1, nrow = 1)
ggsave("sdt.png", sdt, width = 2.1, height = 2, dpi = 300)