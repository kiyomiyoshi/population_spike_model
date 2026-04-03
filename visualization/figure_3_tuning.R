library(tidyverse)
library(patchwork)
library(ggpp)

# tuning curves
Rmax <- 115
C50 <- 19.3
n   <-  2.9

color_wheel <- function(alpha = 1) {
  h <- c(seq(160, 0, -1), seq(180, 161, -1))[-161]
  h <- h / 180
  s <- rep(1, 180)
  v <- rep(1, 180)
  hsv(h, s, v, alpha)
}

colors <- color_wheel()
n_neurons <- 180
orientations <- seq(1, 180, by = 1)
preferred_orientations <-  seq(1, 180, length.out = n_neurons)
max_firing_rate <- Rmax * 25^n / (25^n + C50^n)
tuning_width <- 20

tuning_curves <- matrix(0, nrow = n_neurons, ncol = length(orientations))
for (i in 1:n_neurons) {
  tuning_curves[i, ] <- max_firing_rate * exp(-0.5 * (pmin(abs(orientations - preferred_orientations[i]), 
                                                           180 - abs(orientations - preferred_orientations[i])) ^ 2) / tuning_width ^ 2)
}

FiringRate <- c()
for (i in 1:n_neurons) {
  FiringRate <- c(FiringRate, tuning_curves[i, ])
}

df_tuning <- data.frame(Orientation = rep(orientations, n_neurons),
                        FiringRate = FiringRate,
                        Neuron = rep(rep(orientations, n_neurons), each = length(orientations)),
                        Color = rep(colors, each = length(orientations)))

g4 <- ggplot(df_tuning[df_tuning$Neuron %in% seq(1, n_neurons, by = 18), ], 
             aes(x = Orientation, y = FiringRate, color = Color)) +
  geom_line() + scale_color_identity() +
  annotate("text_npc",
           npcx = 0.5,
           npcy = 1,
           label = "High α", 
           size = 3.5) +
  theme_classic(base_size = 10) + 
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  scale_y_continuous(lim = c(0, 120), breaks = c(0, 50, 100)) +
  guides(color = F) + ylab("Spikes")
g4

df_tuning$FiringRate <- df_tuning$FiringRate + 20
g5 <- ggplot(df_tuning[df_tuning$Neuron %in% seq(1, n_neurons, by = 18), ], 
             aes(x = Orientation, y = FiringRate, color = Color)) +
  geom_line() + scale_color_identity() +
  annotate("text_npc",
           npcx = 0.5,
           npcy = 1,
           label = "Low α", 
           size = 3.5) +
  annotate("text_npc",
           npcx = 0.5,
           npcy = 0.88,
           label = "increased baseline response", 
           size = 2.5) +
  theme_classic(base_size = 10) + 
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  scale_y_continuous(lim = c(0, 120), breaks = c(0, 50, 100)) +
  guides(color = F) + ylab("Spikes")
g5

f1 <- list(g5, g4)
f1 <- lapply(f1, function(p) {
  p + theme(
    legend.position = "none",
    plot.margin = margin(t = 5, r = 0, b = 0, l = 0, unit = "pt")
  )
})

figure_3_tuning <- wrap_plots(f1, nrow = 2)
ggsave("figure_3_tuning.png", figure_3_tuning, width = 2, height = 3.2, dpi = 300)