library(tidyverse)
library(cowplot)
library(patchwork)
library(ggpp)

### Naka-Rushton function ###
Contrast <- seq(0, 10, 0.1)
Rmax <- 115
C50 <- 19.3
n   <-  2.9
Max_firing <- Rmax * Contrast^n / (Contrast^n + C50^n)
nr <- data.frame(Contrast, Max_firing)

points_df <- nr %>%
  dplyr::filter(Contrast %in% c(3, 4, 5, 6))

g1 <- ggplot(nr, aes(x = Contrast, y = Max_firing)) +
  geom_line(alpha = 0.85, linewidth = 0.7) +
  geom_point(
    data = points_df,
    size = 2
  ) +
  scale_x_continuous(
    limits = c(1, 6),
    breaks = c(0, 2, 4, 6)
  ) +
  coord_cartesian(ylim = c(0, 5)) +
  scale_y_continuous(breaks = seq(0, 5, by = 1)) +
  labs(
    x = "Contrast (%)",
    y = "Tuning curve peak"
  ) +
  theme_classic(base_size = 11)

g1

### Population spike model ###
color_wheel <- function(alpha = 1) {
  h <- c(seq(160, 0, -1), seq(180, 161, -1))[-161]
  h <- h / 180
  s <- rep(1, 180)
  v <- rep(1, 180)
  hsv(h, s, v, alpha)
}

colors <- color_wheel()

# parameters
n_neurons <- 180                                              
orientations <- seq(1, 180, by = 1)                           
preferred_orientations <-  seq(1, 180, length.out = n_neurons)
tuning_width    <- 75

# contrast = 3%
contrast <- 3
max_firing_rate <- Rmax * contrast^n / (contrast^n + C50^n)
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

g2 <- ggplot(df_tuning[df_tuning$Neuron %in% seq(1, n_neurons, by = 18), ], 
             aes(x = Orientation, y = FiringRate, color = Color)) +
  geom_line() + scale_color_identity() +
  annotate("text_npc",
           npcx = 0.5,
           npcy = 1,
           label = "Contrast = 3%", 
           size = 3.5) +
  theme_classic(base_size = 11) + 
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  coord_cartesian(ylim = c(0, 5)) +
  scale_y_continuous(breaks = seq(0, 5, by = 1)) +
  guides(color = F) + ylab("Spikes")
g2

# contrast = 4%
contrast <- 4
max_firing_rate <- Rmax * contrast^n / (contrast^n + C50^n)
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

g3 <- ggplot(df_tuning[df_tuning$Neuron %in% seq(1, n_neurons, by = 18), ], 
             aes(x = Orientation, y = FiringRate, color = Color)) +
  geom_line() + scale_color_identity() +
  annotate("text_npc",
           npcx = 0.5,
           npcy = 1,
           label = "Contrast = 4%", 
           size = 3.5) +
  theme_classic(base_size = 11) + 
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  coord_cartesian(ylim = c(0, 5)) +
  scale_y_continuous(breaks = seq(0, 5, by = 1)) +
  guides(color = F) + ylab("Spikes")
g3

# contrast = 5%
contrast <- 5
max_firing_rate <- Rmax * contrast^n / (contrast^n + C50^n)
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
           label = "Contrast = 5%", 
           size = 3.5) +
  theme_classic(base_size = 11) + 
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  coord_cartesian(ylim = c(0, 5)) +
  scale_y_continuous(breaks = seq(0, 5, by = 1)) +
  guides(color = F) + ylab("Spikes")
g4

# contrast = 6%
contrast <- 6
max_firing_rate <- Rmax * contrast^n / (contrast^n + C50^n)
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

g5 <- ggplot(df_tuning[df_tuning$Neuron %in% seq(1, n_neurons, by = 18), ], 
             aes(x = Orientation, y = FiringRate, color = Color)) +
  geom_line() + scale_color_identity() +
  annotate("text_npc",
           npcx = 0.5,
           npcy = 1,
           label = "Contrast = 6%", 
           size = 3.5) +
  theme_classic(base_size = 11) + 
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  coord_cartesian(ylim = c(0, 5)) +
  scale_y_continuous(breaks = seq(0, 5, by = 1)) +
  guides(color = F) + ylab("Spikes")
g5

plot_list <- list(g1, g2, g3, g4, g5)

plots <- lapply(plot_list, function(p) {
  p + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"))
})

plot <- wrap_plots(plots, nrow = 1) + plot_annotation(title = "90 vs. 180")
ggsave("tuning_90_180.png", plot, width = 10, height = 2, dpi = 300)