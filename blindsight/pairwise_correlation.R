library(tidyverse)
library(doParallel)
library(foreach)
library(cowplot)
library(patchwork)
library(plotly)
library(magick)
library(webshot2)
library(minpack.lm)

cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

### Naka-Rushton function ###
Contrast <- seq(0, 100, 0.1)
Rmax <- 115
C50 <- 19.3
n   <-  2.9
Max_firing <- Rmax * Contrast^n / (Contrast^n + C50^n)
nr <- data.frame(Contrast, Max_firing)

g1 <- ggplot(nr, aes(x = Contrast, y = Max_firing)) +
  geom_line(alpha = 0.85, linewidth = 0.7) +
  annotate(
    "text", 
    x = 50,              
    y = 135,             
    label = "Naka-Rushton", 
    vjust = 1,           
    hjust = 0.5,         
    size = 3.5             
  ) +
  scale_x_continuous(
    limits = c(1, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  coord_cartesian(ylim = c(0, 130)) +
  scale_y_continuous(breaks = seq(0, 120, by = 30)) +
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

# High gain variability
n_neurons <- 180                                              
orientations <- seq(1, 180, by = 1)                           
preferred_orientations <-  seq(1, 180, length.out = n_neurons)
max_firing_rate <- Rmax * 15^n / (15^n + C50^n) # 15% contrast
tuning_width    <- 20                                         
n_trials <- 30000

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
  annotate(
    "text", 
    x = 90,              
    y = 135,             
    label = "Contrast = 15%", 
    vjust = 1,           
    hjust = 0.5,         
    size = 3.5             
  ) +
  theme_classic(base_size = 11) + 
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  coord_cartesian(ylim = c(0, 130)) +
  scale_y_continuous(breaks = seq(0, 120, by = 30)) +
  guides(color = F) + ylab("Spikes")
g2

input_90 <- rnorm(n_trials, 90, 0) # zero stimulus encoding noise
responses_90 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_90[i]) 
  sigma_g <- 0.2
  g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
  mu <- g * tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 90, trial = i)
}

input_100 <- rnorm(n_trials, 100, 0)
responses_100 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_100[i]) 
  sigma_g <- 0.2
  g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
  mu <- g * tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 100, trial = i)
}

df_responses_1 <- as.data.frame(rbind(responses_90, responses_100))
colnames(df_responses_1) <- c("Spikes", "Neuron", "Stimulus", "Trial")

# Low gain variability
n_neurons <- 180                                               
orientations <- seq(1, 180, by = 1)                            
preferred_orientations <-  seq(1, 180, length.out = n_neurons) 
max_firing_rate <- Rmax * 15^n / (15^n + C50^n)                                 
tuning_width <-    20                                  

tuning_curves <- matrix(0, nrow = n_neurons, ncol = length(orientations))
for (i in 1:n_neurons) {
  tuning_curves[i, ] <- max_firing_rate * exp(-0.5 * (pmin(abs(orientations - preferred_orientations[i]), 
                                                           180 - abs(orientations - preferred_orientations[i])) ^ 2) / tuning_width ^ 2)
}

input_90 <- rnorm(n_trials, 90, 0) # zero stimulus encoding noise
responses_90 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_90[i]) 
  sigma_g <- 0.01
  g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
  mu <- g * tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 90, trial = i)
}

input_100 <- rnorm(n_trials, 100, 0)
responses_100 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_100[i]) 
  sigma_g <- 0.01
  g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
  mu <- g * tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 100, trial = i)
}

df_responses_2 <- as.data.frame(rbind(responses_90, responses_100))
colnames(df_responses_2) <- c("Spikes", "Neuron", "Stimulus", "Trial")

# stopCluster(cl)

### Visualization ###
# across-trial spike distributions
# 2d density
df_responses_1$GV <- "high"
df_responses_2$GV <- "low"
df <- rbind(df_responses_1, df_responses_2)

g3 <- df %>%
  dplyr::filter(Neuron %in% c(90, 91)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `90`, y = `91`,
             color = GV)) +
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.5)) + # parametric density
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.7)) + 
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.9)) + 
  labs(x = "Neuron 90", y = "Neuron 91") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
g3

df %>%
  dplyr::filter(Neuron %in% c(90, 91)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  group_by(GV, Stimulus) %>%
  summarise(
    r = cor(`90`, `91`, use = "complete.obs"),
    .groups = "drop") 


g4 <- df %>%
  dplyr::filter(Neuron %in% c(90, 100)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `90`, y = `100`,
             color = GV)) +
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.5)) + # parametric density
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.7)) + 
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.9)) + 
  labs(x = "Neuron 90", y = "Neuron 100") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
g4

df %>%
  dplyr::filter(Neuron %in% c(90, 100)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  group_by(GV, Stimulus) %>%
  summarise(
    r = cor(`90`, `100`, use = "complete.obs"),
    .groups = "drop") 


g5 <- df %>%
  dplyr::filter(Neuron %in% c(90, 120)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `90`, y = `120`,
             color = GV)) +
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.5)) + # parametric density
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.7)) + 
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.9)) + 
  labs(x = "Neuron 90", y = "Neuron 120") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
g5

df %>%
  dplyr::filter(Neuron %in% c(90, 120)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  group_by(GV, Stimulus) %>%
  summarise(
    r = cor(`90`, `120`, use = "complete.obs"),
    .groups = "drop") 


g6 <- df %>%
  dplyr::filter(Neuron %in% c(90, 140)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `90`, y = `140`,
             color = GV)) +
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.5)) + # parametric density
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.7)) + 
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.9)) + 
  labs(x = "Neuron 90", y = "Neuron 140") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
g6

df %>%
  dplyr::filter(Neuron %in% c(90, 140)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  group_by(GV, Stimulus) %>%
  summarise(
    r = cor(`90`, `140`, use = "complete.obs"),
    .groups = "drop") 


g7 <- df %>%
  dplyr::filter(Neuron %in% c(140, 145)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `140`, y = `145`,
             color = GV)) +
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.5)) + # parametric density
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.7)) + 
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.9)) + 
  labs(x = "Neuron 140", y = "Neuron 145") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
g7

df %>%
  dplyr::filter(Neuron %in% c(140, 145)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  group_by(GV, Stimulus) %>%
  summarise(
    r = cor(`140`, `145`, use = "complete.obs"),
    .groups = "drop") 


# save figures
plot_list <- list(g3, g4, g5, g6, g7)
plots <- lapply(plot_list, function(p) {
  p + theme(
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )
})
pc <- wrap_plots(plots, nrow = 1)
ggsave("pairwise_correlation.png", pc, width = 11.25, height = 2.5, dpi = 300)