library(tidyverse)
library(doParallel)
library(foreach)
library(cowplot)
library(patchwork)

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
  scale_x_continuous(
    limits = c(1, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  coord_cartesian(ylim = c(0, 130)) +
  scale_y_continuous(breaks = seq(0, 120, by = 60)) +
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

# High contrast
n_neurons <- 180                                               # Number of neurons in the population
orientations <- seq(1, 180, by = 1)                            # Possible orientations (1 to 180 degrees)
preferred_orientations <-  seq(1, 180, length.out = n_neurons) # Preferred orientation of each neuron
max_firing_rate <- 115                                         # Maximum firing rate of each neuron
tuning_width    <- 20                                          # Tuning width (standard deviation) of each neuron's response curve
n_trials <- 10000

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
  theme_classic(base_size = 12) + 
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  ylim(0, 120) + 
  guides(color = F) + ylab("Spikes")
g2

input_90 <- rnorm(n_trials, 90, 0) # zero stimulus encoding noise
responses_90 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_90[i]) 
  mu <- tuning_curves[, idx]
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
  mu <- tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 100, trial = i)
}

df_responses_1 <- as.data.frame(rbind(responses_90, responses_100))
colnames(df_responses_1) <- c("Spikes", "Neuron", "Stimulus", "Trial")

# Low contrast
n_neurons <- 180                                               
orientations <- seq(1, 180, by = 1)                            
preferred_orientations <-  seq(1, 180, length.out = n_neurons) 
max_firing_rate <- Rmax * 20^n / (20^n + C50^n) # 20% contrast                                 
tuning_width <-    20                                     

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
  theme_classic(base_size = 12) + 
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  ylim(0, 120) + 
  guides(color = F) + ylab("Spikes")
g3

input_90 <- rnorm(n_trials, 90, 0) # zero stimulus encoding noise
responses_90 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_90[i]) 
  mu <- tuning_curves[, idx]
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
  mu <- tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 100, trial = i)
}

df_responses_2 <- as.data.frame(rbind(responses_90, responses_100))
colnames(df_responses_2) <- c("Spikes", "Neuron", "Stimulus", "Trial")

# stopCluster(cl)

### Visualization ###
df_responses_1$Contrast <- "100"
df_responses_2$Contrast <- "20"
df <- rbind(df_responses_1, df_responses_2)

df %>%
  group_by(Stimulus, Contrast, Trial) %>%
  summarise(Sum_spikes = sum(Spikes), .groups = "drop") -> df_sum

df_sum %>%
  group_by(Stimulus, Contrast) %>%
  summarise(Mean_sum_spikes = mean(Sum_spikes),
            SD_sum_spikes = sd(Sum_spikes), .groups = "drop")

# tuning curves
g4 <- ggplot(df_sum, aes(x = Sum_spikes, color = Contrast, fill = Contrast)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  theme_classic() +
  labs(
    x = "Sum of spikes",
    y = "Density",
    color = "Contrast",
    fill = "Contrast")
g4

g5 <- df %>%
  filter(Neuron %in% c(90, 100)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `90`, y = `100`, color = factor(Stimulus))) +
  geom_point(size = 3, alpha = 0.2) +
  labs(x = "Neuron 90", y = "Neuron 100") +
  coord_cartesian(xlim = c(0, 150), ylim = c(0, 150)) +
  theme_minimal() +
  facet_wrap(. ~ Contrast)
g5