# Contrastについてなだらかな彩度変化に変更
# 3Dプロットは少数の点の散布図に変更し、90 vs. 100を異なる色相にする

library(tidyverse)
library(doParallel)
library(foreach)
library(cowplot)
library(patchwork)
library(plotly)

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
n_trials <- 50000

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
max_firing_rate <- Rmax * 15^n / (15^n + C50^n) # 15% contrast                                 
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
  ggplot(aes(x = `90`, y = `100`,
             color = Contrast)) +
  stat_ellipse(aes(group = interaction(Contrast, Stimulus)), level = c(0.5)) + # parametric density
  stat_ellipse(aes(group = interaction(Contrast, Stimulus)), level = c(0.7)) + 
  stat_ellipse(aes(group = interaction(Contrast, Stimulus)), level = c(0.9)) + 
  labs(x = "Neuron 90", y = "Neuron 100") +
  coord_cartesian(xlim = c(0,150), ylim = c(0,150)) +
  scale_color_manual(values = c("#2C2C7A", "#E69F00")) +
  theme_classic() +
  theme(
    legend.position = c(0.4, 1),
    legend.justification = c(1, 1)
  )

g6 <- df %>%
  filter(Neuron %in% c(90, 100)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `90`, y = `100`, color = Contrast)) +
  geom_density_2d(aes(linetype = factor(Stimulus)), linewidth = 1) + # nonparametric density
  labs(x = "Neuron 90", y = "Neuron 100",
       color = "Contrast", linetype = "Stimulus") +
  coord_cartesian(xlim = c(0, 150), ylim = c(0, 150)) +
  scale_color_manual(values = c("#2C2C7A", "#E69F00")) +
  theme_classic() +
  theme(
    legend.position = c(0.4, 1),
    legend.justification = c(1, 1)
  ) +
  theme(
    legend.position = c(0.4, 1),
    legend.justification = c(1, 1)
  )

df3 <- df %>% 
  filter(Neuron %in% c(80, 90, 100)) %>% 
  pivot_wider(id_cols = c(Contrast, Stimulus, Trial), 
              names_from = Neuron, values_from = Spikes)

u <- seq(0, 2*pi, length.out = 50)
v <- seq(0, pi, length.out = 50)
x_sphere <- outer(cos(u), sin(v))
y_sphere <- outer(sin(u), sin(v))
z_sphere <- outer(rep(1,length(u)), cos(v)) 
unit_sphere <- rbind(as.vector(x_sphere), as.vector(y_sphere), as.vector(z_sphere))

ellipsoids <- df3 %>% 
  group_by(Stimulus, Contrast) %>% 
  group_map(~{
    
    mu <- colMeans(.x[,c("80","90","100")])
    
    sigma <- cov(.x[,c("80","90","100")])
    if(det(sigma) <= 0){
      sigma <- sigma + diag(1e-6,3)
    }
    
    eig <- eigen(sigma)
    
    coords <- eig$vectors %*% diag(sqrt(eig$values)) %*% unit_sphere + mu
    
    list(
      x = matrix(coords[1,], nrow=length(u)),
      y = matrix(coords[2,], nrow=length(u)),
      z = matrix(coords[3,], nrow=length(u)),
      Stimulus = .y$Stimulus,
      Contrast = .y$Contrast
    )
    
  })

p <- plot_ly()

for (e in ellipsoids) {
  
  col <- ifelse(e$Contrast == 20, "Reds", "Blues")  # Contrastで色分け
  
  p <- p %>% add_surface(
    x = e$x,
    y = e$y,
    z = e$z,
    colorscale = col,
    showscale = FALSE,
    opacity = 0.6,
    name = paste("Stim:", e$Stimulus, "Con:", e$Contrast),
    hoverinfo = "text",
    text = paste("Stimulus:", e$Stimulus,
                 "<br>Contrast:", e$Contrast)
  )
}

p <- p %>%
  layout(
    scene = list(
      xaxis = list(title = "Neuron 80", range = c(0,150)),
      yaxis = list(title = "Neuron 90", range = c(0,150)),
      zaxis = list(title = "Neuron 100", range = c(0,150))
    )
  )

p