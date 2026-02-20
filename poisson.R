library(tidyverse)

# For visualization
color_wheel <- function(alpha = 1) {
  h <- c(seq(320, 0, -1), seq(360, 321, -1))[-321]
  h <- h / 360
  s <- rep(1, 360)
  v <- rep(1, 360)
  hsv(h, s, v, alpha)
}

colors <- color_wheel()

### Attended ###
# Parameters
n_neurons <- 360                                               # Number of neurons in the population
orientations <- seq(1, 360, by = 1)                            # Possible orientations (0 to 359 degrees)
preferred_orientations <-  seq(1, 360, length.out = n_neurons) # Preferred orientation of each neuron
max_firing_rate <- 60                                          # Maximum firing rate of each neuron
tuning_width <- 40                                             # Tuning width (standard deviation) of each neuron's response curve
spontaneous_firing_rate <- 0                                   

# Tuning curves
tuning_curves <- matrix(0, nrow = n_neurons, ncol = length(orientations))
for (i in 1:n_neurons) {
  tuning_curves[i, ] <- max_firing_rate * exp(-0.5 * (pmin(abs(orientations - preferred_orientations[i]), 
                                                           360 - abs(orientations - preferred_orientations[i])) ^ 2) / tuning_width ^ 2)
}

FiringRate <- c()
for (i in 1:n_neurons) {
  FiringRate <- c(FiringRate, tuning_curves[i, ])
}

df_tuning <- data.frame(Orientation = rep(orientations, n_neurons),
                        FiringRate = FiringRate,
                        Neuron = rep(rep(orientations, n_neurons), each = length(orientations)),
                        Color = rep(colors, each = length(orientations)))

p1 <- ggplot(df_tuning[df_tuning$Neuron %in% seq(1, n_neurons, by = 60), ], 
             aes(x = Orientation, y = FiringRate, color = Color)) +
  geom_line() + scale_color_identity() +
  theme_minimal() + 
  ylim(0, 60) + guides(color = F) + ylab("Firing rate")

# Manifold
responses <- c()

input_90 <- rnorm(1000, 90, 0) # stimulus input plus internal noise
for (i in 1:1000) {
  resp <- rpois(n_neurons, lambda = tuning_curves[, input_90[i] + 1]) + spontaneous_firing_rate
  responses <- rbind(responses, cbind(resp, seq(1, 360, 1), 90, i))
}

input_120 <- rnorm(1000, 120, 0) # stimulus input plus internal noise
for (i in 1:1000) {
  resp <- rpois(n_neurons, lambda = tuning_curves[, input_120[i] + 1]) + spontaneous_firing_rate
  responses <- rbind(responses, cbind(resp, seq(1, 360, 1), 120, i))
}

df_responses_1 <- as.data.frame(responses)
colnames(df_responses_1) <- c("Spikes", "Neuron", "Stimulus", "Trial")

### Unattended ###
# Parameters (unattended)
n_neurons <- 360                                               
orientations <- seq(1, 360, by = 1)                           
preferred_orientations <-  seq(1, 360, length.out = n_neurons)
max_firing_rate <- 60 * 0.8                                          
tuning_width <- 40 * 1.25                                            
spontaneous_firing_rate <- 0                                  

# Tuning curves
tuning_curves <- matrix(0, nrow = n_neurons, ncol = length(orientations))
for (i in 1:n_neurons) {
  tuning_curves[i, ] <- max_firing_rate * exp(-0.5 * (pmin(abs(orientations - preferred_orientations[i]), 
                                                           360 - abs(orientations - preferred_orientations[i])) ^ 2) / tuning_width ^ 2)
}

FiringRate <- c()
for (i in 1:n_neurons) {
  FiringRate <- c(FiringRate, tuning_curves[i, ])
}

df_tuning <- data.frame(Orientation = rep(orientations, n_neurons),
                        FiringRate = FiringRate,
                        Neuron = rep(rep(orientations, n_neurons), each = length(orientations)),
                        Color = rep(colors, each = length(orientations)))

p2 <- ggplot(df_tuning[df_tuning$Neuron %in% seq(1, n_neurons, by = 60), ], 
             aes(x = Orientation, y = FiringRate, color = Color)) +
  geom_line() + scale_color_identity() +
  theme_minimal() + 
  ylim(0, 60) + guides(color = F) + ylab("Firing rate")

# Manifold
responses <- c()

input_90 <- rnorm(1000, 90, 0) # stimulus input plus internal noise
for (i in 1:1000) {
  resp <- rpois(n_neurons, lambda = tuning_curves[, input_90[i] + 1]) + spontaneous_firing_rate # important
  responses <- rbind(responses, cbind(resp, seq(1, 360, 1), 90, i))
}

input_120 <- rnorm(1000, 120, 0) # stimulus input plus internal noise
for (i in 1:1000) {
  resp <- rpois(n_neurons, lambda = tuning_curves[, input_120[i] + 1]) + spontaneous_firing_rate # important
  responses <- rbind(responses, cbind(resp, seq(1, 360, 1), 120, i))
}

df_responses_2 <- as.data.frame(responses)
colnames(df_responses_2) <- c("Spikes", "Neuron", "Stimulus", "Trial")

### Visualization ###
p <- cowplot::plot_grid(p1, p2)
p

df_responses_1$Attention <- "attended"
df_responses_2$Attention <- "unattended"
df <- rbind(df_responses_1, df_responses_2)

df %>%
  group_by(Stimulus, Attention, Trial) %>%
  summarise(Sum_spikes = sum(Spikes), .groups = "drop") -> df_sum

df_sum %>%
  group_by(Stimulus, Attention) %>%
  summarise(Mean_sum_spikes = mean(Sum_spikes),
            SD_sum_spikes = sd(Sum_spikes), .groups = "drop") 