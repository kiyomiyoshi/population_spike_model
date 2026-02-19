library(tidyverse)
library(doParallel)
library(foreach)

cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

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
tuning_width    <- 40                                          # Tuning width (standard deviation) of each neuron's response curve
n_trials <- 10000

# Spontaneous firing need to be considered
# spontaneous_firing_rate <- 0
# mu_stim  <- tuning_curves[, input_91[i] + 1]
# mu_total <- mu_stim + spontaneous_firing_rate
# size <- mu_total / 0.05   # Fano = 1.05
# resp <- rnbinom(n_neurons, size = size, mu = mu_total)

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
  ylim(0, 65) + guides(color = F) + ylab("Firing rate")

responses <- c()

input_90 <- rnorm(n_trials, 90, 0) # stimulus encoding noise
responses_90 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_90[i]) + 1
  mu <- tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 90, trial = i)
}

input_91 <- rnorm(n_trials, 91, 0)
responses_91 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_91[i]) + 1
  mu <- tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 91, trial = i)
}

df_responses_1 <- as.data.frame(rbind(responses_90, responses_91))
colnames(df_responses_1) <- c("Spikes", "Neuron", "Stimulus", "Trial")


### Unattended ###
# Parameters (unattended)
n_neurons <- 360                                               
orientations <- seq(1, 360, by = 1)                            
preferred_orientations <-  seq(1, 360, length.out = n_neurons) 
max_firing_rate <- 60 * 0.8   # 0.95                                   
tuning_width <-    40 * 1.25  # 1.05                                     

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
  ylim(0, 65) + guides(color = F) + ylab("Firing rate")

input_90 <- rnorm(n_trials, 90, 0) # stimulus encoding noise
responses_90 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_90[i]) + 1
  mu <- tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  size <- 5 * mu # Fano factor = 1.2
  resp <- rnbinom(n_neurons, size = size, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 90, trial = i)
}

input_91 <- rnorm(n_trials, 91, 0)
responses_91 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_91[i]) + 1
  mu <- tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  size <- 5 * mu     # Fano factor = 1.2
  resp <- rnbinom(n_neurons, size = size, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 91, trial = i)
}

df_responses_2 <- as.data.frame(rbind(responses_90, responses_91))
colnames(df_responses_2) <- c("Spikes", "Neuron", "Stimulus", "Trial")


###
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

# tuning curves
g1 <- cowplot::plot_grid(p1, p2)
g1

g2 <- ggplot(df_sum, aes(x = Sum_spikes, color = Attention, fill = Attention)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  theme_classic() +
  labs(
    x = "Sum of spikes",
    y = "Density",
    color = "Attention",
    fill = "Attention")
g2

df %>%
  filter(Neuron %in% c(90, 91)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `90`, y = `91`, color = factor(Stimulus))) +
  geom_point(size = 3, alpha = 0.2) +
  labs(x = "Neuron 90", y = "Neuron 91") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  theme_minimal() +
  facet_wrap(. ~ Attention) -> p3
p3

# sanity check for fano factor 
fano_summary <- df %>%
  group_by(Attention, Stimulus, Neuron) %>%
  summarise(
    mean_spikes = mean(Spikes),
    var_spikes  = var(Spikes),
    fano = var_spikes / mean_spikes,
    .groups = "drop"
  ) %>%
  filter(mean_spikes > 0) %>%
  group_by(Attention, Stimulus) %>%
  summarise(
    mean_fano = mean(fano, na.rm = TRUE), # average across neurons
    sd_fano   = sd(fano, na.rm = TRUE),   # sd across neurons
    .groups = "drop"
  )

fano_summary

# decoding with logistic regression
df_wide <- reshape(
  df,
  idvar = c("Trial", "Attention", "Stimulus"),
  timevar = "Neuron",
  direction = "wide"
)

colnames(df_wide) <- gsub("Spikes.", "N", colnames(df_wide), fixed = TRUE)
df_wide$Stimulus_bin <- ifelse(df_wide$Stimulus == 90, 1, 0)

run_logistic_cv_parallel <- function(data, k = 5, ncores = detectCores() - 1) {
  
  predictors <- grep("^N", colnames(data), value = TRUE)
  n <- nrow(data)
  
  set.seed(123)
  folds <- sample(rep(1:k, length.out = n))
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  acc <- foreach(i = 1:k, .combine = c) %dopar% {
    
    train <- data[folds != i, ]
    test  <- data[folds == i, ]
    
    model <- glm(
      Stimulus_bin ~ .,
      data = train[, c("Stimulus_bin", predictors)],
      family = binomial
    )
    
    prob <- predict(model, newdata = test, type = "response")
    pred <- ifelse(prob > 0.5, 1, 0)
    
    mean(pred == test$Stimulus_bin)
  }
  
  stopCluster(cl)
  
  mean(acc)
}

df_att   <- subset(df_wide, Attention == "attended")
df_unatt <- subset(df_wide, Attention == "unattended")

acc_att   <- run_logistic_cv_parallel(df_att, k = 5, ncores = detectCores() - 1)
acc_unatt <- run_logistic_cv_parallel(df_unatt, k = 5, ncores = detectCores() - 1)

acc_att
acc_unatt