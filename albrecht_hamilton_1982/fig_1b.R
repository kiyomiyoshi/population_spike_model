# Simulation based on Albrecht & Hamilton (1982, Fig. 1B)

library(tidyverse)
library(doParallel)
library(foreach)

Contrast <- c(2, 4, 8, 16, 24, 40, 56, 80, 96)
Rmax <- 115 # fixed at 115 for comparability
C50 <- 10.5
n   <-  1.7

Max_firing <- Rmax * Contrast^n / (Contrast^n + C50^n)
df_ah <- data.frame(Contrast, Max_firing)

g1 <- ggplot(df_ah, aes(x = Contrast, y = Max_firing)) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(0, 100)) +
  theme_minimal()
g1            

g2 <- ggplot(df_ah, aes(x = Contrast, y = Max_firing)) +
  geom_point(size = 2) +
  scale_x_log10() +
  theme_minimal()
g2 

### Simulation ###
cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

n_neurons <- 180                                               # Number of neurons in the population
orientations <- seq(1, 180, by = 1)                            # Possible orientations (0 to 359 degrees)
preferred_orientations <-  seq(1, 180, length.out = n_neurons) # Preferred orientation of each neuron
max_firing_seq <- df_ah$Max_firing                             # Maximum firing rate of each neuron
contrast <-       df_ah$Contrast
tuning_width    <- 20                                          # Tuning width (standard deviation) of each neuron's response curve
n_trials <- 10000

all_results <- list()

for (max_firing_rate in max_firing_seq) {
  
  cat("Running max firing rate =", max_firing_rate, "\n")
  
  tuning_curves <- matrix(0, nrow = n_neurons, ncol = length(orientations))
  
  for (i in 1:n_neurons) {
    tuning_curves[i, ] <- max_firing_rate * exp(
      -0.5 * (
        pmin(abs(orientations - preferred_orientations[i]),
             180 - abs(orientations - preferred_orientations[i]))^2
      ) / tuning_width^2
    )
  }
  
  # ---- 90 ----
  responses_90 <- foreach(
    i = 1:n_trials,
    .combine = rbind,
    .packages = "stats"
  ) %dopar% {
    
    input_90 <- rnorm(n_trials, 90, 0)
    idx <- round(input_90[i]) + 1
    idx <- pmin(pmax(idx, 1), 180)
    
    mu <- tuning_curves[, idx]
    mu[mu == 0] <- 1e-8
    
    resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
    
    cbind(resp,
          neuron = seq_len(n_neurons),
          stim = 90,
          trial = i,
          maxFR = max_firing_rate)
  }
  
  # ---- 91 ----
  responses_91 <- foreach(
    i = 1:n_trials,
    .combine = rbind,
    .packages = "stats"
  ) %dopar% {
    
    input_91 <- rnorm(n_trials, 91, 0)
    idx <- round(input_91[i]) + 1
    idx <- pmin(pmax(idx, 1), 180)
    
    mu <- tuning_curves[, idx]
    mu[mu == 0] <- 1e-8
    
    resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
    
    cbind(resp,
          neuron = seq_len(n_neurons),
          stim = 91,
          trial = i,
          maxFR = max_firing_rate)
  }
  
  df_tmp <- as.data.frame(rbind(responses_90, responses_91))
  colnames(df_tmp) <- c("Spikes", "Neuron", "Stimulus", "Trial", "Contrast")
  
  all_results[[as.character(max_firing_rate)]] <- df_tmp
}

df <- do.call(rbind, all_results)
df$Contrast <- contrast[match(df$Contrast, max_firing_seq)]

### Visualization ###
df %>%
  group_by(Stimulus, Contrast, Trial) %>%
  summarise(Sum_spikes = sum(Spikes), .groups = "drop") -> df_sum

g3 <- ggplot(df_sum, aes(x = Sum_spikes, color = factor(Contrast))) +
  geom_density(alpha = 0.2, linewidth = 1) +
  theme_classic() +
  labs(
    x = "Sum of spikes",
    y = "Density",
    color = "Contrast",
    fill =  "Contrast"
  )
g3

df_sum %>%
  group_by(Contrast) %>%
  summarise(Mean_sum_spikes = mean(Sum_spikes),
            SD_sum_spikes = sd(Sum_spikes), .groups = "drop") -> df_summary

df_summary %>%
  ggplot(aes(x = Contrast, y = Mean_sum_spikes)) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(0, 100)) +
  theme_minimal() -> g4
g4

df_summary %>%
  ggplot(aes(x = Contrast, y = Mean_sum_spikes)) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() -> g5
g5

df_summary %>%
  ggplot(aes(x = Contrast, y = SD_sum_spikes)) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(0, 100)) +
  theme_minimal() -> g6
g6

df_summary %>%
  ggplot(aes(x = Contrast, y = SD_sum_spikes)) +
  geom_point(size = 2) +
  scale_x_log10() +
  theme_minimal() -> g7
g7

df %>%
  filter(Neuron %in% c(90, 91)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `90`, y = `91`, color = factor(Stimulus))) +
  geom_point(size = 2, alpha = 0.2) +
  labs(x = "Neuron 90", y = "Neuron 91") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  theme_minimal() +
  facet_wrap(. ~ Contrast) -> g8
# g8

# sanity check for fano factor 
df %>%
  group_by(Contrast, Stimulus, Neuron) %>%
  summarise(
    mean_spikes = mean(Spikes),
    var_spikes  = var(Spikes),
    fano = var_spikes / mean_spikes,
    .groups = "drop"
  ) %>%
  filter(mean_spikes > 0) %>%
  group_by(Contrast, Stimulus) %>%
  summarise(
    mean_mean_spikes = mean(mean_spikes, na.rm = TRUE),
    mean_var_spikes = mean(var_spikes),
    mean_fano = mean(fano, na.rm = TRUE), # average across neurons
    sd_fano   = sd(fano, na.rm = TRUE),   # sd across neurons
    .groups = "drop")

### Linear Fisher information ###
df_wide <- df %>%
  pivot_wider(
    id_cols = c(Trial, Contrast, Stimulus),
    names_from = Neuron,
    values_from = Spikes)

df_wide$Stimulus_bin <- ifelse(df_wide$Stimulus == 90, 1, 0)

compute_lfi <- function(data, lambda = 1e-6) {
  
  neurons <- grep("^[0-9]+$", colnames(data), value = TRUE)  
  
  data0 <- data[data$Stimulus_bin == 0, neurons]
  data1 <- data[data$Stimulus_bin == 1, neurons]
  
  mu0 <- colMeans(data0)
  mu1 <- colMeans(data1)
  
  df <- mu1 - mu0
  
  Sigma0 <- cov(data0)
  Sigma1 <- cov(data1)
  Sigma  <- (Sigma0 + Sigma1) / 2
  
  p <- length(neurons)              # 逆行列を安定化させる処理 (ニューロン間の相関を少し弱める)
  Sigma <- Sigma + lambda * diag(p) # Ridge正則化, Tikhonov regularization (小さい固有値を持ち上げる)
  
  I <- as.numeric(t(df) %*% solve(Sigma) %*% df) # LFI
  I_per_neuron <- I / p
  
  return(I_per_neuron)
}

maxFR_values <- sort(unique(df_wide$Contrast))

lfi_results <- data.frame(
  Contrast = numeric(),
  LFI_per_neuron = numeric()
)

for (maxFR in maxFR_values) {
  df_sub <- df_wide[df_wide$Contrast == maxFR, ]
  lfi_value <- compute_lfi(df_sub)
  lfi_results <- rbind(lfi_results,
                       data.frame(Contrast = maxFR,
                                  LFI_per_neuron = lfi_value))
}

lfi_results$Contrast <- contrast
print(lfi_results)

lfi_results %>%
  ggplot(aes(x = Contrast, y = LFI_per_neuron)) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(0, 100)) +
  theme_minimal() -> g9
g9

### Decoding with logistic regression ###
run_logistic_cv_parallel <- function(data, k = 5, ncores = detectCores() - 1) {
  
  predictors <- grep("^[0-9]+$", colnames(data), value = TRUE)  
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

decoding_results <- data.frame(
  Contrast = numeric(),
  Decoding_accuracy = numeric()
)

for (maxFR in maxFR_values) {
  df_sub <- df_wide[df_wide$Contrast == maxFR, ]
  accuracy <- run_logistic_cv_parallel(df_sub, k = 5, ncores = detectCores() - 1)
  decoding_results <- rbind(decoding_results,
                            data.frame(Contrast = maxFR,
                                       Decoding_accuracy = accuracy))
}

decoding_results <- mutate(decoding_results, D_prime = qnorm(Decoding_accuracy) - qnorm(1 - Decoding_accuracy),
                           Contrast = contrast)
decoding_results

decoding_results %>%
  ggplot(aes(x = Contrast, y = D_prime)) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(0, 100)) +
  theme_minimal() -> g10
g10

decoding_results %>%
  ggplot(aes(x = Contrast, y = D_prime)) +
  geom_point(size = 2) +
  scale_x_log10() +
  theme_minimal() -> g11
g11

g <- cowplot::plot_grid(g1, g2, g10, g11, g4, g5, nrow = 3)
ggsave("fig_1b.png", g, width = 4, height = 5.6, dpi = 300)