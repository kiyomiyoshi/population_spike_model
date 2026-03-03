# Simulation assuming Naka-Rushton function on Albrecht & Hamilton (1982, Fig. 4)

library(tidyverse)
library(doParallel)
library(foreach)
library(cowplot)


#################### High Attention ####################
Contrast <- c(0.5, 1, 4, seq(10, 100, 15))
Rmax <- 115
C50 <- 19.3
n   <-  2.9
R0 <- Rmax * 0.03 # spontaneous firing based on Geisler & Albrecht (1997)

Max_firing <- Rmax * Contrast^n / (Contrast^n + C50^n)
df_ah <- data.frame(Contrast, Max_firing)

h1 <- ggplot(df_ah, aes(x = Contrast, y = Max_firing)) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 130)) +
  ylab("Tuning curve peak") +
  theme_minimal()
h1            

h2 <- ggplot(df_ah, aes(x = Contrast, y = Max_firing)) +
  geom_point(size = 2) +
  scale_x_log10() +
  ylab("Tuning curve peak") +
  theme_minimal()
h2 

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
    
    resp <- rnbinom(n_neurons, size = 1000000, mu = mu) + rnbinom(n_neurons, size = 1000000, mu = R0)
    
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
    
    resp <- rnbinom(n_neurons, size = 1000000, mu = mu) + rnbinom(n_neurons, size = 1000000, mu = R0)
    
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
write.csv(df, 
        file = "high_attention.csv", 
        row.names = FALSE)


### Sanity check for fano factor ### 
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


### Contrast discrimination ###
df %>%
  group_by(Stimulus, Contrast, Trial) %>%
  summarise(Sum_spikes = sum(Spikes), .groups = "drop") -> df_sum

h3 <- ggplot(df_sum, aes(x = Sum_spikes, color = factor(Contrast))) +
  geom_density(alpha = 0.2, linewidth = 1) +
  theme_classic() +
  labs(
    x = "Sum of spikes",
    y = "Density",
    color = "Contrast",
    fill =  "Contrast"
  ) + 
  theme(
    legend.position = c(0.98, 1.05),
    legend.justification = c(1, 1),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8)
  ) +
  coord_cartesian(xlim = c(0, 10000), ylim = c(0, 0.03)) +
  ggtitle("High attention (Fano factor = 1.0)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
h3

df_sum %>%
  group_by(Contrast) %>%
  summarise(Mean_sum_spikes = mean(Sum_spikes),
            SD_sum_spikes = sd(Sum_spikes), .groups = "drop") -> df_summary
colnames(df_summary) <- c("Contrast", "mu", "sigma")

dC <-  diff(df_summary$Contrast)
dmu <- diff(df_summary$mu)

# μ'
mu_prime <- c(
  dmu[1] / dC[1],                                 　　　　　# 前方差分（最初）
  (dmu[-1]/dC[-1] + dmu[-length(dmu)]/dC[-length(dC)]) / 2, # 中央差分
  dmu[length(dmu)] / dC[length(dC)]                         # 後方差分（最後）
)
df_summary$mu_prime <- mu_prime

# μ'(C)ΔC ≈ σ
# μ'(C): コントラストの微小変化に対するμの変化
# ΔC:    コントラストの変化量
# μの変化がσと同じ大きさになったとき丁度弁別できる(d'=1に相当するJND)
df_summary$DeltaC_pred <- df_summary$sigma / df_summary$mu_prime # JND
df_summary <- mutate(df_summary, Weber_ratio_pred = DeltaC_pred / Contrast) # ΔC / C (Weber ratio)
df_summary

df_plot <- df_summary %>%
  slice(2:(n() - 1))

h4 <- ggplot(df_plot, aes(x = Contrast, y = DeltaC_pred)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  scale_x_log10() +
  scale_y_log10(limits = c(0.2, 35)) +
  ylab("ΔC (JND of d' = 1)") +
  xlab("Contrast")

h5 <- ggplot(df_plot, aes(x = Contrast, y = Weber_ratio_pred)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  coord_cartesian(ylim = c(0, 3)) +
  ylab("ΔC/C (Weber fraction)") +
  xlab("Contrast")

h6 <- df_summary %>%
  ggplot(aes(x = Contrast, y = mu)) +
  geom_point(size = 2) +
  # scale_x_log10() +
  # scale_y_log10() +
  coord_cartesian(ylim = c(0, 10000)) +
  theme_minimal() +
  ylab("Mean total spikes")

### Orientation discrimination (LFI) ###
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
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 0.07)) +
  ylab("LFI per neuron") +
  theme_minimal() -> h7
h7


### Orientation discrimination (decoding with logistic regression) ###
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
  ylab("d' for orientation") +
  coord_cartesian(xlim = c(0, 3.5)) +
  theme_minimal() -> h8
h8

decoding_results %>%
  ggplot(aes(x = Contrast, y = D_prime)) +
  geom_point(size = 2) +
  scale_x_log10() +
  ylab("d' for orientation") +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_minimal() -> h9
h9


### Subjective inflation ###
df_sum %>%
  dplyr::filter(Contrast == 1 | Contrast == 4) %>%
  summarise(Criterion = mean(Sum_spikes)) %>%
  as.numeric() -> criterion

h10 <- ggplot(subset(df_sum, df_sum$Contrast == 1 | df_sum$Contrast == 4), 
            aes(x = Sum_spikes, color = factor(Contrast))) +
  geom_vline(xintercept = criterion, linetype = "dashed") +
  geom_density(alpha = 0.2, linewidth = 1) +
  theme_classic() +
  labs(
    x = "Sum of spikes",
    y = "Density",
    color = "Contrast",
    fill =  "Contrast"
  ) + 
  theme(
    legend.position = c(0.98, 1.05),
    legend.justification = c(1, 1),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8)
  ) +
  coord_cartesian(xlim = c(500, 800), ylim = c(0, 0.02))
h10
  
h11 <- df_sum %>%
  dplyr::filter(Contrast == 1 | Contrast == 4) %>%
  mutate(Yes = Sum_spikes > criterion) %>%
  group_by(Contrast) %>%
  summarise(P_yes = mean(Yes)) %>%
  ggplot(aes(x = Contrast, y = P_yes, color = factor(Contrast))) +
  geom_point(size = 3) +
  ylim(0, 1) +
  theme_classic() +
  labs(
    x = "Contrast",
    y = "P(Yes)",
    color = "Contrast"
  ) + 
  theme(legend.position = "None")
h11



#################### Low Attention ####################
Contrast <- c(0.5, 1, 4, seq(10, 100, 15))
Rmax <- 115
C50 <- 19.3
n   <-  2.9
R0 <- Rmax * 0.03

Max_firing <- Rmax * Contrast^n / (Contrast^n + C50^n)
df_ah <- data.frame(Contrast, Max_firing)

l1 <- ggplot(df_ah, aes(x = Contrast, y = Max_firing)) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 130)) +
  ylab("Tuning curve peak") +
  theme_minimal()
l1            

l2 <- ggplot(df_ah, aes(x = Contrast, y = Max_firing)) +
  geom_point(size = 2) +
  scale_x_log10() +
  ylab("Tuning curve peak") +
  theme_minimal()
l2 


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
   
    resp <- rnbinom(n_neurons, size = 5 * mu, mu = mu) + rnbinom(n_neurons, size = 5 * R0, mu = R0) # Fano factor = 1.2 (var = μ + μ^2/size)1000000
    
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
    
    resp <- rnbinom(n_neurons, size = 5 * mu, mu = mu) + rnbinom(n_neurons, size = 5 * R0, mu = R0)
    
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
# write.csv(df, 
#           file = "low_attention.csv", 
#           row.names = FALSE)


### Sanity check for fano factor ### 
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


### Contrast discrimination ###
df %>%
  group_by(Stimulus, Contrast, Trial) %>%
  summarise(Sum_spikes = sum(Spikes), .groups = "drop") -> df_sum

l3 <- ggplot(df_sum, aes(x = Sum_spikes, color = factor(Contrast))) +
  geom_density(alpha = 0.2, linewidth = 1) +
  theme_classic() +
  labs(
    x = "Sum of spikes",
    y = "Density",
    color = "Contrast",
    fill =  "Contrast"
  ) + 
  theme(
    legend.position = c(0.98, 1.05),
    legend.justification = c(1, 1),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8)
  ) +
  coord_cartesian(xlim = c(0, 10000), ylim = c(0, 0.03)) +
  ggtitle("Low attention (Fano factor = 1.2)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
l3

df_sum %>%
  group_by(Contrast) %>%
  summarise(Mean_sum_spikes = mean(Sum_spikes),
            SD_sum_spikes = sd(Sum_spikes), .groups = "drop") -> df_summary
colnames(df_summary) <- c("Contrast", "mu", "sigma")

dC <-  diff(df_summary$Contrast)
dmu <- diff(df_summary$mu)

# μ'
mu_prime <- c(
  dmu[1] / dC[1],                                 　　　　　# 前方差分（最初）
  (dmu[-1]/dC[-1] + dmu[-length(dmu)]/dC[-length(dC)]) / 2, # 中央差分
  dmu[length(dmu)] / dC[length(dC)]                         # 後方差分（最後）
)
df_summary$mu_prime <- mu_prime

# μ'(C)ΔC ≈ σ
# μ'(C): コントラストの微小変化に対するμの変化
# ΔC:    コントラストの変化量
# μの変化がσと同じ大きさになったとき丁度弁別できる(d'=1に相当するJND)
df_summary$DeltaC_pred <- df_summary$sigma / df_summary$mu_prime # JND
df_summary <- mutate(df_summary, Weber_ratio_pred = DeltaC_pred / Contrast) # ΔC / C (Weber ratio)
df_summary

df_plot <- df_summary %>%
  slice(2:(n() - 1))

l4 <- ggplot(df_plot, aes(x = Contrast, y = DeltaC_pred)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  scale_x_log10() +
  scale_y_log10(limits = c(0.2, 35)) +
  ylab("ΔC (JND of d' = 1)") +
  xlab("Contrast")

l5 <- ggplot(df_plot, aes(x = Contrast, y = Weber_ratio_pred)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  coord_cartesian(ylim = c(0, 3)) +
  ylab("ΔC/C (Weber fraction)") +
  xlab("Contrast")

l6 <- df_summary %>%
  ggplot(aes(x = Contrast, y = mu)) +
  geom_point(size = 2) +
  # scale_x_log10() +
  # scale_y_log10() +
  coord_cartesian(ylim = c(0, 10000)) +
  theme_minimal() +
  ylab("Mean total spikes")


### Orientation discrimination (LFI) ###
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
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 0.07)) +
  ylab("LFI per neuron") +
  theme_minimal() -> l7
l7


### Orientation discrimination (decoding with logistic regression) ###
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
  ylab("d' for orientation") +
  coord_cartesian(xlim = c(0, 3.5)) +
  theme_minimal() -> l8
l8

decoding_results %>%
  ggplot(aes(x = Contrast, y = D_prime)) +
  geom_point(size = 2) +
  scale_x_log10() +
  ylab("d' for orientation") +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_minimal() -> l9
l9

### Subjective inflation ###
l10 <- ggplot(subset(df_sum, df_sum$Contrast == 1 | df_sum$Contrast == 4), 
              aes(x = Sum_spikes, color = factor(Contrast))) +
  geom_vline(xintercept = criterion, linetype = "dashed") +
  geom_density(alpha = 0.2, linewidth = 1) +
  theme_classic() +
  labs(
    x = "Sum of spikes",
    y = "Density",
    color = "Contrast",
    fill =  "Contrast"
  ) + 
  theme(
    legend.position = c(0.98, 1.05),
    legend.justification = c(1, 1),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8)
  ) +
  coord_cartesian(xlim = c(500, 800), ylim = c(0, 0.02))
l10

l11 <- df_sum %>%
  dplyr::filter(Contrast == 1 & Contrast == 4) %>%
  mutate(Yes = Sum_spikes > criterion) %>%
  group_by(Contrast) %>%
  summarise(P_yes = mean(Yes)) %>%
  ggplot(aes(x = Contrast, y = P_yes, color = factor(Contrast))) +
  geom_point(size = 3) +
  ylim(0, 1) +
  theme_classic() +
  labs(
    x = "Contrast",
    y = "P(Yes)",
    color = "Contrast"
  ) + 
  theme(legend.position = "None")
l11

### Figures ###
p1 <- plot_grid(h3, ncol = 1)
p2 <- plot_grid(h1, h6, ncol = 2)
p3 <- plot_grid(h4, h5, ncol = 2)
p4 <- plot_grid(h7, h9, ncol = 2)
p5 <- plot_grid(h10, h11, ncol = 2)
high_alpha <- plot_grid(p1, p2, p3, p4, p5, ncol = 1, rel_heights = c(1, 1))

p1 <- plot_grid(l3, ncol = 1)
p2 <- plot_grid(l1, l6, ncol = 2)
p3 <- plot_grid(l4, l5, ncol = 2)
p4 <- plot_grid(l7, l9, ncol = 2)
p5 <- plot_grid(l10, l11, ncol = 2)
low_alpha <- plot_grid(p1, p2, p3, p4, p5, ncol = 1, rel_heights = c(1, 1))
subjective_inflation <- plot_grid(low_alpha, high_alpha, ncol = 2)
ggsave("subjective_inflation.png", subjective_inflation, width = 12, height = 11.25, dpi = 300)