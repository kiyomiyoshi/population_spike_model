# Simulation assuming Naka-Rushton function on Albrecht & Hamilton (1982, Fig. 4)

library(tidyverse)
library(doParallel)
library(foreach)
library(cowplot)

#################### High Attention ####################
Contrast <- c(1, 4, seq(10, 100, 15))
Rmax <- 115
C50 <- 19.3
n   <-  2.9
R0 <- Rmax * 0.03 # spontaneous firing based on Geisler & Albrecht (1997)

Max_firing <- Rmax * Contrast^n / (Contrast^n + C50^n)
df <- data.frame(Contrast, Max_firing)

### Simulation ###
cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

n_neurons <- 180                                              # Number of neurons in the population
orientations <- seq(1, 180, by = 1)                           # Possible orientations (0 to 359 degrees)
preferred_orientations <- seq(1, 180, length.out = n_neurons) # Preferred orientation of each neuron
max_firing_seq <- df$Max_firing                               # Maximum firing rate of each neuron
contrast <- df$Contrast
tuning_width <- 20                                            # Tuning width (standard deviation) of each neuron's response curve
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

df_high <- do.call(rbind, all_results)
df_high$Contrast <- contrast[match(df_high$Contrast, max_firing_seq)]

### Contrast discrimination ###
df_high %>%
  group_by(Stimulus, Contrast, Trial) %>%
  summarise(Sum_spikes = sum(Spikes), .groups = "drop") -> df_sum_high

df_sum_high %>%
  group_by(Contrast) %>%
  summarise(Mean_sum_spikes = mean(Sum_spikes),
            SD_sum_spikes = sd(Sum_spikes), .groups = "drop") -> df_summary_high
colnames(df_summary_high) <- c("Contrast", "mu", "sigma")

dC <-  diff(df_summary_high$Contrast)
dmu <- diff(df_summary_high$mu)

# μ'
mu_prime <- c(
  dmu[1] / dC[1],                                 　　　　　# 前方差分（最初）
  (dmu[-1]/dC[-1] + dmu[-length(dmu)]/dC[-length(dC)]) / 2, # 中央差分
  dmu[length(dmu)] / dC[length(dC)]                         # 後方差分（最後）
)
df_summary_high$mu_prime <- mu_prime

# μ'(C)ΔC ≈ σ
# μ'(C): コントラストの微小変化に対するμの変化
# ΔC:    コントラストの変化量
# μの変化がσと同じ大きさになったとき丁度弁別できる(d'=1に相当するJND)
df_summary_high$DeltaC_pred <- df_summary_high$sigma / df_summary_high$mu_prime # JND
df_summary_high <- mutate(df_summary_high, Weber_ratio_pred = DeltaC_pred / Contrast) # ΔC / C (Weber ratio)

### Orientation discrimination (LFI) ###
df_high_wide <- df_high %>%
  pivot_wider(
    id_cols = c(Trial, Contrast, Stimulus),
    names_from = Neuron,
    values_from = Spikes)

df_high_wide$Stimulus_bin <- ifelse(df_high_wide$Stimulus == 90, 1, 0)

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

maxFR_values <- sort(unique(df_high_wide$Contrast))

lfi_results_high <- data.frame(
  Contrast = numeric(),
  LFI_per_neuron = numeric()
)

for (maxFR in maxFR_values) {
  df_sub <- df_high_wide[df_high_wide$Contrast == maxFR, ]
  lfi_value <- compute_lfi(df_sub)
  lfi_results_high <- rbind(lfi_results_high,
                            data.frame(Contrast = maxFR,
                                       LFI_per_neuron = lfi_value))
}

lfi_results_high$Contrast <- contrast

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

decoding_results_high <- data.frame(
  Contrast = numeric(),
  Decoding_accuracy = numeric()
)

for (maxFR in maxFR_values) {
  df_sub <- df_high_wide[df_high_wide$Contrast == maxFR, ]
  accuracy <- run_logistic_cv_parallel(df_sub, k = 5, ncores = detectCores() - 1)
  decoding_results_high <- rbind(decoding_results_high,
                                 data.frame(Contrast = maxFR,
                                            Decoding_accuracy = accuracy))
}

decoding_results_high <- mutate(decoding_results_high, D_prime = qnorm(Decoding_accuracy) - qnorm(1 - Decoding_accuracy),
                                Contrast = contrast)


#################### Low Attention ####################
Contrast <- c(1, 4, seq(10, 100, 15))
Rmax <- 115
C50 <- 19.3
n   <-  2.9
R0 <- Rmax * 0.03

Max_firing <- Rmax * Contrast^n / (Contrast^n + C50^n)
df <- data.frame(Contrast, Max_firing)

### Simulation ###
cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

n_neurons <- 180                                               # Number of neurons in the population
orientations <- seq(1, 180, by = 1)                            # Possible orientations (0 to 359 degrees)
preferred_orientations <-  seq(1, 180, length.out = n_neurons) # Preferred orientation of each neuron
max_firing_seq <- df$Max_firing                                # Maximum firing rate of each neuron
contrast <-       df$Contrast
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
   
    resp <- rnbinom(n_neurons, size = mu/0.3, mu = mu) + rnbinom(n_neurons, size = R0/0.3, mu = R0) # Fano factor = 1.2 (var = μ + μ^2/size)
   
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
    
    resp <- rnbinom(n_neurons, size = mu/0.3, mu = mu) + rnbinom(n_neurons, size = R0/0.3, mu = R0)
  
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

df_low <- do.call(rbind, all_results)
df_low$Contrast <- contrast[match(df_low$Contrast, max_firing_seq)]

### Contrast discrimination ###
df_low %>%
  group_by(Stimulus, Contrast, Trial) %>%
  summarise(Sum_spikes = sum(Spikes), .groups = "drop") -> df_sum_low

df_sum_low %>%
  group_by(Contrast) %>%
  summarise(Mean_sum_spikes = mean(Sum_spikes),
            SD_sum_spikes = sd(Sum_spikes), .groups = "drop") -> df_summary_low
colnames(df_summary_low) <- c("Contrast", "mu", "sigma")

dC <-  diff(df_summary_low$Contrast)
dmu <- diff(df_summary_low$mu)

# μ'
mu_prime <- c(
  dmu[1] / dC[1],                                 　　　　
  (dmu[-1]/dC[-1] + dmu[-length(dmu)]/dC[-length(dC)]) / 2,
  dmu[length(dmu)] / dC[length(dC)]                       
)
df_summary_low$mu_prime <- mu_prime

df_summary_low$DeltaC_pred <- df_summary_low$sigma / df_summary_low$mu_prime # JND
df_summary_low <- mutate(df_summary_low, Weber_ratio_pred = DeltaC_pred / Contrast) # ΔC / C (Weber ratio)

### Orientation discrimination (LFI) ###
df_low_wide <- df_low %>%
  pivot_wider(
    id_cols = c(Trial, Contrast, Stimulus),
    names_from = Neuron,
    values_from = Spikes)

df_low_wide$Stimulus_bin <- ifelse(df_low_wide$Stimulus == 90, 1, 0)

maxFR_values <- sort(unique(df_low_wide$Contrast))

lfi_results_low <- data.frame(
  Contrast = numeric(),
  LFI_per_neuron = numeric()
)

for (maxFR in maxFR_values) {
  df_sub <- df_low_wide[df_low_wide$Contrast == maxFR, ]
  lfi_value <- compute_lfi(df_sub)
  lfi_results_low <- rbind(lfi_results_low,
                           data.frame(Contrast = maxFR,
                                      LFI_per_neuron = lfi_value))
}

lfi_results_low$Contrast <- contrast

### Orientation discrimination (decoding with logistic regression) ###
decoding_results_low <- data.frame(
  Contrast = numeric(),
  Decoding_accuracy = numeric()
)

for (maxFR in maxFR_values) {
  df_sub <- df_low_wide[df_low_wide$Contrast == maxFR, ]
  accuracy <- run_logistic_cv_parallel(df_sub, k = 5, ncores = detectCores() - 1)
  decoding_results_low <- rbind(decoding_results_low,
                                data.frame(Contrast = maxFR,
                                           Decoding_accuracy = accuracy))
}

decoding_results_low <- mutate(decoding_results_low, D_prime = qnorm(Decoding_accuracy) - qnorm(1 - Decoding_accuracy),
                               Contrast = contrast)

### Fano factor sanity check ### 
df_high$Attention <- "high"
df_low$Attention  <- "low"
df_integrated <- rbind(df_high, df_low) 

df_integrated %>%
  group_by(Contrast, Stimulus, Neuron, Attention) %>%
  summarise(
    mean_spikes = mean(Spikes),
    var_spikes  = var(Spikes),
    fano = var_spikes / mean_spikes,
    .groups = "drop"
  ) %>%
  filter(mean_spikes > 0) %>%
  group_by(Contrast, Stimulus, Attention) %>%
  summarise(
    mean_mean_spikes = mean(mean_spikes, na.rm = TRUE),
    mean_var_spikes = mean(var_spikes),
    mean_fano = mean(fano, na.rm = TRUE), # average across neurons
    sd_fano   = sd(fano, na.rm = TRUE),   # sd across neurons
    .groups = "drop") -> ff
print(ff, n = 99999)

# write.csv(df_integrated, 
#           file = "samaha_effect_df_integrated.csv", 
#           row.names = FALSE)

### Figures ###
g1 <- ggplot(df, aes(x = Contrast, y = Max_firing)) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 130)) +
  ylab("Tuning curve peak") +
  theme_minimal()

g2 <- ggplot(df, aes(x = Contrast, y = Max_firing)) +
  geom_point(size = 2) +
  scale_x_log10() +
  ylab("Tuning curve peak") +
  theme_minimal()

g3 <- ggplot(df_sum_high, aes(x = Sum_spikes, color = factor(Contrast))) +
  geom_density(alpha = 0.2, linewidth = 1) +
  theme_minimal() +
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
  ggtitle("High attention:Fano factor = 1.0") +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

g4 <- ggplot(df_sum_low, aes(x = Sum_spikes, color = factor(Contrast))) +
  geom_density(alpha = 0.2, linewidth = 1) +
  theme_minimal() +
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
  ggtitle("Low attention: Fano factor = 1.3") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))

df_summary_high$Attention <- "high"
df_summary_low$Attention  <- "low"
df_summary <- rbind(df_summary_high, df_summary_low)

g5 <- ggplot(df_summary, aes(x = Contrast, y = DeltaC_pred, color = Attention)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  scale_x_log10() +
  scale_y_log10(limits = c(0.2, 45)) +
  ylab("ΔC (JND of d' = 1)") +
  xlab("Contrast") +
  theme(
    legend.position = c(0.4, 0.9),
    legend.justification = c(1, 1)
  )

g6 <- ggplot(df_summary, aes(x = Contrast, y = Weber_ratio_pred, color = Attention)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 1.5)) +
  ylab("ΔC/C (Weber fraction)") +
  xlab("Contrast") +
  theme(
    legend.position = c(0.4, 0.9),
    legend.justification = c(1, 1)
  )

g7 <- df_summary %>%
  ggplot(aes(x = Contrast, y = mu, color = Attention)) +
  geom_point(size = 2) +
  # scale_x_log10() +
  # scale_y_log10() +
  coord_cartesian(ylim = c(0, 7500)) +
  theme_minimal() +
  ylab("Mean total spikes") +
  theme(
    legend.position = c(0.9, 0.5),
    legend.justification = c(1, 1)
  )

lfi_results_high$Attention <- "high"
lfi_results_low$Attention  <- "low"
lfi_results <- rbind(lfi_results_high, lfi_results_low)

g8 <- lfi_results %>%
  ggplot(aes(x = Contrast, y = LFI_per_neuron, color = Attention)) +
  geom_point(size = 2) +
  scale_x_log10() +
  ylab("LFI per neuron for orientation") +
  theme_minimal() +
  theme(
    legend.position = c(0.4, 0.9),
    legend.justification = c(1, 1)
  )

decoding_results_high$Attention <- "high"
decoding_results_low$Attention  <- "low"
decoding_results <- rbind(decoding_results_high, decoding_results_low)

g9 <- decoding_results %>%
  ggplot(aes(x = Contrast, y = D_prime, color = Attention)) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(0, 100)) +
  ylab("d' for orientation") +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.5),
    legend.justification = c(1, 1)
  )

g10 <- decoding_results %>%
  ggplot(aes(x = Contrast, y = D_prime, color = Attention)) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(0, 100)) +
  scale_x_log10() +
  ylab("d' for orientation") +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.5),
    legend.justification = c(1, 1)
  )

df_sum_high$Attention <- "high"
df_sum_low$Attention  <- "low"
df_sum <- rbind(df_sum_high, df_sum_low)

df_sum %>%
  dplyr::filter(Contrast == 1 | Contrast == 4) %>%
  summarise(Criterion = mean(Sum_spikes)) %>%
  as.numeric() -> criterion

g11 <- ggplot(subset(df_sum, df_sum$Contrast == 1 | df_sum$Contrast == 4), 
            aes(x = Sum_spikes, color = factor(Contrast))) +
  geom_vline(xintercept = criterion, linetype = "dashed") +
  geom_density(alpha = 0.2, linewidth = 1) +
  theme_minimal() +
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
  coord_cartesian(xlim = c(500, 800), ylim = c(0, 0.02)) +
  facet_wrap(.~ Attention, nrow = 2)
g11
  
g12 <- df_sum %>%
  dplyr::filter(Contrast == 1 | Contrast == 4) %>%
  mutate(Yes = Sum_spikes > criterion) %>%
  group_by(Contrast, Attention) %>%
  summarise(P_yes = mean(Yes)) %>%
  ggplot(aes(x = Contrast, y = P_yes, color = factor(Contrast), shape = Attention)) +
  geom_point(size = 3) +
  ylim(0, 1) +
  theme_minimal() +
  labs(
    x = "Contrast",
    y = "P(Yes)",
    color = "Contrast",
    shape = "Attention"
  ) +
  theme(
    legend.position = c(0.95, 0.3),
    legend.justification = c(1, 0.5)
  ) +
  guides(color = "none")

p1 <- plot_grid(g4, g3,  ncol = 2)
p2 <- plot_grid(g1, g7,  ncol = 2)
p3 <- plot_grid(g5, g6,  ncol = 2)
p4 <- plot_grid(g10, g8, ncol = 2)
p5 <- plot_grid(g11, g12, ncol = 2)
subjective_inflation <- plot_grid(p1, p2, p3, p4, p5, ncol = 1, rel_heights = c(1, 1))
ggsave("subjective_inflation.png", subjective_inflation, width = 7, height = 11.25, dpi = 300)