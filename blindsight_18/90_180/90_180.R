# Simulation assuming Naka-Rushton function on Albrecht & Hamilton (1982, Fig. 4)

library(tidyverse)
library(doParallel)
library(foreach)
library(cowplot)
library(patchwork)

#################### High Gain Variability ####################
Contrast <- c(3, 4, 5, 6)
Rmax <- 115
C50 <- 19.3
n   <-  2.9
R0 <- Rmax * 0.03 # spontaneous firing based on Geisler & Albrecht (1997)

Max_firing <- Rmax * Contrast^n / (Contrast^n + C50^n)
df <- data.frame(Contrast, Max_firing)

### Simulation ###
cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

n_neurons <- 18                                               # Number of neurons in the population
orientations <- seq(1, 180, by = 1)                           # Possible orientations (1 to 180 degrees)
preferred_orientations <- seq(1, 180, length.out = n_neurons) # Preferred orientation of each neuron
max_firing_seq <- df$Max_firing                               # Maximum firing rate of each neuron
contrast <- df$Contrast
tuning_width <- 75                                            # Tuning width (standard deviation) of each neuron's response curve
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
    idx <- round(input_90[i]) 
    idx <- pmin(pmax(idx, 1), 180)
    
    sigma_g <- 0.2
    g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
    
    mu <- g * tuning_curves[, idx]
    mu[mu == 0] <- 1e-8
    
    resp <- rnbinom(n_neurons, size = 1000000, mu = mu) + rnbinom(n_neurons, size = 1000000, mu = R0)
    
    cbind(resp,
          neuron = seq_len(n_neurons),
          stim = 90,
          trial = i,
          maxFR = max_firing_rate)
  }
  
  # ---- 180 ----
  responses_180 <- foreach(
    i = 1:n_trials,
    .combine = rbind,
    .packages = "stats"
  ) %dopar% {
    
    input_180 <- rnorm(n_trials, 180, 0)
    idx <- round(input_180[i]) 
    idx <- pmin(pmax(idx, 1), 180)
    
    sigma_g <- 0.2
    g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
    
    mu <- g * tuning_curves[, idx]
    mu[mu == 0] <- 1e-8
    
    resp <- rnbinom(n_neurons, size = 1000000, mu = mu) + rnbinom(n_neurons, size = 1000000, mu = R0)
    
    cbind(resp,
          neuron = seq_len(n_neurons),
          stim = 180,
          trial = i,
          maxFR = max_firing_rate)
  }
  
  df_tmp <- as.data.frame(rbind(responses_90, responses_180))
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


#################### Low Gain Variability ####################
Contrast <- c(3, 4, 5, 6)
Rmax <- 115
C50 <- 19.3
n   <-  2.9
R0 <- Rmax * 0.03

Max_firing <- Rmax * Contrast^n / (Contrast^n + C50^n)
df <- data.frame(Contrast, Max_firing)

### Simulation ###
cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

n_neurons <- 18                                               # Number of neurons in the population
orientations <- seq(1, 180, by = 1)                           # Possible orientations (0 to 359 degrees)
preferred_orientations <- seq(1, 180, length.out = n_neurons) # Preferred orientation of each neuron
max_firing_seq <- df$Max_firing                               # Maximum firing rate of each neuron
contrast <- df$Contrast
tuning_width <- 75                                           # Tuning width (standard deviation) of each neuron's response curve
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
    idx <- round(input_90[i]) 
    idx <- pmin(pmax(idx, 1), 180)
    
    sigma_g <- 0.01
    g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
    
    mu <- g * tuning_curves[, idx]
    mu[mu == 0] <- 1e-8
    
    resp <- rnbinom(n_neurons, size = 1000000, mu = mu) + rnbinom(n_neurons, size = 1000000, mu = R0)
    
    cbind(resp,
          neuron = seq_len(n_neurons),
          stim = 90,
          trial = i,
          maxFR = max_firing_rate)
  }
  
  # ---- 180 ----
  responses_180 <- foreach(
    i = 1:n_trials,
    .combine = rbind,
    .packages = "stats"
  ) %dopar% {
    
    input_180 <- rnorm(n_trials, 180, 0)
    idx <- round(input_180[i]) 
    idx <- pmin(pmax(idx, 1), 180)
    
    sigma_g <- 0.01
    g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
    
    mu <- g * tuning_curves[, idx]
    mu[mu == 0] <- 1e-8
    
    resp <- rnbinom(n_neurons, size = 1000000, mu = mu) + rnbinom(n_neurons, size = 1000000, mu = R0)
    
    cbind(resp,
          neuron = seq_len(n_neurons),
          stim = 180,
          trial = i,
          maxFR = max_firing_rate)
  }
  
  df_tmp <- as.data.frame(rbind(responses_90, responses_180))
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
df_high$GV <- "high"
df_low$GV  <- "low"
df_integrated <- rbind(df_high, df_low) 

df_integrated %>%
  group_by(Contrast, Stimulus, Neuron, GV) %>%
  summarise(
    mean_spikes = mean(Spikes),
    var_spikes  = var(Spikes),
    fano = var_spikes / mean_spikes,
    .groups = "drop"
  ) %>%
  filter(mean_spikes > 0) %>%
  group_by(Contrast, Stimulus, GV) %>%
  summarise(
    mean_mean_spikes = mean(mean_spikes, na.rm = TRUE),
    mean_var_spikes = mean(var_spikes),
    mean_fano = mean(fano, na.rm = TRUE), # average across neurons
    sd_fano   = sd(fano, na.rm = TRUE),   # sd across neurons
    .groups = "drop") -> ff
print(ff, n = 99999)

cor_results <- df_integrated %>%
  mutate(Spikes = as.numeric(Spikes)) %>%
  group_by(GV, Stimulus, Contrast) %>%
  group_modify(~ {
    data_wide <- .x %>%
      select(Trial, Neuron, Spikes) %>%
      pivot_wider(names_from = Neuron, values_from = Spikes) %>%
      column_to_rownames("Trial")
    
    # 2. 【重要】すべての値が同じ（分散が0）列を特定して削除
    # 全トライアルで0スパイク、または発火数が一定のニューロンを除外
    keep_cols <- apply(data_wide, 2, function(x) sd(x, na.rm = TRUE) > 0)
    data_filtered <- data_wide[, keep_cols, drop = FALSE]
    
    # 有効なペアが残っていない場合は空のデータフレームを返す
    if (ncol(data_filtered) < 2) return(data.frame())
    
    cor_matrix <- cor(as.matrix(data_filtered), method = "pearson", use = "pairwise.complete.obs")
    
    cor_matrix %>%
      as.data.frame() %>%
      rownames_to_column("Neuron1") %>%
      pivot_longer(-Neuron1, names_to = "Neuron2", values_to = "Correlation")
  }) %>%
  ungroup() %>%
  filter(Neuron1 < Neuron2)

### Figures ###
g0 <- cor_results %>%
  ggplot() +
  geom_histogram(aes(x = Correlation, fill = GV)) +
  facet_wrap(. ~ Contrast) +
  theme_classic(base_size = 11) +
  coord_cartesian(xlim = c(0, 0.6)) +
  scale_x_continuous(breaks = c(0, 0.3, 0.6), labels = c("0", "0.3", "0.6")) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  labs(
    x = "Pairwise spike correlation",
    y = "Count"
  ) +
  theme(legend.position = c(0.9, 0.8))
g0

g1 <- ggplot(df, aes(x = Contrast, y = Max_firing)) +
  geom_point(size = 2.2, alpha = 0.85) +
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

g2 <- ggplot(df, aes(x = Contrast, y = Max_firing)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_x_log10(
    limits = c(1, 100),
    breaks = c(1, 3, 10, 30, 100)
  ) +
  coord_cartesian(ylim = c(0, 130)) +
  scale_y_continuous(breaks = seq(0, 120, by = 60)) +
  labs(
    x = "Contrast (%)",
    y = "Tuning curve peak"
  ) +
  theme_classic(base_size = 11)

df_sum_high$GV <- "high"
df_sum_low$GV  <- "low"
df_sum <- rbind(df_sum_high, df_sum_low)
df_sum %>%
  dplyr::filter(Contrast == 55 | Contrast == 70) %>%
  summarise(Criterion = mean(Sum_spikes)) %>%
  as.numeric() -> criterion

default_colors <- scales::hue_pal()(9)
two_colors <- default_colors[6:7]

df_sub <- subset(df_sum, Contrast %in% c(55, 70) & GV == "high")
df_sub$Contrast <- factor(df_sub$Contrast, levels = c(55, 70))

g_inset_high <- ggplot(df_sub, aes(x = Sum_spikes, color = factor(Contrast))) +
  geom_vline(xintercept = criterion, linetype = "dashed", linewidth = 0.3) +
  geom_density(alpha = 0.85, linewidth = 0.7) +
  scale_x_continuous(limits = c(5800, 6800), breaks = c(5800, 6800)) +
  scale_y_continuous(limits = c(0, 0.02), breaks = c(0, 0.01, 0.02)) + 
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_text(size = 6)
  ) +
  scale_color_manual(values = two_colors)

g3 <- ggplot(df_sum_high, aes(x = Sum_spikes, color = factor(Contrast))) +
  geom_density(alpha = 0.85, linewidth = 1) +
  coord_cartesian(xlim = c(0, 200), ylim = c(0.0015, 0.05)) +
  scale_x_continuous(breaks = c(0, 100, 200)) +
  theme_classic() +
  labs(
    x = "Total spikes",
    y = "Density",
    color = "Contrast"
  ) + 
  theme(
    legend.position = c(1.07, 1.05),
    legend.justification = c(1, 1),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8)
  ) +
  ggtitle("High gain variability:\nsigma_g = 0.2") +
  theme(plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5))

# g3 <- g3 + annotation_custom(
#   grob = ggplotGrob(g_inset_high),
#   xmin = 2000, xmax = 8000,
#   ymin = 0.01, ymax = 0.03
# )

legend_contrast <- get_legend(
  g3 + theme(legend.position = "right")
)

png("legend_contrast.png", width = 600, height = 400, res = 300)
grid::grid.draw(legend_contrast)
dev.off()

df_sub <- subset(df_sum, Contrast %in% c(55, 70) & GV == "low")
df_sub$Contrast <- factor(df_sub$Contrast, levels = c(55, 70))

g_inset_low <- ggplot(df_sub, aes(x = Sum_spikes, color = factor(Contrast))) +
  geom_vline(xintercept = criterion, linetype = "dashed", linewidth = 0.3) +
  geom_density(alpha = 0.85, linewidth = 0.7) +
  scale_x_continuous(limits = c(5800, 6800), breaks = c(5800, 6800)) +
  scale_y_continuous(limits = c(0, 0.02), breaks = c(0, 0.01, 0.02)) + 
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_text(size = 6)
  ) +
  scale_color_manual(values = two_colors)

g4 <- ggplot(df_sum_low, aes(x = Sum_spikes, color = factor(Contrast))) +
  geom_density(alpha = 0.85, linewidth = 1) +
  coord_cartesian(xlim = c(0, 200), ylim = c(0.0015, 0.05)) +
  scale_x_continuous(breaks = c(0, 100, 200)) +
  theme_classic() +
  labs(
    x = "Total spikes",
    y = "Density",
    color = "Contrast"
  ) + 
  theme(
    legend.position = c(1.07, 1.05),
    legend.justification = c(1, 1),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8)
  ) +
  ggtitle("Low gain variability:\nsigma_g = 0.01") +
  theme(plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5))

# g4 <- g4 + annotation_custom(
#   grob = ggplotGrob(g_inset_low),
#   xmin = 2000, xmax = 8000,
#   ymin = 0.01, ymax = 0.03
# )

df_summary_high$GV <- "high"
df_summary_low$GV  <- "low"
df_summary <- rbind(df_summary_high, df_summary_low)

g5 <- ggplot(df_summary, aes(x = Contrast, y = DeltaC_pred, color = GV)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_x_continuous(
    limits = c(0, 6),
    breaks = c(0, 2, 4, 6)
  ) +
  scale_y_log10(limits = c(0.2, 3)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  labs(
    x = "Contrast (%)",
    y = "ΔC (JND of d' = 1)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = c(0.45, 1),
    legend.justification = c(1, 1),
    legend.key = element_blank()
  )

legend_gv <- get_legend(
  g5 + theme(legend.position = "right")
)

png("legend_gv.png", width = 600, height = 400, res = 300)
grid::grid.draw(legend_gv)
dev.off()

g6 <- ggplot(df_summary, aes(x = Contrast, y = Weber_ratio_pred, color = GV)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_x_continuous(
    limits = c(0, 6),
    breaks = c(0, 2, 4, 6)
  ) +
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  labs(
    x = "Contrast (%)",
    y = "ΔC/C (Weber fraction)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.key = element_blank()
  )

g7 <- df_summary %>%
  ggplot(aes(x = Contrast, y = mu, color = GV)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_x_continuous(
    limits = c(0, 6),
    breaks = c(0, 2, 4, 6)
  ) +
  coord_cartesian(ylim = c(0, 200)) +
  scale_y_continuous(breaks = seq(0, 200, by = 100)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  labs(
    x = "Contrast (%)",
    y = "Mean total spikes"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = c(1, 0.55),
    legend.justification = c(1, 1),
    legend.key = element_blank()
  )

lfi_results_high$GV <- "high"
lfi_results_low$GV  <- "low"
lfi_results <- rbind(lfi_results_high, lfi_results_low)

g8 <- lfi_results %>%
  ggplot(aes(x = Contrast, y = LFI_per_neuron, color = GV)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_x_continuous(
    limits = c(0, 6),
    breaks = c(0, 2, 4, 6)
  ) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  labs(
    x = "Contrast (%)",
    y = "Orientation LFI"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = c(0.45, 1),
    legend.justification = c(1, 1),
    legend.key = element_blank()
  )

decoding_results_high$GV <- "high"
decoding_results_low$GV  <- "low"
decoding_results <- rbind(decoding_results_high, decoding_results_low)

g9 <- decoding_results %>%
  ggplot(aes(x = Contrast, y = D_prime, color = GV)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_x_continuous(
    limits = c(0, 6),
    breaks = c(0, 2, 4, 6)
  ) +
  scale_y_continuous(
    limits = c(-0.1, 3),
    breaks = seq(0, 3, by = 1)
  ) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  labs(
    x = "Contrast (%)",
    y = "Orientation d'"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = c(0.2, 0.95),
    legend.justification = c(1, 1),
    legend.key = element_blank()
  )

g10 <- decoding_results %>%
  ggplot(aes(x = Contrast, y = D_prime, color = GV)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_x_log10(
    limits = c(1, 10),
    breaks = c(1, 3, 10)
  ) +
  scale_y_continuous(
    limits = c(-0.1, 4),
    breaks = seq(0, 4, by = 1)
  ) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  labs(
    x = "Contrast (%)",
    y = "Orientation d'"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = c(0.4, 1),
    legend.justification = c(1, 1),
    legend.key = element_blank()
  )

plot_list <- list(g4, g3, g7, g5, g6, g9)
plots_no_legend <- lapply(plot_list, function(p) {
  p + theme(
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )
})

blindsight_18 <- wrap_plots(plots_no_legend, ncol = 3)
ggsave("90_180.png", blindsight_18, width = 6, height = 4, dpi = 300)


# LFI decomposition
compute_lfi_decomposition <- function(data, lambda = 1e-6) {
  
  neurons <- grep("^[0-9]+$", colnames(data), value = TRUE)
  p <- length(neurons)
  
  data0 <- data[data$Stimulus_bin == 0, neurons]
  data1 <- data[data$Stimulus_bin == 1, neurons]
  
  df <- colMeans(data1) - colMeans(data0)
  
  Sigma <- (cov(data0) + cov(data1)) / 2
  Sigma <- Sigma + lambda * diag(p)
  
  w <- solve(Sigma) %*% df
  
  contribution <- as.numeric(df * w)
  
  Sigma_diag <- diag(diag(Sigma))
  w_shuffled <- solve(Sigma_diag) %*% df
  I_shuffled <- as.numeric(t(df) %*% w_shuffled)
  
  model <- glm(Stimulus_bin ~ ., data = data[, c(neurons, "Stimulus_bin")], 
               family = binomial(link = "logit"))
  
  logistic_weights <- coef(model)[-1] 
  
  return(list(
    total_I = sum(contribution),
    neuron_contributions = contribution,
    shuffled_I = I_shuffled,
    correlation_gain = sum(contribution) / I_shuffled,
    lfi_weights = w,  # LFIによる重み
    logistic_weights = logistic_weights  # ロジスティック回帰による重み
  ))
}

maxFR_values <- c(3, 4, 5, 6)
figure_list <- list()

for (i in 1:length(maxFR_values)) {
  
  maxFR <- maxFR_values[i]
  df_sub_high <- df_high_wide[df_high_wide$Contrast == maxFR, ]
  lfi_high <- compute_lfi_decomposition(df_sub_high)
  
  plot_data_high <- data.frame(
    Neuron_ID = 1:length(lfi_high$neuron_contributions),
    Contribution = lfi_high$neuron_contributions
  )
  
  plot_logistic_high <- data.frame(
    Neuron_ID = 1:length(lfi_high$neuron_contributions),
    Weight = lfi_high$logistic_weights
  )
  
  df_sub_low <- df_low_wide[df_low_wide$Contrast == maxFR, ]
  lfi_low <- compute_lfi_decomposition(df_sub_low)
  
  plot_data_low <- data.frame(
    Neuron_ID = 1:length(lfi_low$neuron_contributions),
    Contribution = lfi_low$neuron_contributions
  )
  
  plot_logistic_low <- data.frame(
    Neuron_ID = 1:length(lfi_low$neuron_contributions),
    Weight = lfi_low$logistic_weights
  )
  
  combined_data <- data.frame(
    Neuron_ID = 1:length(lfi_high$neuron_contributions),
    High = lfi_high$neuron_contributions,
    Low = lfi_low$neuron_contributions
  ) %>%
    pivot_longer(
      cols = c(High, Low),
      names_to = "Condition",
      values_to = "Contribution"
    )
  
  gg1 <- ggplot(combined_data, aes(x = Neuron_ID, y = Contribution, color = Condition)) +
    geom_line(size = 1, alpha = 0.8) +
    scale_color_manual(values = c("High" = "#E41A1C", "Low" = "#377EB8")) +
    theme_classic() +
    labs(
      x = "Neuron ID (1 - 180)",
      y = "df * w",
      color = "Condition"
    ) +
    annotate("text", x = 1, y = Inf, hjust = 0, vjust = 1,
             label = paste0(
               "LFI:", round(lfi_high$total_I, 2),
               " | CG:", round(lfi_high$correlation_gain, 2)
             ),
             color = "#E41A1C", size = 4) +
    annotate("text", x = 1, y = Inf, hjust = 0, vjust = 3,
             label = paste0(
               "LFI:", round(lfi_low$total_I, 2),
               " | CG:", round(lfi_low$correlation_gain, 2)
             ),
             color = "#377EB8", size = 4) +
    scale_x_continuous(breaks = seq(0, 180, by = 30)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
    theme(legend.position = "top")
  
  combined_logistic <- data.frame(
    Neuron_ID = 1:length(lfi_high$logistic_weights),
    High = lfi_high$logistic_weights,
    Low = lfi_low$logistic_weights
  ) %>%
    pivot_longer(
      cols = c(High, Low),
      names_to = "Condition",
      values_to = "Weight"
    )
  
  gg2 <- ggplot(combined_logistic, aes(x = Neuron_ID, y = Weight, color = Condition)) +
    geom_line(size = 1, alpha = 0.8) +
    scale_color_manual(values = c("High" = "#E41A1C", "Low" = "#377EB8")) +
    theme_classic() +
    labs(
      x = "Neuron ID (1 - 180)",
      y = "Logistic weight",
      color = "Condition"
    ) +
    scale_x_continuous(breaks = seq(0, 180, by = 30)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
    theme(legend.position = "top")
  
  plot_list <- list(gg1, gg2)
  plots <- lapply(plot_list, function(p) {
    p + theme(
      legend.position = "none",
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
    )
  })
  
  plot <- wrap_plots(plots, row = 1)
  figure_list[[i]] <- plot
  
}

final_plot <- wrap_plots(figure_list, nrow = 2, ncol = 2) +
  plot_annotation(title = "90 vs. 180")
ggsave("lfi_90_180.png", final_plot, width = 9, height = 3.5, dpi = 300)
ggsave("cor_90_180.png", g0, width = 3.5, height = 3.5, dpi = 300)