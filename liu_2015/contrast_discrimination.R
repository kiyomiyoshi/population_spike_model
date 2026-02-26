# Simulation assuming log squared function on Liu et al. (2015, fig. 1A)

library(tidyverse)
library(doParallel)
library(foreach)
library(cowplot)

Contrast <- c(2, 12, 22, 32, 42, 52, 62, 72)
a <- 2.2766
b <- 1.2853

Max_firing <- a * (log(Contrast) + b)^2
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

### Weber law in contrast discrimination ###
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

g4 <- ggplot(df_plot, aes(x = Contrast, y = DeltaC_pred)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  ylab("ΔC (JND of d' = 1)") +
  xlab("Contrast")

g5 <- ggplot(df_plot, aes(x = Contrast, y = Weber_ratio_pred)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  coord_cartesian(ylim = c(0.03, 0.05)) +
  ylab("ΔC/C") +
  xlab("Contrast")

g6 <- df_summary %>%
  ggplot(aes(x = Contrast, y = mu)) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal()

top    <- plot_grid(g3, ncol = 1)
middle <- plot_grid(g1, g6, ncol = 2)
bottom <- plot_grid(g4, g5, ncol = 2)
g <- plot_grid(top, middle, bottom, ncol = 1, rel_heights = c(1, 1))
ggsave("contrast_discrimination_log_squared.png", g, width = 6, height = 8.2, dpi = 300)