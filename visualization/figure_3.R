library(tidyverse)
library(doParallel)
library(foreach)
library(cowplot)
library(patchwork)
library(plotly)
library(magick)
library(htmlwidgets)
library(webshot2)
library(minpack.lm)
library(ggpp)

#################### High Alpha ####################
Contrast <- 25
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
orientations <- seq(1, 180, by = 1)                           # Possible orientations (1 to 180 degrees)
preferred_orientations <- seq(1, 180, length.out = n_neurons) # Preferred orientation of each neuron
max_firing_seq <- df$Max_firing                               # Maximum firing rate of each neuron
contrast <- df$Contrast
tuning_width <- 20                                            # Tuning width (standard deviation) of each neuron's response curve
n_trials <- 30

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
    
    mu <- tuning_curves[, idx]
    mu[mu == 0] <- 1e-8
    
    resp <- rnbinom(n_neurons, size = 1000000, mu = mu) + rnbinom(n_neurons, size = 1000000, mu = R0)
    
    cbind(resp,
          neuron = seq_len(n_neurons),
          stim = 90,
          trial = i,
          maxFR = max_firing_rate)
  }
  
  # ---- 100 ----
  responses_100 <- foreach(
    i = 1:n_trials,
    .combine = rbind,
    .packages = "stats"
  ) %dopar% {
    
    input_100 <- rnorm(n_trials, 100, 0)
    idx <- round(input_100[i]) 
    idx <- pmin(pmax(idx, 1), 180)
    
    mu <- tuning_curves[, idx]
    mu[mu == 0] <- 1e-8
    
    resp <- rnbinom(n_neurons, size = 1000000, mu = mu) + rnbinom(n_neurons, size = 1000000, mu = R0)
    
    cbind(resp,
          neuron = seq_len(n_neurons),
          stim = 100,
          trial = i,
          maxFR = max_firing_rate)
  }
  
  df_tmp <- as.data.frame(rbind(responses_90, responses_100))
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


#################### Low Alpha ####################
Contrast <- 25
Rmax <- 115
C50 <- 19.3
n   <-  2.9
R0 <- Rmax * 0.13

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
n_trials <- 30

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
    
    mu <- tuning_curves[, idx]
    mu[mu == 0] <- 1e-8
    
    resp <- rnbinom(n_neurons, size = 1000000, mu = mu) + rnbinom(n_neurons, size = 1000000, mu = R0)
    
    cbind(resp,
          neuron = seq_len(n_neurons),
          stim = 90,
          trial = i,
          maxFR = max_firing_rate)
  }
  
  # ---- 100 ----
  responses_100 <- foreach(
    i = 1:n_trials,
    .combine = rbind,
    .packages = "stats"
  ) %dopar% {
    
    input_100 <- rnorm(n_trials, 100, 0)
    idx <- round(input_100[i]) 
    idx <- pmin(pmax(idx, 1), 180)
    
    mu <- tuning_curves[, idx]
    mu[mu == 0] <- 1e-8
    
    resp <- rnbinom(n_neurons, size = 1000000, mu = mu) + rnbinom(n_neurons, size = 1000000, mu = R0)
    
    cbind(resp,
          neuron = seq_len(n_neurons),
          stim = 100,
          trial = i,
          maxFR = max_firing_rate)
  }
  
  df_tmp <- as.data.frame(rbind(responses_90, responses_100))
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
df_high$Alpha <- "high"
df_low$Alpha  <- "low"
df_integrated <- rbind(df_high, df_low) 

df_integrated %>%
  group_by(Contrast, Stimulus, Neuron, Alpha) %>%
  summarise(
    mean_spikes = mean(Spikes),
    var_spikes  = var(Spikes),
    fano = var_spikes / mean_spikes,
    .groups = "drop"
  ) %>%
  filter(mean_spikes > 0) %>%
  group_by(Contrast, Stimulus, Alpha) %>%
  summarise(
    mean_mean_spikes = mean(mean_spikes, na.rm = TRUE),
    mean_var_spikes = mean(var_spikes),
    mean_fano = mean(fano, na.rm = TRUE), # average across neurons
    sd_fano   = sd(fano, na.rm = TRUE),   # sd across neurons
    .groups = "drop") -> ff
print(ff, n = 99999)

### Visualization ###
# 2d density
g1 <- df_integrated %>%
  dplyr::filter(Neuron %in% c(90, 100)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `90`, y = `100`,
             color = Alpha)) +
  stat_ellipse(aes(group = interaction(Alpha, Stimulus)), level = c(0.5)) + # parametric density
  stat_ellipse(aes(group = interaction(Alpha, Stimulus)), level = c(0.7)) + 
  stat_ellipse(aes(group = interaction(Alpha, Stimulus)), level = c(0.9)) + 
  labs(x = "Neuron 90", y = "Neuron 100") +
  coord_cartesian(xlim = c(0, 150), ylim = c(0, 150)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = c(0.3, 1),
    legend.justification = c(1, 1)
  )

g2 <- df_integrated %>%
  filter(Neuron %in% c(90, 100)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `90`, y = `100`, color = Alpha)) +
  geom_density_2d(aes(linetype = factor(Stimulus)), linewidth = 0.7) + # nonparametric density
  labs(x = "Neuron 90", y = "Neuron 100",
       color = "Alpha", linetype = "Stimulus") +
  coord_cartesian(xlim = c(0, 150), ylim = c(0, 150)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = c(0.3, 1),
    legend.justification = c(1, 1)
  ) +
  theme(
    legend.position = c(0.3, 1),
    legend.justification = c(1, 1)
  )

# sum of spikes
g3 <- df_integrated %>%
  group_by(Stimulus, Alpha, Trial) %>%
  summarise(Sum_spikes = sum(Spikes), .groups = "drop") %>%
  ggplot(aes(x = Sum_spikes, color = Alpha, fill = Alpha)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  coord_cartesian(xlim = c(4000, 7000)) +
  scale_x_continuous(breaks = c(4000, 5000, 6000, 7000)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic(base_size = 11) +
  labs(
    x = "Total spikes",
    y = "Density",
    color = "Alpha",
    fill = "Alpha") +
  theme(
    legend.position = c(0.15, 1),
    legend.justification = c(1, 1)
  )

# 3d density
u <- seq(0, 2*pi, length.out = 50)
v <- seq(0, pi, length.out = 50)
x_sphere <- outer(cos(u), sin(v))
y_sphere <- outer(sin(u), sin(v))
z_sphere <- outer(rep(1,length(u)), cos(v)) 
unit_sphere <- rbind(as.vector(x_sphere), as.vector(y_sphere), as.vector(z_sphere))

df_integrated %>% 
  filter(Neuron %in% c(80, 90, 100)) %>% 
  pivot_wider(id_cols = c(Alpha, Stimulus, Trial), 
              names_from = Neuron, values_from = Spikes) %>%
  group_by(Stimulus, Alpha) %>% 
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
      Alpha = .y$Alpha
    )
  }) -> ellipsoids

p1 <- plot_ly()

for (e in ellipsoids) {
  
  col <- ifelse(e$Alpha == 25, "Reds", "Blues")
  
  p1 <- p1 %>% add_surface(
    x = e$x,
    y = e$y,
    z = e$z,
    colorscale = col,
    showscale = FALSE,
    opacity = 0.6,
    name = paste("Stim:", e$Stimulus, "Con:", e$Alpha),
    hoverinfo = "text",
    text = paste("Stimulus:", e$Stimulus,
                 "<br>Alpha:", e$Alpha)
  )
}

p1 <- p1 %>%
  layout(
    scene = list(
      xaxis = list(title = "Neuron 80", range =  c(40, 120)),
      yaxis = list(title = "Neuron 90", range =  c(40, 120)),
      zaxis = list(title = "Neuron 100", range = c(40, 120))
    )
  )
p1

# 3d scatter plot
df_integrated %>% 
  filter(Neuron %in% c(80, 90, 100)) %>% 
  pivot_wider(id_cols = c(Alpha, Stimulus, Trial), 
              names_from = Neuron, values_from = Spikes) -> df_scatter

stim_list <- unique(df_scatter$Stimulus)
symbols <- c("circle","cross")
colors <- c("#E41A1C", "#377EB8")
# colors <- scales::hue_pal()(9)[6:7]

p2 <- plot_ly()

for (i in seq_along(stim_list)) {
  d <- df_scatter %>% filter(Stimulus == stim_list[i])
  p2 <- p2 %>% add_trace(
    data = d,
    x = ~`80`,
    y = ~`90`,
    z = ~`100`,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 3,
      symbol = symbols[i]
    ),
    color = ~factor(Alpha),
    colors = colors,
    name = paste("Stimulus", stim_list[i])
  )
}

p2 <- p2 %>%
  layout(
    scene = list(
      xaxis = list(title = "Neuron 80",  range = c(40, 120)),
      yaxis = list(title = "Neuron 90",  range = c(40, 120)),
      zaxis = list(title = "Neuron 100", range = c(40, 120)),
      aspectmode = "cube"   
    )
  )

p2 <- p2 %>%
  layout(
    scene = list(
      xaxis = list(
        title = "Neuron 80",
        range = c(40, 120),
        tickvals = seq(40, 120, by = 20),
        color = "black", 
        tickcolor = "black",
        titlefont = list(color = "black")
      ),
      yaxis = list(
        title = "Neuron 90",
        range = c(40, 120),
        tickvals = seq(40, 120, by = 20),
        color = "black",
        tickcolor = "black",
        titlefont = list(color = "black")
      ),
      zaxis = list(
        title = "Neuron 100",
        range = c(40, 120),
        tickvals = seq(40, 120, by = 20),
        color = "black",
        tickcolor = "black",
        titlefont = list(color = "black")
      )
    )
  )
p2

# discrimination plane
df_scatter$stim_bin <- as.numeric(as.factor(df_scatter$Stimulus)) - 1
fit <- glm(stim_bin ~ `80` + `90` + `100`, data = df_scatter, family = binomial)

x_seq <- seq(40, 120, length.out = 50)
y_seq <- seq(40, 120, length.out = 50)
z_seq <- seq(40, 120, length.out = 50)

grid3d <- expand.grid(`80` = x_seq, `90` = y_seq, `100` = z_seq)
grid3d$prob <- predict(fit, newdata = grid3d, type = "response")
decision_points <- grid3d %>% filter(abs(prob - 0.5) < 0.02)

p3 <- p2 %>%
  add_trace(
    data = decision_points,
    x = ~`80`, y = ~`90`, z = ~`100`,
    type = "mesh3d",
    color = I("navy"),    # <- I() ensures literal color
    opacity = 0.3,
    name = "Decision boundary",
    showscale = FALSE
  )
p3

# visibility plane
z2 <- outer(x_seq, y_seq, function(x, y) 230 - x - y)

p4 <- p3 %>% add_surface(
  x = x_seq,
  y = y_seq,
  z = z2,
  opacity = 0.3,
  showscale = FALSE,
  name = "x+y+z=230 plane",
  surfacecolor = matrix(rep(1, length(x_seq)*length(y_seq)),
                        nrow = length(x_seq),
                        ncol = length(y_seq)),
  colorscale = list(c(0, 1), c("pink", "pink"))
)
p4

# save 3d figure
p4 <- p4 %>%
  layout(
    font = list(size = 22),
    title = list(font = list(size = 28)),
    margin = list(l = 0, r = 0, b = 0, t = 40),
    scene = list(
      camera = list(eye = list(x = 3, y = -3, z = 2.2)),
      xaxis = list(titlefont = list(size = 16), tickfont = list(size = 12)),
      yaxis = list(titlefont = list(size = 16), tickfont = list(size = 12)),
      zaxis = list(titlefont = list(size = 16), tickfont = list(size = 12))
    )
  )

saveWidget(p4, "p4.html", selfcontained = TRUE)
saveWidget(p4, "p4.png", selfcontained = TRUE)

webshot2::webshot(
  "p4.html",
  "p4.png",
  vwidth  = 800,
  vheight = 560,
  zoom = 10
)