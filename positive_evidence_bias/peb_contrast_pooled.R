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
library(mvtnorm)
library(viridis)
library(data.table)
library(mgcv)

cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

# Naka-Rushton function
Rmax <- 115
C50 <- 19.3
n   <-  2.9
R0 <- Rmax * 0.03

# parameters
n_neurons <- 180
orientations <- c(45, 135)
preferred_orientations <- seq(1, 180, length.out = n_neurons)
tuning_width <- 20
n_trials <- 5000

# stimulus 
contrasts <- seq(1, 10, 1)
stim_angles <- c(45, 135)

results_list <- list()

for (cc in contrasts) {
  
  max_firing_rate <- Rmax * cc^n / (cc^n + C50^n)
  
  tuning_curves <- matrix(0, nrow = n_neurons, ncol = length(orientations))
  
  for (i in 1:n_neurons) {
    tuning_curves[i, ] <- max_firing_rate * exp(
      -0.5 * (pmin(abs(orientations - preferred_orientations[i]),
                   180 - abs(orientations - preferred_orientations[i]))^2) /
        tuning_width^2
    )
  }
  
  for (stim_angle in stim_angles) {
    
    input_stim <- rnorm(n_trials, stim_angle, 0)
    
    responses <- foreach(i = 1:n_trials,
                         .combine = rbind,
                         .packages = "stats") %dopar% {
                           idx <- ifelse(round(input_stim[i]) == 45, 1, 2)
                           
                           g <- 1
                           mu <- g * tuning_curves[, idx]
                           mu[mu == 0] <- 1e-8
                           
                           resp <- rnbinom(n_neurons, size = 1000000, mu = mu) +
                             rnbinom(n_neurons, size = 1000000, mu = R0)
                           
                           cbind(resp,
                                 neuron = seq_len(n_neurons),
                                 stim = stim_angle,
                                 contrast = cc,
                                 trial = i)
                         }
    key <- paste(cc, stim_angle, sep = "_")
    results_list[[key]] <- as.data.frame(responses)
  }
}

df_responses <- do.call(rbind, results_list)
colnames(df_responses) <- c("Spikes", "Neuron", "Stimulus", "Contrast", "Trial")
df_responses <- df_responses %>%
  mutate(Target = ifelse(Stimulus == 45, 1, 0))
colSums(is.na(df_responses))
any(is.na(df_responses))
# fwrite(df_responses, "df_responses.csv")
# df_responses <- fread("df_responses.csv", header = T)

# LFI
compute_lfi_decomposition <- function(data, lambda = 1e-6) {
  
  neurons <- grep("^[0-9]+$", colnames(data), value = TRUE)
  p <- length(neurons)
  
  data0 <- data[data$Target == 0, neurons]
  data1 <- data[data$Target == 1, neurons]
  
  df <- colMeans(data1) - colMeans(data0)
  
  Sigma <- (cov(data0) + cov(data1)) / 2
  Sigma <- Sigma + lambda * diag(p)
  
  # 1. 最適なデコーディング重み w = Sigma^-1 * df
  w <- solve(Sigma) %*% df
  
  # 2. 各ニューロンの貢献度
  contribution <- as.numeric(df * w)
  
  # 3. 独立仮定時の情報量
  Sigma_diag <- diag(diag(Sigma))
  w_shuffled <- solve(Sigma_diag) %*% df
  I_shuffled <- as.numeric(t(df) %*% w_shuffled)
  
  return(list(
    total_I = sum(contribution),
    neuron_contributions = contribution,
    shuffled_I = I_shuffled,
    correlation_gain = sum(contribution) / I_shuffled,
    lfi_weights = w
    ))
}

df_wide <- df_responses %>%
  pivot_wider(
    id_cols = c(Stimulus, Contrast, Trial, Target),
    names_from = Neuron,
    values_from = Spikes
  )
colSums(is.na(df_wide))
any(is.na(df_wide))

lfi <- compute_lfi_decomposition(df_wide)
lfi
plot(cumsum(sort(lfi$neuron_contributions, decreasing = TRUE)))

plot_data <- data.frame(
  Neuron_ID = 1:length(lfi$neuron_contributions),
  Contribution = lfi$neuron_contributions
)

g1 <- ggplot(plot_data, aes(x = Neuron_ID, y = Contribution)) +
  geom_bar(stat = "identity", fill = "#E41A1C") +
  theme_classic() +
  labs(
    subtitle = paste("Total LFI:", round(lfi$total_I, 2), 
                     " | Correlation Gain:", round(lfi$correlation_gain, 2)),
    x = "Neuron ID (1 - 180)",
    y = "Contribution (df * w)"
  ) +
  scale_x_continuous(breaks = seq(0, 180, by = 45)) +
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) 
g1

# across-trial spike distributions
# 2d density
df <- df_responses

g2 <- df %>%
  dplyr::filter(Neuron %in% c(45, 135)) %>%
  tidyr::pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `45`, y = `135`, colour = factor(Target))) +
  stat_ellipse(level = 0.5) +
  stat_ellipse(level = 0.7) +
  stat_ellipse(level = 0.9) +
  scale_x_continuous(limits = c(0, 12), breaks = seq(0, 12, 4)) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 4)) +
  scale_colour_discrete(
    name = "Stimulus",
    labels = c("135 deg", "45 deg")
  ) +
  labs(
    x = "Neuron 45",
    y = "Neuron 135"
  ) +
  coord_fixed() +
  theme_classic(base_size = 11) +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_blank()
  )
g2

# proportion of 45 deg stimulus
g3 <- df %>%
  dplyr::filter(Neuron %in% c(45, 135)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `45`, y = `135`, z = Target)) +
  stat_summary_2d(
    fun = mean,
    bins = 15
  ) +
  scale_fill_viridis_c(
    limits = c(0, 1)
  ) +
  coord_equal() +
  labs(
    fill = "Proportion(45 deg)",
    x = "Neuron 45",
    y = "Neuron 135"
  ) +
  theme_minimal()
g3

# confidence as posterior probability
params <- df %>%
  filter(Neuron %in% c(45, 135)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  group_by(Stimulus) %>%
  summarise(
    mu_x = mean(`45`,  na.rm = TRUE),
    mu_y = mean(`135`, na.rm = TRUE),
    cov_mat = list(cov(cbind(`45`, `135`), use = "complete.obs")),
    .groups = "drop"
  )
params$cov_mat[[1]]
eigen(params$cov_mat[[1]])

sigma_45 <-  sqrt(params$cov_mat[[1]][1, 1])
sigma_135 <- sqrt(params$cov_mat[[1]][2, 2])
correlation <- params$cov_mat[[1]][1, 2] / (sigma_45 * sigma_135)
correlation
mu_45  <- params[1, 2:3]
mu_135 <- params[2, 2:3]
cov_mat_45  <- params$cov_mat[[1]]
cov_mat_135 <- params$cov_mat[[2]]

likelihood <- function(x, y, mu_x, mu_y, cov_mat) {
  return(dmvnorm(cbind(x, y), mean = c(mu_x, mu_y), sigma = cov_mat))
}

prior_45  <- 0.5
prior_135 <- 0.5

x_range <- seq(0, 12, length.out = 100)
y_range <- seq(0, 12, length.out = 100)

grid <- expand.grid(x = x_range, y = y_range)
grid$likelihood_45   <- mapply(likelihood, grid$x, grid$y, MoreArgs = list(mu_x = as.numeric(mu_45[1, 1]),  mu_y = as.numeric(mu_45[1, 2]),  cov_mat = cov_mat_45))
grid$likelihood_135  <- mapply(likelihood, grid$x, grid$y, MoreArgs = list(mu_x = as.numeric(mu_135[1, 1]), mu_y = as.numeric(mu_135[1, 2]), cov_mat = cov_mat_135))
grid$posterior_45 <- (grid$likelihood_45 * prior_45) / (grid$likelihood_45 * prior_45 + grid$likelihood_135 * prior_135)

g4 <- ggplot(grid, aes(x = x, y = y, fill = posterior_45)) +
  geom_tile() +
  scale_fill_viridis(
    option = "C",
    limits = c(0, 1)
  ) +
  labs(
    x = "Neuron 45",
    y = "Neuron 135",
    fill = "P(Stimulus = 45 | x)"
  ) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 12), breaks = seq(0, 12, 4)) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 4)) +
  theme(legend.position = "right")
g4

# pairwise spike correlations
cor_results_1 <- df %>%
  mutate(Spikes = as.numeric(Spikes)) %>%
  group_by(Stimulus, Contrast) %>%
  group_modify(~ {
    data_wide <- .x %>%
      dplyr::select(Trial, Neuron, Spikes) %>%
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

g5 <- cor_results_1 %>%
  ggplot() +
  geom_histogram(aes(x = Correlation)) +
  facet_wrap(. ~ Contrast) +
  theme_classic(base_size = 11) +
  labs(
    x = "Pairwise spike correlation",
    y = "Count"
  ) +
  theme(legend.position = c(0.9, 0.8))
g5

# 3d density
u <- seq(0, 2*pi, length.out = 50)
v <- seq(0, pi, length.out = 50)
x_sphere <- outer(cos(u), sin(v))
y_sphere <- outer(sin(u), sin(v))
z_sphere <- outer(rep(1,length(u)), cos(v)) 
unit_sphere <- rbind(as.vector(x_sphere), as.vector(y_sphere), as.vector(z_sphere))

df %>% 
  filter(Neuron %in% c(45, 90, 135)) %>% 
  pivot_wider(id_cols = c(Stimulus, Contrast, Trial), 
              names_from = Neuron, values_from = Spikes) %>%
  group_by(Stimulus) %>% 
  group_map(~{
    mu <- colMeans(.x[,c("45","90","135")])
    sigma <- cov(.x[,c("45","90","135")])
    if(det(sigma) <= 0){
      sigma <- sigma + diag(1e-6,3)
    }
    eig <- eigen(sigma)
    coords <- eig$vectors %*% diag(sqrt(eig$values)) %*% unit_sphere + mu
    list(
      x = matrix(coords[1,], nrow=length(u)),
      y = matrix(coords[2,], nrow=length(u)),
      z = matrix(coords[3,], nrow=length(u)),
      Stimulus = .y$Stimulus
    )
  }) -> ellipsoids

p1 <- plot_ly()

for (e in ellipsoids) {
  
  col <- "#E41A1C" # "#377EB8"
  
  p1 <- p1 %>% add_surface(
    x = e$x,
    y = e$y,
    z = e$z,
    colorscale = list(c(0,"white"), c(1,col)),
    showscale = FALSE,
    opacity = 0.6,
    name = paste("Stim:", e$Stimulus),
    hoverinfo = "text",
    text = paste("Stimulus:", e$Stimulus)
  )
}

p1 <- p1 %>%
  layout(
    scene = list(
      xaxis = list(title = "Neuron 45",  range = c(0, 15)),
      yaxis = list(title = "Neuron 90",  range = c(0, 15)),
      zaxis = list(title = "Neuron 135", range = c(0, 15)),
      aspectmode = "cube"   
    )
  )

p1

# 3d scatter plot
df %>% 
  filter(Neuron %in% c(45, 90, 135)) %>% 
  pivot_wider(id_cols = c(Stimulus, Contrast, Trial), 
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
    x = ~`45`,
    y = ~`90`,
    z = ~`135`,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 3,
      symbol = symbols[i]
    ),
    # color = ~factor(GV),
    colors = colors,
    name = paste("Stimulus", stim_list[i])
  )
}

p2 <- p2 %>%
  layout(
    scene = list(
      xaxis = list(title = "Neuron 45",  range = c(0, 40)),
      yaxis = list(title = "Neuron 90",  range = c(0, 40)),
      zaxis = list(title = "Neuron 135", range = c(0, 40)),
      aspectmode = "cube"   
    )
  )

p2 <- p2 %>%
  layout(
    scene = list(
      xaxis = list(
        title = "Neuron 45",
        range = c(0, 40),
        color = "black", 
        tickcolor = "black",
        titlefont = list(color = "black")
      ),
      yaxis = list(
        title = "Neuron 90",
        range = c(0, 40),
        color = "black",
        tickcolor = "black",
        titlefont = list(color = "black")
      ),
      zaxis = list(
        title = "Neuron 135",
        range = c(0, 40),
        color = "black",
        tickcolor = "black",
        titlefont = list(color = "black")
      )
    )
  )

p2

# 3d scatter plot with logistic regression
df_scatter$stim_bin <- as.numeric(as.factor(df_scatter$Stimulus)) - 1
fit <- glm(stim_bin ~ `45` + `90` + `135`, data = df_scatter, family = binomial)

x_seq <- seq(0, 40, length.out = 50)
y_seq <- seq(0, 40, length.out = 50)
z_seq <- seq(0, 40, length.out = 50)

grid3d <- expand.grid(`45` = x_seq, `90` = y_seq, `135` = z_seq)
grid3d$prob <- predict(fit, newdata = grid3d, type = "response")
decision_points <- grid3d %>% filter(abs(prob - 0.5) < 0.02)

p3 <- p2 %>%
  add_trace(
    data = decision_points,
    x = ~`45`, y = ~`90`, z = ~`135`,
    type = "mesh3d",
    color = I("navy"),    # <- I() ensures literal color
    opacity = 0.3,
    name = "Decision boundary",
    showscale = FALSE
  )

p3

# save figures
plot_list <- list(g1, g2, g3, g4)
plots <- lapply(plot_list, function(p) {
  p + theme(
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )
})
peb_contrast_pooled <- wrap_plots(plots, ncol = 4)
ggsave("peb_contrast_pooled.png", peb_contrast_pooled, width = 12, height = 3, dpi = 2000)

p2 <- p2 %>%
  plotly::layout(
    scene = list(
      camera = list(
        eye = list(x = 2.2, y = 2.2, z = 1.2)
      )
    )
  )

p3 <- p3 %>%
  plotly::layout(
    scene = list(
      camera = list(
        eye = list(x = 2.2, y = 2.2, z = 1.2)
      )
    )
  )

saveWidget(p2, "scatter_3d_peb_cp.html", selfcontained = TRUE)
saveWidget(p3, "scatter_3d_peb_logistic_cp.html", selfcontained = TRUE)