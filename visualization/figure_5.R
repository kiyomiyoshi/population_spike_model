library(tidyverse)
library(doParallel)
library(foreach)
library(cowplot)
library(patchwork)
library(plotly)
library(magick)
library(webshot2)
library(minpack.lm)

set.seed(125)

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
  annotate(
    "text", 
    x = 50,              
    y = 135,             
    label = "Naka-Rushton", 
    vjust = 1,           
    hjust = 0.5,         
    size = 3.5             
  ) +
  scale_x_continuous(
    limits = c(1, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  coord_cartesian(ylim = c(0, 130)) +
  scale_y_continuous(breaks = seq(0, 120, by = 30)) +
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

# High gain variability
n_neurons <- 180                                              
orientations <- seq(1, 180, by = 1)                           
preferred_orientations <-  seq(1, 180, length.out = n_neurons)
max_firing_rate <- Rmax * 18^n / (18^n + C50^n) # 15% contrast
tuning_width    <- 20                                         
n_trials <- 30

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
  annotate(
    "text", 
    x = 90,              
    y = 135,             
    label = "Contrast = 15%", 
    vjust = 1,           
    hjust = 0.5,         
    size = 3.5             
  ) +
  theme_classic(base_size = 11) + 
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  coord_cartesian(ylim = c(0, 130)) +
  scale_y_continuous(breaks = seq(0, 120, by = 30)) +
  guides(color = F) + ylab("Spikes")
g2

input_90 <- rnorm(n_trials, 90, 0) # zero stimulus encoding noise
responses_90 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_90[i]) 
  sigma_g <- 0.2
  g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
  mu <- g * tuning_curves[, idx]
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
  sigma_g <- 0.2
  g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
  mu <- g * tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 100, trial = i)
}

df_responses_1 <- as.data.frame(rbind(responses_90, responses_100))
colnames(df_responses_1) <- c("Spikes", "Neuron", "Stimulus", "Trial")

# Low gain variability
n_neurons <- 180                                               
orientations <- seq(1, 180, by = 1)                            
preferred_orientations <-  seq(1, 180, length.out = n_neurons) 
max_firing_rate <- Rmax * 18^n / (18^n + C50^n)                                 
tuning_width <-    20                                  

tuning_curves <- matrix(0, nrow = n_neurons, ncol = length(orientations))
for (i in 1:n_neurons) {
  tuning_curves[i, ] <- max_firing_rate * exp(-0.5 * (pmin(abs(orientations - preferred_orientations[i]), 
                                                           180 - abs(orientations - preferred_orientations[i])) ^ 2) / tuning_width ^ 2)
}

input_90 <- rnorm(n_trials, 90, 0) # zero stimulus encoding noise
responses_90 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_90[i]) 
  sigma_g <- 0.01
  g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
  mu <- g * tuning_curves[, idx]
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
  sigma_g <- 0.01
  g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
  mu <- g * tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 100, trial = i)
}

df_responses_2 <- as.data.frame(rbind(responses_90, responses_100))
colnames(df_responses_2) <- c("Spikes", "Neuron", "Stimulus", "Trial")

# stopCluster(cl)

### Visualization ###
# across-trial spike distributions
# 2d density
df_responses_1$GV <- "high"
df_responses_2$GV <- "low"
df_integrated <- rbind(df_responses_1, df_responses_2)

g3 <- df_integrated %>%
  dplyr::filter(Neuron %in% c(90, 100)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `90`, y = `100`,
             color = GV)) +
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.5)) + # parametric density
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.7)) + 
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.9)) + 
  labs(x = "Neuron 90", y = "Neuron 100") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
g3

g4 <- df_integrated %>%
  filter(Neuron %in% c(90, 100)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `90`, y = `100`, color = GV)) +
  geom_density_2d(aes(linetype = factor(Stimulus)), linewidth = 0.7) + # nonparametric density
  labs(x = "Neuron 90", y = "Neuron 100",
       color = "GV", linetype = "Stimulus") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
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
g4

# sum of spikes
g5 <- df_integrated %>%
  group_by(Stimulus, GV, Trial) %>%
  summarise(Sum_spikes = sum(Spikes), .groups = "drop") %>%
  ggplot(aes(x = Sum_spikes, color = GV, fill = GV)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  coord_cartesian(xlim = c(0, 4000)) +
  scale_x_continuous(breaks = c(0, 1000, 2000, 3000, 4000)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic(base_size = 11) +
  labs(
    x = "Total spikes",
    y = "Density",
    color = "GV",
    fill = "GV") +
  theme(
    legend.position = c(1.1, 1),
    legend.justification = c(1, 1)
  )
g5

# Pairwise spike correlations
cor_results <- df_integrated %>%
  mutate(Spikes = as.numeric(Spikes)) %>%
  group_by(GV, Stimulus) %>%
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

print(head(cor_results))

g6 <- cor_results %>%
  ggplot() + geom_histogram(aes(x = Correlation, fill = GV)) +
  facet_wrap(. ~ GV) +
  theme_classic(base_size = 11) +
  coord_cartesian(xlim = c(0, 0.6)) +
  scale_x_continuous(breaks = c(0, 0.3, 0.6), labels = c("0", "0.3", "0.6")) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  theme(legend.position = "none") +
  labs(
    x = "Pairwise spike correlation",
    y = "Count")
g6

# 3d density
u <- seq(0, 2*pi, length.out = 50)
v <- seq(0, pi, length.out = 50)
x_sphere <- outer(cos(u), sin(v))
y_sphere <- outer(sin(u), sin(v))
z_sphere <- outer(rep(1,length(u)), cos(v)) 
unit_sphere <- rbind(as.vector(x_sphere), as.vector(y_sphere), as.vector(z_sphere))

df_integrated %>% 
  filter(Neuron %in% c(90, 95, 100)) %>% 
  pivot_wider(id_cols = c(GV, Stimulus, Trial), 
              names_from = Neuron, values_from = Spikes) %>%
  group_by(Stimulus, GV) %>% 
  group_map(~{
    mu <- colMeans(.x[,c("90","95","100")])
    sigma <- cov(.x[,c("90","95","100")])
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
      GV = .y$GV
    )
  }) -> ellipsoids

p1 <- plot_ly()

for (e in ellipsoids) {
  
  col <- ifelse(e$GV == 25, "Reds", "Blues")
  
  p1 <- p1 %>% add_surface(
    x = e$x,
    y = e$y,
    z = e$z,
    colorscale = col,
    showscale = FALSE,
    opacity = 0.6,
    name = paste("Stim:", e$Stimulus, "Con:", e$GV),
    hoverinfo = "text",
    text = paste("Stimulus:", e$Stimulus,
                 "<br>GV:", e$GV)
  )
}

p1 <- p1 %>%
  layout(
    scene = list(
      xaxis = list(title = "Neuron 90", range =  c(0, 100)),
      yaxis = list(title = "Neuron 95", range =  c(0, 100)),
      zaxis = list(title = "Neuron 100", range = c(0, 100))
    )
  )
p1

# 3d scatter plot
df_integrated %>% 
  filter(Neuron %in% c(90, 95, 100)) %>% 
  pivot_wider(id_cols = c(GV, Stimulus, Trial), 
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
    x = ~`90`,
    y = ~`95`,
    z = ~`100`,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 3,
      symbol = symbols[i]
    ),
    color = ~factor(GV),
    colors = colors,
    name = paste("Stimulus", stim_list[i])
  )
}

p2 <- p2 %>%
  layout(
    scene = list(
      xaxis = list(title = "Neuron 90",  range = c(0, 100)),
      yaxis = list(title = "Neuron 95",  range = c(0, 100)),
      zaxis = list(title = "Neuron 100", range = c(0, 100)),
      aspectmode = "cube"   
    )
  )

p2 <- p2 %>%
  layout(
    scene = list(
      xaxis = list(
        title = "Neuron 90",
        range = c(0, 100),
        tickvals = seq(0, 100, by = 20),
        color = "black", 
        tickcolor = "black",
        titlefont = list(color = "black")
      ),
      yaxis = list(
        title = "Neuron 95",
        range = c(0, 100),
        tickvals = seq(0, 100, by = 20),
        color = "black",
        tickcolor = "black",
        titlefont = list(color = "black")
      ),
      zaxis = list(
        title = "Neuron 100",
        range = c(0, 100),
        tickvals = seq(0, 100, by = 20),
        color = "black",
        tickcolor = "black",
        titlefont = list(color = "black")
      )
    )
  )

# discrimination plane
df_scatter$stim_bin <- as.numeric(as.factor(df_scatter$Stimulus)) - 1
fit <- glm(stim_bin ~ `90` + `95` + `100`, data = df_scatter, family = binomial)

x_seq <- seq(0, 100, length.out = 50)
y_seq <- seq(0, 100, length.out = 50)
z_seq <- seq(0, 100, length.out = 50)

grid3d <- expand.grid(`90` = x_seq, `95` = y_seq, `100` = z_seq)
grid3d$prob <- predict(fit, newdata = grid3d, type = "response")
decision_points <- grid3d %>% filter(abs(prob - 0.5) < 0.02)

p3 <- p2 %>%
  add_trace(
    data = decision_points,
    x = ~`90`, y = ~`95`, z = ~`100`,
    type = "mesh3d",
    color = I("navy"),    # <- I() ensures literal color
    opacity = 0.3,
    name = "Decision boundary",
    showscale = FALSE
  )

# visibility plane
z2 <- outer(x_seq, y_seq, function(x, y) 150 - x - y)

p4 <- p3 %>% add_surface(
  x = x_seq,
  y = y_seq,
  z = z2,
  opacity = 0.5,
  showscale = FALSE,
  name = "x+y+z=150 plane",
  surfacecolor = matrix(rep(1, length(x_seq)*length(y_seq)),
                        nrow = length(x_seq),
                        ncol = length(y_seq)),
  colorscale = list(c(0, 1), c("pink", "pink"))
)

# save 3d figure
p4 <- p4 %>%
  layout(
    font = list(size = 22),
    title = list(font = list(size = 28)),
    margin = list(l = 0, r = 0, b = 0, t = 40),
    scene = list(
      camera = list(eye = list(x = 1, y = -4, z = 2.5)),
      xaxis = list(titlefont = list(size = 16), tickfont = list(size = 12)),
      yaxis = list(titlefont = list(size = 16), tickfont = list(size = 12)),
      zaxis = list(titlefont = list(size = 16), tickfont = list(size = 12))
    )
  )
p4

saveWidget(p4, "figure_5.html", selfcontained = TRUE)
saveWidget(p4, "figure_5.png", selfcontained = TRUE)

webshot2::webshot(
  "figure_5.html",
  "figure_5.png",
  vwidth  = 800 * 1.1,
  vheight = 560 * 1.1,
  zoom = 10
)