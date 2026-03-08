library(tidyverse)
library(doParallel)
library(foreach)
library(cowplot)
library(patchwork)
library(plotly)
library(magick)
library(webshot2)
library(minpack.lm)

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

# High contrast
n_neurons <- 180                                              
orientations <- seq(1, 180, by = 1)                           
preferred_orientations <-  seq(1, 180, length.out = n_neurons)
max_firing_rate <- Rmax * 25^n / (25^n + C50^n) # 25% contrast
tuning_width    <- 20                                         
n_trials <- 3000

tuning_curves <- matrix(0, nrow = n_neurons, ncol = length(orientations))
for (i in 1:n_neurons) {
  tuning_curves[i, ] <- max_firing_rate * exp(-0.5 * (pmin(abs(orientations - preferred_orientations[i]), 
                                                           180 - abs(orientations - preferred_orientations[i])) ^ 2) / tuning_width ^ 2)
}

input_orientation <- 135
responses_high    <- rpois(n_neurons, lambda = tuning_curves[, input_orientation])
df_responses_high <- data.frame(Spikes = responses_high, Neuron = seq(1:n_neurons))

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
    label = "Contrast = 25%", 
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

input_45 <- rnorm(n_trials, 45, 0) # zero stimulus encoding noise
responses_45 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_45[i]) 
  mu <- tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 45, trial = i)
}

input_135 <- rnorm(n_trials, 135, 0)
responses_135 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_135[i]) 
  mu <- tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 135, trial = i)
}

df_responses_1 <- as.data.frame(rbind(responses_45, responses_135))
colnames(df_responses_1) <- c("Spikes", "Neuron", "Stimulus", "Trial")

# Low contrast
n_neurons <- 180                                               
orientations <- seq(1, 180, by = 1)                            
preferred_orientations <-  seq(1, 180, length.out = n_neurons) 
max_firing_rate <- Rmax * 15^n / (15^n + C50^n) # 15% contrast                                 
tuning_width <-    20                                     

tuning_curves <- matrix(0, nrow = n_neurons, ncol = length(orientations))
for (i in 1:n_neurons) {
  tuning_curves[i, ] <- max_firing_rate * exp(-0.5 * (pmin(abs(orientations - preferred_orientations[i]), 
                                                           180 - abs(orientations - preferred_orientations[i])) ^ 2) / tuning_width ^ 2)
}

input_orientation <- 135
responses_low    <- rpois(n_neurons, lambda = tuning_curves[, input_orientation])
df_responses_low <- data.frame(Spikes = responses_low, Neuron = seq(1:n_neurons))

FiringRate <- c()
for (i in 1:n_neurons) {
  FiringRate <- c(FiringRate, tuning_curves[i, ])
}

df_tuning <- data.frame(Orientation = rep(orientations, n_neurons),
                        FiringRate = FiringRate,
                        Neuron = rep(rep(orientations, n_neurons), each = length(orientations)),
                        Color = rep(colors, each = length(orientations)))

g3 <- ggplot(df_tuning[df_tuning$Neuron %in% seq(1, n_neurons, by = 18), ], 
             aes(x = Orientation, y = FiringRate, color = Color)) +
  geom_line() + scale_color_identity() +
  annotate(
    "text", 
    x = 90,              
    y = 135,             
    label = "Contrast = 10%", 
    vjust = 1,           
    hjust = 0.5,         
    size = 3.5             
  ) +
  theme_classic(base_size = 11) + 
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  coord_cartesian(ylim = c(0, 130)) +
  scale_y_continuous(breaks = seq(0, 120, by = 30)) +
  guides(color = F) + ylab("Spikes")
g3

input_45 <- rnorm(n_trials, 45, 0) # zero stimulus encoding noise
responses_45 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_45[i]) 
  mu <- tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 45, trial = i)
}

input_135 <- rnorm(n_trials, 135, 0)
responses_135 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_135[i]) 
  mu <- tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 135, trial = i)
}

df_responses_2 <- as.data.frame(rbind(responses_45, responses_135))
colnames(df_responses_2) <- c("Spikes", "Neuron", "Stimulus", "Trial")

# stopCluster(cl)

### Visualization ###
# within-trial spike distributions
df_responses_high$Contrast <- "high"
df_responses_low$Contrast  <- "low"
df_responses <- rbind(df_responses_high, df_responses_low)

gaussian <- function(x, A, mu, sigma){
  A * exp(-(x - mu)^2 / (2 * sigma^2))
}

fit_gaussian <- function(data){
  
  nlsLM(
    Spikes ~ gaussian(Neuron, A, mu, sigma),
    data = data,
    start = list(
      A = max(data$Spikes),
      mu = data$Neuron[which.max(data$Spikes)],
      sigma = sd(data$Neuron)
    )
  )
}

(fit_high <- fit_gaussian(df_responses_high))
(fit_low  <- fit_gaussian(df_responses_low))
(curvature_high <- summary(fit_high)$parameters[1]/summary(fit_high)$parameters[3]^2) # peak curvature
(curvature_low  <- summary(fit_low)$parameters[1]/summary(fit_low)$parameters[3]^2)

x_seq <- seq(min(df_responses$Neuron), max(df_responses$Neuron), length.out = 200)
pred_high <- data.frame(
  Neuron = x_seq,
  Spikes = predict(fit_high, newdata = data.frame(Neuron = x_seq))
)
pred_low <- data.frame(
  Neuron = x_seq,
  Spikes = predict(fit_low, newdata = data.frame(Neuron = x_seq))
)

g4 <- ggplot(df_responses_high, aes(x = Neuron, y = Spikes)) +
  geom_point(size = 0.7, color = "grey") +
  geom_line(data = pred_high, aes(x = Neuron, y = Spikes),
            color = "#E41A1C", linewidth = 1) +
  annotate(
    "text", 
    x = 90,              
    y = 135,             
    label = paste0("Peak curvature: ", round(curvature_high, 2)), 
    vjust = 1,           
    hjust = 0.5,         
    size = 3.5             
  ) +
  theme_classic(base_size = 11) + 
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  coord_cartesian(ylim = c(0, 130)) +
  scale_y_continuous(breaks = seq(0, 120, by = 30))

g5 <- ggplot(df_responses_low, aes(x = Neuron, y = Spikes)) +
  geom_point(size = 0.7, color = "grey") +
  geom_line(data = pred_low, aes(x = Neuron, y = Spikes),
            color = "#1f78b4", linewidth = 1) +
  annotate(
    "text", 
    x = 90,              
    y = 135,             
    label = paste0("Peak curvature: ", round(curvature_low, 2)), 
    vjust = 1,           
    hjust = 0.5,         
    size = 3.5             
  ) +
  theme_classic(base_size = 11) + 
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  coord_cartesian(ylim = c(0, 130)) +
  scale_y_continuous(breaks = seq(0, 120, by = 30))

# across-trial spike distributions
# 2d density
df_responses_1$Contrast <- "high"
df_responses_2$Contrast <- "low"
df <- rbind(df_responses_1, df_responses_2)

g6 <- df %>%
  dplyr::filter(Neuron %in% c(45, 135)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `45`, y = `135`,
             color = Contrast)) +
  stat_ellipse(aes(group = interaction(Contrast, Stimulus)), level = c(0.5)) + # parametric density
  stat_ellipse(aes(group = interaction(Contrast, Stimulus)), level = c(0.7)) + 
  stat_ellipse(aes(group = interaction(Contrast, Stimulus)), level = c(0.9)) + 
  labs(x = "Neuron 45", y = "Neuron 135") +
  coord_cartesian(xlim = c(0, 150), ylim = c(0, 150)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = c(0.3, 1),
    legend.justification = c(1, 1)
  )

g7 <- df %>%
  filter(Neuron %in% c(45, 135)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `45`, y = `135`, color = Contrast)) +
  geom_density_2d(aes(linetype = factor(Stimulus)), linewidth = 0.7) + # nonparametric density
  labs(x = "Neuron 45", y = "Neuron 135",
       color = "Contrast", linetype = "Stimulus") +
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
g8 <- df %>%
  group_by(Stimulus, Contrast, Trial) %>%
  summarise(Sum_spikes = sum(Spikes), .groups = "drop") %>%
  ggplot(aes(x = Sum_spikes, color = Contrast, fill = Contrast)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  coord_cartesian(xlim = c(1000, 4000)) +
  scale_x_continuous(breaks = c(1000, 2000, 3000, 4000)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic(base_size = 11) +
  labs(
    x = "Total spikes",
    y = "Density",
    color = "Contrast",
    fill = "Contrast") +
  theme(
    legend.position = c(0.3, 1),
    legend.justification = c(1, 1)
  )

# 3d density
u <- seq(0, 2*pi, length.out = 50)
v <- seq(0, pi, length.out = 50)
x_sphere <- outer(cos(u), sin(v))
y_sphere <- outer(sin(u), sin(v))
z_sphere <- outer(rep(1,length(u)), cos(v)) 
unit_sphere <- rbind(as.vector(x_sphere), as.vector(y_sphere), as.vector(z_sphere))

df %>% 
  filter(Neuron %in% c(45, 90, 135)) %>% 
  pivot_wider(id_cols = c(Contrast, Stimulus, Trial), 
              names_from = Neuron, values_from = Spikes) %>%
  group_by(Stimulus, Contrast) %>% 
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
      Stimulus = .y$Stimulus,
      Contrast = .y$Contrast
    )
  }) -> ellipsoids

p1 <- plot_ly()

for (e in ellipsoids) {
  
  col <- ifelse(e$Contrast == 25, "Reds", "Blues")
  
  p1 <- p1 %>% add_surface(
    x = e$x,
    y = e$y,
    z = e$z,
    colorscale = col,
    showscale = FALSE,
    opacity = 0.6,
    name = paste("Stim:", e$Stimulus, "Con:", e$Contrast),
    hoverinfo = "text",
    text = paste("Stimulus:", e$Stimulus,
                 "<br>Contrast:", e$Contrast)
  )
}

p1 <- p1 %>%
  layout(
    scene = list(
      xaxis = list(title = "Neuron 45", range = c(0,100)),
      yaxis = list(title = "Neuron 90", range = c(0,100)),
      zaxis = list(title = "Neuron 135", range = c(0,100))
    )
  )
p1

# 3d scatter plot
df %>% 
  filter(Neuron %in% c(45, 90, 135)) %>% 
  pivot_wider(id_cols = c(Contrast, Stimulus, Trial), 
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
    color = ~factor(Contrast),
    colors = colors,   # ←ここで指定
    name = paste("Stimulus", stim_list[i])
  )
}

p2 <- p2 %>%
  layout(
    scene = list(
      xaxis = list(title = "Neuron 45", range = c(0,100)),
      yaxis = list(title = "Neuron 90", range = c(0,100)),
      zaxis = list(title = "Neuron 135", range = c(0,100))
    )
  )

p2 <- p2 %>%
  layout(
    scene = list(
      xaxis = list(
        title = "Neuron 45",
        range = c(0,100),
        color = "black", 
        tickcolor = "black",
        titlefont = list(color = "black")
      ),
      yaxis = list(
        title = "Neuron 90",
        range = c(0,100),
        color = "black",
        tickcolor = "black",
        titlefont = list(color = "black")
      ),
      zaxis = list(
        title = "Neuron 135",
        range = c(0,100),
        color = "black",
        tickcolor = "black",
        titlefont = list(color = "black")
      )
    )
  )

p2

# 3d scatter plot with logistic regression
contrast_list <- unique(df_scatter$Contrast)
planes <- list()

for (c in contrast_list) {
  d <- df_scatter %>% 
    filter(Contrast == c)
  d$stim_bin <- as.numeric(as.factor(d$Stimulus)) - 1
  fit <- glm(stim_bin ~ `45` + `90` + `135`,
             data = d,
             family = binomial)
  b <- coef(fit)
  x_seq <- seq(0,100,length.out = 30)
  y_seq <- seq(0,100,length.out = 30)
  grid <- expand.grid(x = x_seq, y = y_seq)
  z <- -(b[1] + b[2]*grid$x + b[3]*grid$y) / b[4]
  z <- matrix(z, nrow = length(x_seq), ncol = length(y_seq))
  planes[[as.character(c)]] <- list(
    x = x_seq,
    y = y_seq,
    z = z
  )
}

for (c in contrast_list) {
  p3 <- p2 %>% add_surface(
    x = planes[[as.character(c)]]$x,
    y = planes[[as.character(c)]]$y,
    z = planes[[as.character(c)]]$z,
    opacity = 0.3,
    showscale = FALSE,
    name = paste("Decision plane Contrast", c)
  )
}

p3

# save figures
plot_list_1 <- list(g1, g3, g2)
plots_no_legend_1 <- lapply(plot_list_1, function(p) {
  p + theme(
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )
})
ps_1 <- wrap_plots(plots_no_legend_1, ncol = 3)
ggsave("ps_1.png", ps_1, width = 6, height = 2, dpi = 300)

plot_list_2 <- list(g5, g4, g8)
plots_no_legend_2 <- lapply(plot_list_2, function(p) {
  p + theme(
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )
})
ps_2 <- wrap_plots(plots_no_legend_2, ncol = 3)
ggsave("ps_2.png", ps_2, width = 6, height = 2, dpi = 300)

htmlwidgets::saveWidget(p2, "scatter_3d.html")
htmlwidgets::saveWidget(p3, "scatter_3d_logistic.html")