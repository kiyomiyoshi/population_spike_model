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

cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

# Naka-Rushton function
Rmax <- 115
C50 <- 19.3
n   <-  2.9

# High gain variability
n_neurons <- 180                                              
orientations <- seq(1, 180, by = 1)                           
preferred_orientations <-  seq(1, 180, length.out = n_neurons)
contrast <- 5
max_firing_rate <- Rmax * contrast^n / (contrast^n + C50^n)
tuning_width    <- 75                                         
n_trials <- 3000

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
  sigma_g <- 0.2
  g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
  mu <- g * tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 90, trial = i)
}

input_180 <- rnorm(n_trials, 180, 0)
responses_180 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_180[i]) 
  sigma_g <- 0.2
  g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
  mu <- g * tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 180, trial = i)
}

df_responses_1 <- as.data.frame(rbind(responses_90, responses_180))
colnames(df_responses_1) <- c("Spikes", "Neuron", "Stimulus", "Trial")

# Low gain variability
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

input_180 <- rnorm(n_trials, 180, 0)
responses_180 <- foreach(
  i = 1:n_trials,
  .combine = rbind,
  .packages = "stats"
) %dopar% {
  idx <- round(input_180[i]) 
  sigma_g <- 0.01
  g <- rgamma(1, shape = 1 / sigma_g^2, scale = sigma_g^2)
  mu <- g * tuning_curves[, idx]
  mu[mu == 0] <- 1e-8
  resp <- rnbinom(n_neurons, size = 1000000, mu = mu)
  cbind(resp, neuron = seq_len(n_neurons), stim = 180, trial = i)
}

df_responses_2 <- as.data.frame(rbind(responses_90, responses_180))
colnames(df_responses_2) <- c("Spikes", "Neuron", "Stimulus", "Trial")

# stopCluster(cl)

### Visualization ###
# across-trial spike distributions
# 2d density
df_responses_1$GV <- "high"
df_responses_2$GV <- "low"
df <- rbind(df_responses_1, df_responses_2)

g1 <- df %>%
  dplyr::filter(Neuron %in% c(90, 180)) %>%
  pivot_wider(names_from = Neuron, values_from = Spikes) %>%
  ggplot(aes(x = `90`, y = `180`,
             color = GV)) +
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.5)) + # parametric density
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.7)) + 
  stat_ellipse(aes(group = interaction(GV, Stimulus)), level = c(0.9)) + 
  labs(x = "Neuron 90", y = "Neuron 180") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
  scale_x_continuous(breaks = c(0, 5, 10)) +
  scale_y_continuous(breaks = c(0, 5, 10)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
g1

# sum of spikes
g2 <- df %>%
  group_by(Stimulus, GV, Trial) %>%
  summarise(Sum_spikes = sum(Spikes), .groups = "drop") %>%
  ggplot(aes(x = Sum_spikes, color = GV, fill = GV)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  coord_cartesian(xlim = c(0, 600)) +
  scale_x_continuous(breaks = c(0, 200, 400, 600)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic(base_size = 11) +
  labs(
    x = "Total spikes",
    y = "Density",
    color = "GV",
    fill = "GV") +
  theme(
    legend.position = c(1.1, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(color = NA),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )
g2

# Pairwise spike correlations
# 3d density
u <- seq(0, 2*pi, length.out = 50)
v <- seq(0, pi, length.out = 50)
x_sphere <- outer(cos(u), sin(v))
y_sphere <- outer(sin(u), sin(v))
z_sphere <- outer(rep(1,length(u)), cos(v)) 
unit_sphere <- rbind(as.vector(x_sphere), as.vector(y_sphere), as.vector(z_sphere))

df %>% 
  filter(Neuron %in% c(90, 135, 180)) %>% 
  pivot_wider(id_cols = c(GV, Stimulus, Trial), 
              names_from = Neuron, values_from = Spikes) %>%
  group_by(Stimulus, GV) %>% 
  group_map(~{
    mu <- colMeans(.x[,c("90","135","180")])
    sigma <- cov(.x[,c("90","135","180")])
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
  
  col <- ifelse(e$GV == "high", "#E41A1C", "#377EB8")
  
  p1 <- p1 %>% add_surface(
    x = e$x,
    y = e$y,
    z = e$z,
    colorscale = list(c(0,"white"), c(1,col)),
    showscale = FALSE,
    opacity = 0.6,
    name = paste("Stim:", e$Stimulus, "GV:", e$GV),
    hoverinfo = "text",
    text = paste("Stimulus:", e$Stimulus,
                 "<br>GV:", e$GV)
  )
}

p1 <- p1 %>%
  layout(
    scene = list(
      xaxis = list(title = "Neuron 90",  range = c(0, 10)),
      yaxis = list(title = "Neuron 135", range = c(0, 10)),
      zaxis = list(title = "Neuron 180", range = c(0, 10)),
      aspectmode = "cube"   
    )
  )
p1

# 3d scatter plot
df %>% 
  filter(Neuron %in% c(90, 135, 180)) %>% 
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
    y = ~`135`,
    z = ~`180`,
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
      xaxis = list(title = "Neuron 90",  range = c(0, 10)),
      yaxis = list(title = "Neuron 135", range = c(0, 10)),
      zaxis = list(title = "Neuron 180", range = c(0, 10)),
      aspectmode = "cube"   
    )
  )

p2 <- p2 %>%
  layout(
    scene = list(
      xaxis = list(
        title = "Neuron 90",
        range = c(0, 10),
        color = "black", 
        tickcolor = "black",
        titlefont = list(color = "black")
      ),
      yaxis = list(
        title = "Neuron 135",
        range = c(0, 10),
        color = "black",
        tickcolor = "black",
        titlefont = list(color = "black")
      ),
      zaxis = list(
        title = "Neuron 180",
        range = c(0, 10),
        color = "black",
        tickcolor = "black",
        titlefont = list(color = "black")
      )
    )
  )

p2

# 3d scatter plot with logistic regression
df_scatter$stim_bin <- as.numeric(as.factor(df_scatter$Stimulus)) - 1
fit <- glm(stim_bin ~ `90` + `135` + `180`, data = df_scatter, family = binomial)

x_seq <- seq(0, 10, length.out = 30)
y_seq <- seq(0, 10, length.out = 30)
z_seq <- seq(0, 10, length.out = 30)

grid3d <- expand.grid(`90` = x_seq, `135` = y_seq, `180` = z_seq)
grid3d$prob <- predict(fit, newdata = grid3d, type = "response")
decision_points <- grid3d %>% filter(abs(prob - 0.5) < 0.02)

p3 <- p1 %>%
  add_trace(
    data = decision_points,
    x = ~`90`, y = ~`135`, z = ~`180`,
    type = "mesh3d",
    color = 'navy',
    opacity = 0.5,
    name = "Decision boundary"
  )

p3

# save figures
plot_list <- list(g1, g2)
plots <- lapply(plot_list, function(p) {
  p + theme(
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )
})
manifold <- wrap_plots(plots, ncol = 2)
ggsave("manifold_90_180.png", manifold, width = 4, height = 1.95, dpi = 300)

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

saveWidget(p2, "scatter_3d_90_180.html", selfcontained = TRUE)
saveWidget(p3, "scatter_3d_logistic_90_180.html", selfcontained = TRUE)
webshot("scatter_3d_90_180.html", "scatter_3d_90_180.png", vwidth = 800, vheight = 600)
webshot("scatter_3d_logistic_90_180.html", "scatter_3d_logistic_90_180.png", vwidth = 800, vheight = 600)