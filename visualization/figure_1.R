library(tidyverse)
library(cowplot)
library(patchwork)
library(ggpp)

set.seed(123)

color_wheel <- function(alpha = 1) {
  h <- c(seq(160, 0, -1), seq(180, 161, -1))[-161]
  h <- h / 180
  s <- rep(1, 180)
  v <- rep(1, 180)
  hsv(h, s, v, alpha)
}

colors <- color_wheel()

n_neurons <- 180
orientations <- seq(1, 180, by = 1)
preferred_orientations <-  seq(1, 180, length.out = n_neurons)
max_firing_rate <- 115
tuning_width <- 20

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

p1 <- ggplot(df_tuning[df_tuning$Neuron %in% seq(1, n_neurons, by = 18), ], 
             aes(x = Orientation, y = FiringRate, color = Color)) +
  geom_line() + scale_color_identity() +
  annotate("text_npc",
           npcx = 0.5,
           npcy = 1,
           label = "Tuning curves", 
           size = 3.5) +
  theme_classic(base_size = 12) + 
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  ylim(0, 160) + 
  guides(color = F) + ylab("Spikes")
p1

# Spikes
input_orientation <- 90

responses <- rpois(n_neurons, lambda = tuning_curves[, input_orientation])
df_responses <- data.frame(Spikes = responses, Neuron = seq(1:n_neurons))

responses_augmented <- rpois(n_neurons, lambda = tuning_curves[, input_orientation]) + rpois(n_neurons, lambda = 20)
df_responses_augmented <- data.frame(Spikes = responses_augmented, Neuron = seq(1:n_neurons))

p2 <- ggplot(df_responses) + geom_col(aes(x = Neuron, y = Spikes), color = "grey34") +
  theme_classic(base_size = 12) + 
  annotate("text_npc",
           npcx = 0.5,
           npcy = 1,
           label = "Baseline excitability", 
           size = 3.5) +
  ylim(0, 165) + ylab("Spikes")

p3 <- ggplot(df_responses_augmented) + geom_col(aes(x = Neuron, y = Spikes), color = "grey34") +
  theme_classic(base_size = 12) + 
  annotate("text_npc",
           npcx = 0.5,
           npcy = 1,
           label = "Elevated excitability", 
           size = 3.5) +
  ylim(0, 165) + ylab("Spikes")

# Likelihood
likelihood <- numeric(length(orientations))
for (i in 1:length(orientations)) {                                               # 尤度（i度の刺激下で観測発火パターンが生じる確率）
  likelihood[i] <- sum(dpois(responses, lambda = tuning_curves[, i], log = TRUE)) # i度刺激下での360個のニューロンのλがチューニングカーブから決まる
}                                                                                 # これら360個のポアソン分布から観測発火パターンが生じる確率を求める
df_likelihood <- data.frame(Orientation = orientations, Density = likelihood, Type = "Likelihood")
max_likelihood_orientation <- orientations[which.max(likelihood)]
cat("The maximum likelihood orientation is", max_likelihood_orientation, "degrees\n")

likelihood_augmented <- numeric(length(orientations))
for (i in 1:length(orientations)) {                               
  likelihood_augmented[i] <- sum(dpois(responses_augmented, lambda = tuning_curves[, i], log = TRUE))
} 
df_likelihood_augmented <- data.frame(Orientation = orientations, Density = likelihood_augmented, Type = "Likelihood")
max_likelihood_orientation_augmented <- orientations[which.max(likelihood_augmented)]
cat("The maximum likelihood orientation is", max_likelihood_orientation_augmented, "degrees\n")

# Bayesian inference
prior_mean <- 90 
prior_sd <- 100    
prior <- dnorm(orientations, mean = prior_mean, sd = prior_sd)
df_prior <- data.frame(Orientation = orientations, Density = prior, Type = "Prior")

posterior_log <- likelihood + log(prior)
posterior_unnormalized <- exp(posterior_log - max(posterior_log))
posterior <- posterior_unnormalized / sum(posterior_unnormalized)
df_posterior <- data.frame(Orientation = orientations, Density = posterior, Type = "Posterior")
map_orientation <- orientations[which.max(posterior)]
cat("The maximum a posteriori (MAP) orientation is", map_orientation, "degrees\n")

posterior_log_augmented <- likelihood_augmented + log(prior)
posterior_unnormalized_augmented <- exp(posterior_log_augmented - max(posterior_log_augmented))
posterior_augmented <- posterior_unnormalized_augmented / sum(posterior_unnormalized_augmented)
df_posterior_augmented <- data.frame(Orientation = orientations, Density = posterior_augmented, Type = "Posterior")
map_orientation_augmented <- orientations[which.max(posterior_augmented)]
cat("The maximum a posteriori (MAP) orientation is", map_orientation_augmented, "degrees\n")

p4 <- ggplot(df_posterior, aes(x = Orientation, y = Density)) +
  geom_line() +
  labs(x = "Orientation", y = "Density") +
  annotate("text_npc",
           npcx = 0.5,
           npcy = 1,
           label = "Bayesian inference", 
           size = 3.5) +
  theme_classic(base_size = 12) + 
  coord_cartesian(xlim = c(85, 95)) +
  scale_x_continuous(breaks = c(85, 90, 95)) + 
  scale_y_continuous(limits = c(0, 1.13), breaks = c(0, 0.5, 1))


p5 <- ggplot(df_posterior_augmented, aes(x = Orientation, y = Density)) +
  geom_line() +
  theme_classic(base_size = 12) + 
  labs(x = "Orientation", y = "Density") +
  annotate("text_npc",
           npcx = 0.5,
           npcy = 1,
           label = "Posterior unchanged", 
           size = 3.5) +
  coord_cartesian(xlim = c(85, 95)) +
  scale_x_continuous(breaks = c(85, 90, 95)) + 
  scale_y_continuous(limits = c(0, 1.13), breaks = c(0, 0.5, 1))


### sdt ###
x_vals <- seq(-6, 6, length.out = 1000)
df <- tibble(
  x = x_vals,
  density_1 = dnorm(x_vals, mean = -0.75, sd = 1.5),
  density_2 = dnorm(x_vals, mean =  0.75, sd = 1.5)
) %>%
  pivot_longer(cols = starts_with("density"),
               names_to = "distribution", values_to = "density")

labels <- c("―CW―", "―CCW―") 
x_pos  <- c(-0.60, 1.38)

p6 <- ggplot(df, aes(x = x, y = density, color = distribution)) +
  geom_segment(x =  0,   xend =  0,   y = 0, yend = 0.40, linetype = "solid",  color = "black",  size = 0.6) +
  geom_segment(x = -1.5, xend = -1.5, y = 0, yend = 0.35, linetype = "dashed", color = "grey70", size = 0.5) +
  geom_segment(x = -3.0, xend = -3.0, y = 0, yend = 0.35, linetype = "dashed", color = "grey70", size = 0.5) +
  geom_segment(x =  1.5, xend =  1.5, y = 0, yend = 0.35, linetype = "dashed", color = "grey70", size = 0.5) +
  geom_segment(x =  3.0, xend =  3.0, y = 0, yend = 0.35, linetype = "dashed", color = "grey70", size = 0.5) +
  geom_line(size = 0.7, alpha = 1) +
  annotate("text", x = -2.4,    y = 0.39, label = "CCW", size = 3.5) +
  annotate("text", x =  2.0,    y = 0.39, label = "CW", size = 3.5) +
  annotate("text", x = -0.75, y = 0.32, label = 1, size = 3.5) +
  annotate("text", x = -2.25, y = 0.32, label = 2, size = 3.5) +
  annotate("text", x = -3.75, y = 0.32, label = 3, size = 3.5) +
  annotate("text", x =  0.75, y = 0.32, label = 1, size = 3.5) +
  annotate("text", x =  2.25, y = 0.32, label = 2, size = 3.5) +
  annotate("text", x =  3.75, y = 0.32, label = 3, size = 3.5) +
  # scale_color_manual(values = c("density_1" = color_1, "density_2" = color_2)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2), limits = c(-6, 6)) +
  scale_y_continuous(limits = c(0, 0.40)) +
  labs(x = "Signal strength", y = "Density") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x  = element_blank())

f1 <- list(p1, p2, p4, p3, p5, 
           p6, p6, p6, p6, p6)
f1 <- lapply(f1, function(p) {
  p + theme(
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )
})

figure_1 <- wrap_plots(f1, nrow = 2)
ggsave("figure_1.png", figure_1, width = 10, height = 3, dpi = 300)