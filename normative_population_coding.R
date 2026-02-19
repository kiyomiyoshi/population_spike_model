library(tidyverse)

# Parameters
n_neurons <- 36
orientations <- seq(1, 360, by = 1)
preferred_orientations <-  seq(1, 360, length.out = n_neurons)
max_firing_rate <- 60
tuning_width <- 40

# Tuning curves
color_wheel <- function(alpha = 1) {
  h <- c(seq(320, 0, -1), seq(360, 321, -1))[-321]
  h <- h / 360
  s <- rep(1, 360)
  v <- rep(1, 360)
  hsv(h, s, v, alpha)
}

colors <- color_wheel()

tuning_curves <- matrix(0, nrow = n_neurons, ncol = length(orientations))
for (i in 1:n_neurons) {
  tuning_curves[i, ] <- max_firing_rate * exp(-0.5 * (pmin(abs(orientations - preferred_orientations[i]), 
                                                           360 - abs(orientations - preferred_orientations[i])) ^ 2) / tuning_width ^ 2)
}

FiringRate <- c()
for (i in 1:n_neurons) {
  FiringRate <- c(FiringRate, tuning_curves[i, ])
}

df_tuning <- data.frame(Orientation = rep(orientations, n_neurons),
                        FiringRate = FiringRate,
                        Neuron = rep(rep(orientations, n_neurons), each = length(orientations)),
                        Color = rep(colors, each = length(orientations)))

p1 <- ggplot(df_tuning[df_tuning$Neuron %in% seq(1, n_neurons, by = 6), ], 
             aes(x = Orientation, y = FiringRate, color = Color)) +
  geom_line() + scale_color_identity() +
  theme_minimal() + 
  ylim(0, 65) + guides(color = F) + ylab("Firing rate")

# Spikes
input_orientation <- 120

responses <- rpois(n_neurons, lambda = tuning_curves[, input_orientation + 1])
df_responses <- data.frame(Spikes = responses, Neuron = seq(1:n_neurons))

responses_augmented <- rpois(n_neurons, lambda = tuning_curves[, input_orientation + 1]) + rpois(n_neurons, lambda = 10)
df_responses_augmented <- data.frame(Spikes = responses_augmented, Neuron = seq(1:n_neurons))

p2 <- ggplot(df_responses) + geom_col(aes(x = Neuron, y = Spikes), color = "grey34") +
  theme_minimal() + ylim(0, 80) + ylab("Spikes")

p3 <- ggplot(df_responses_augmented) + geom_col(aes(x = Neuron, y = Spikes), color = "grey34") +
  theme_minimal() + ylim(0, 80) + ylab("Spikes")

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
prior_mean <- 120 
prior_sd <- 60    
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
  labs(x = "Orientation (degrees)", y = "Posterior Probability") +
  theme_minimal() + xlim(100, 140) + ylim(0, 0.35)

p5 <- ggplot(df_posterior_augmented, aes(x = Orientation, y = Density)) +
  geom_line() +
  labs(x = "Orientation (degrees)", y = "Posterior Probability") +
  theme_minimal() + xlim(100, 140) + ylim(0, 0.35)

g1 <- cowplot::plot_grid(p2, p3)
g2 <- cowplot::plot_grid(p4, p5)

g1
g2