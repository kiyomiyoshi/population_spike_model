library(MASS)
library(mvtnorm)

n <- 10000
mu0 <- c(1, 0); Sigma0 <- matrix(c(2.5, 0.9, 0.9, 1), 2, 2)
mu1 <- c(0, 1); Sigma1 <- matrix(c(1,   0.9, 0.9, 2.5), 2, 2)

x0 <- mvrnorm(n, mu0, Sigma0)
x1 <- mvrnorm(n, mu1, Sigma1)

X <- rbind(x0, x1)
y <- c(rep(0, n), rep(1, n))
df <- data.frame(x1 = X[, 1], x2 = X[, 2], y = factor(y))

# type-1 decision
model1 <- glm(y ~ x1 + x2, data = df, family = binomial)
x_seq <- seq(min(df$x1) - 1, max(df$x1) + 1, length.out = 200)
y_seq <- seq(min(df$x2) - 1, max(df$x2) + 1, length.out = 200)
grid <- expand.grid(x1 = x_seq, x2 = y_seq)
grid_df <- expand.grid(x1 = x_seq, x2 = y_seq) %>%
  mutate(prob = predict(model1, newdata = ., type = "response"))

g1 <- ggplot() +
  geom_raster(data = grid_df, aes(x = x1, y = x2, fill = prob), interpolate = TRUE) +
  scale_fill_gradientn(
    colours = c("#c6dbef", "#f7f7f7", "#fcbba1"),
    limits = c(0, 1),
    name = "P(y=1)"
  ) +
  geom_point(data = df, aes(x = x1, y = x2, color = y), size = 1) +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  
  geom_contour(data = grid_df, aes(x = x1, y = x2, z = prob),
               breaks = 0.5, color = "black", size = 1) +
  
  coord_fixed() +
  labs(x = "Neuron 1", y = "Neuron 2", title = "Type-1 logistic regression") +
  theme_minimal() +
  theme(legend.position = "right")
g1


# type-2 decision
df$prob1 <- predict(model1, type = "response")
df_sub <- df[df$prob1 > 0.5, ]
model2 <- glm(y ~ x1 + x2, data = df_sub, family = binomial)
x_seq <- seq(min(df$x1) - 1, max(df$x1) + 1, length.out = 200)
y_seq <- seq(min(df$x2) - 1, max(df$x2) + 1, length.out = 200)
grid <- expand.grid(x1 = x_seq, x2 = y_seq)
grid$prob2 <- predict(model2, newdata = grid, type = "response")
prob_matrix <- matrix(grid$prob2, nrow = length(x_seq))
prob_masked <- prob_matrix

for (i in 1:length(x_seq)) {
  for (j in 1:length(y_seq)) {
    if (y_seq[j] < x_seq[i]) {
      prob_masked[i, j] <- NA
    }
  }
}

df_sub_plot <- df_sub %>% filter(x2 >= x1)
grid_df <- expand.grid(x1 = x_seq, x2 = y_seq) %>%
  mutate(prob = as.vector(prob_masked)) %>%
  filter(x2 >= x1)

g2 <- ggplot() +
  geom_raster(data = grid_df, aes(x = x1, y = x2, fill = prob), interpolate = TRUE) +
  scale_fill_gradientn(
    colours = c("#c6dbef","#f7f7f7","#fcbba1"),
    limits = c(0, 1),
    name = "P(y=1)"
  ) +
  geom_point(data = df_sub_plot, aes(x = x1, y = x2, color = y), size = 1) +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_contour(data = expand.grid(x1 = x_seq, x2 = y_seq) %>% mutate(prob = as.vector(prob_matrix)),
               aes(x = x1, y = x2, z = prob),
               breaks = 0.5,
               color = "black", size = 1) +
  
  coord_fixed() +
  labs(x = "Neuron 1", y = "Neuron 2", title = "Type-2 logistic regression") +
  theme_minimal() +
  theme(legend.position = "right")
g2

# bayesian decision
prior0 <- 0.5
prior1 <- 0.5

grid$px_C0 <- dmvnorm(as.matrix(grid[, c("x1","x2")]), mean = mu0, sigma = Sigma0)
grid$px_C1 <- dmvnorm(as.matrix(grid[, c("x1","x2")]), mean = mu1, sigma = Sigma1)
grid$prob <- (grid$px_C1 * prior1) / (grid$px_C0*prior0 + grid$px_C1*prior1)

g3 <- ggplot() +
  geom_raster(data = grid, aes(x = x1, y = x2, fill = prob), interpolate = TRUE) +
  scale_fill_gradient2(low = "#c6dbef", mid = "#f7f7f7", high = "#fcbba1", 
                       midpoint = 0.5,
                       limits = c(0, 1),
                       name = "P(y=1)") +
  geom_point(data = df, aes(x = x1, y = x2, color = y), size = 1) +
  # geom_contour(data = grid, aes(z = prob), breaks = 0.5, color = "black", size = 1) +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  labs(x = "Neuron 1", y = "Neuron 2", title = "Bayesian probability") +
  theme_minimal()
g3


plot_list_1 <- list(g1, g2, g3)
plots_1 <- lapply(plot_list_1, function(p) {
  p + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"))
})
logistic_correlated <- patchwork::wrap_plots(plots_1, nrow = 1)
ggsave("logistic_correlated.png", logistic_correlated, width = 12, height = 4, dpi = 300)