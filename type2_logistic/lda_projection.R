# 2つのニューロンが「同じ方向に揺れる」ため、総和方向（信号方向）での変動が増える
# ニューロンが独立なら、一方が大きく揺れても、もう一方が逆方向に揺れることもあるので、総和方向では揺れがある程度打ち消される
# これによって総和方向の分散は抑えられる

library(ggplot2)
library(MASS)

n <- 1000

mu_0 <- c(30, 29) # 非負じゃないとgainをかけたときおかしくなる
mu_1 <- c(29, 30)

Sigma_0 <- matrix(c(1.0, 0, 0, 1.0), nrow = 2)
Sigma_1 <- matrix(c(1.0, 0, 0, 1.0), nrow = 2)

data0 <- mvrnorm(n, mu = mu_0, Sigma = Sigma_0) # class 0
data1 <- mvrnorm(n, mu = mu_1, Sigma = Sigma_1) # class 1

data <- data.frame(
  class = factor(c(rep(0, n), rep(1, n))),
  neuron1 = c(data0[, 1], data1[, 1]),
  neuron2 = c(data0[, 2], data1[, 2])
)

neurons <- c("neuron1", "neuron2")

data0_sub <- data[data$class == 0, neurons]
data1_sub <- data[data$class == 1, neurons]
df <- colMeans(data1_sub) - colMeans(data0_sub)

# 共分散行列と正則化
lambda <- 0.01
Sigma <- (cov(data0) + cov(data1)) / 2
Sigma <- Sigma + lambda * diag(length(neurons))

# 判別軸 w
w <- solve(Sigma) %*% df
w_scaled <- as.numeric(w)

# 平均ベクトル
mean0 <- colMeans(data0)
mean1 <- colMeans(data1)

(lfi <- t(df) %*% solve(Sigma) %*% df)
(dp <- sqrt(lfi))

# 2d plot
arrows_df <- data.frame(
  x = c(mean0[1], mean0[1]),
  y = c(mean0[2], mean0[2]),
  xend = c(mean0[1] + df[1], mean0[1] + w_scaled[1]),
  yend = c(mean0[2] + df[2], mean0[2] + w_scaled[2]),
  type = c("df", "w")
)

g1 <- ggplot(data, aes(x = neuron1, y = neuron2, color = class)) +
  geom_point(size = 2) +
  geom_point(aes(x = mean0[1], y = mean0[2]), color = "blue", size = 2) +
  geom_point(aes(x = mean1[1], y = mean1[2]), color = "red", size = 2) +
  geom_segment(data = arrows_df, 
               aes(x = x, y = y, xend = xend, yend = yend, color = type),
               arrow = arrow(length = unit(0.15,"cm")), lwd = 1.2) +
  coord_fixed() +
  scale_color_manual(values = c("0" = "blue", "1" = "red", "df" = "green", "w" = "purple")) +
  labs(title = "平均値差ベクトル df & 判別軸 w",
       subtitle = paste("d' =", round(dp, 3)),
       x = "Neuron 1", y = "Neuron 2", color = "Legend") +
  xlim(25, 35) + ylim(25, 35) +
  theme_minimal()
g1

# 1d plot
proj0 <- as.numeric(data0 %*% w)
proj1 <- as.numeric(data1 %*% w)

proj_df <- data.frame(
  projection = c(proj0, proj1),
  class = factor(c(rep(0, n), rep(1, n)))
)

mu_proj0 <- mean(proj0)
mu_proj1 <- mean(proj1)
sd_proj <- sqrt((var(proj0) + var(proj1)) / 2) # mean squared sd

(dp_1d <- (mu_proj1 - mu_proj0) / sd_proj)

g2 <- ggplot(proj_df, aes(x = projection, fill = class)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = mu_proj0, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = mu_proj1, color = "red", linetype = "dashed") +
  scale_fill_manual(values = c("0" = "blue", "1" = "red")) +
  labs(title = "判別軸 w への射影後の1次元分布",
       subtitle = paste("d' =", round(dp_1d, 3)),
       x = "Projection onto w",
       y = "Density") +
  xlim(-6, 6) +
  theme_minimal()
g2


### gain variability
sigma_g <- 0.05
gain0 <- rgamma(n, shape = 1 / sigma_g^2, scale = sigma_g^2)
gain1 <- rgamma(n, shape = 1 / sigma_g^2, scale = sigma_g^2)

data0[, 1] <- data0[, 1] * gain0
data0[, 2] <- data0[, 2] * gain0

data1[, 1] <- data1[, 1] * gain1
data1[, 2] <- data1[, 2] * gain1

data <- data.frame(
  class = factor(c(rep(0, n), rep(1, n))),
  neuron1 = c(data0[, 1], data1[, 1]),
  neuron2 = c(data0[, 2], data1[, 2])
)

data0_sub <- data[data$class == 0, neurons] # なぜ必要？　subとは？
data1_sub <- data[data$class == 1, neurons]
df <- colMeans(data1_sub) - colMeans(data0_sub)

# 共分散行列と正則化
lambda <- 0.01
Sigma <- (cov(data0) + cov(data1)) / 2
Sigma <- Sigma + lambda * diag(length(neurons))

# 判別軸 w
w <- solve(Sigma) %*% df
w_scaled <- as.numeric(w)

# 平均ベクトル
mean0 <- colMeans(data0)
mean1 <- colMeans(data1)

(lfi <- t(df) %*% solve(Sigma) %*% df)
(dp <- sqrt(lfi))

# 2d plot
arrows_df <- data.frame(
  x = c(mean0[1], mean0[1]),
  y = c(mean0[2], mean0[2]),
  xend = c(mean0[1] + df[1], mean0[1] + w_scaled[1]),
  yend = c(mean0[2] + df[2], mean0[2] + w_scaled[2]),
  type = c("df", "w")
)

g3 <- ggplot(data, aes(x = neuron1, y = neuron2, color = class)) +
  geom_point(size = 2) +
  geom_point(aes(x = mean0[1], y = mean0[2]), color = "blue", size = 2) +
  geom_point(aes(x = mean1[1], y = mean1[2]), color = "red", size = 2) +
  geom_segment(data = arrows_df, 
               aes(x = x, y = y, xend = xend, yend = yend, color = type),
               arrow = arrow(length = unit(0.15,"cm")), lwd = 1.2) +
  coord_fixed() +
  scale_color_manual(values = c("0" = "blue", "1" = "red", "df" = "green", "w" = "purple")) +
  labs(title = "平均値差ベクトル df & 判別軸 w",
       subtitle = paste("d' =", round(dp, 3)),
       x = "Neuron 1", y = "Neuron 2", color = "Legend") +
  xlim(25, 35) + ylim(25, 35) +
  theme_minimal()
g3

# 1d plot
proj0 <- as.numeric(data0 %*% w)
proj1 <- as.numeric(data1 %*% w)

proj_df <- data.frame(
  projection = c(proj0, proj1),
  class = factor(c(rep(0, n), rep(1, n)))
)

mu_proj0 <- mean(proj0)
mu_proj1 <- mean(proj1)
sd_proj <- sqrt((var(proj0) + var(proj1)) / 2) # mean squared sd

(dp_1d <- (mu_proj1 - mu_proj0) / sd_proj)

g4 <- ggplot(proj_df, aes(x = projection, fill = class)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = mu_proj0, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = mu_proj1, color = "red", linetype = "dashed") +
  scale_fill_manual(values = c("0" = "blue", "1" = "red")) +
  labs(title = "判別軸 w への射影後の1次元分布",
       subtitle = paste("d' =", round(dp_1d, 3)),
       x = "Projection onto w",
       y = "Density") +
  xlim(-6, 6) +
  theme_minimal()
g4

plot_list_1 <- list(g1, g3, g2, g4)
plots_1 <- lapply(plot_list_1, function(p) {
  p + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"))
})
lda_projection <- patchwork::wrap_plots(plots_1, nrow = 2)
ggsave("lda_projection.png", lda_projection, width = 7, height = 6, dpi = 300)