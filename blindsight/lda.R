library(ggplot2)

set.seed(123)
n <- 1000
neuron1_0 <- rnorm(n, mean = 0, sd = 1)
neuron2_0 <- rnorm(n, mean = 0, sd = 2)
neuron1_1 <- rnorm(n, mean = 2, sd = 1)
neuron2_1 <- rnorm(n, mean = 1, sd = 2)

data <- data.frame(
  Stimulus_bin = factor(c(rep(0, n), rep(1, n))),
  neuron1 = c(neuron1_0, neuron1_1),
  neuron2 = c(neuron2_0, neuron2_1)
)

neurons <- c("neuron1", "neuron2")

# 平均差ベクトル df
data0 <- data[data$Stimulus_bin == 0, neurons]
data1 <- data[data$Stimulus_bin == 1, neurons]
df <- colMeans(data1) - colMeans(data0)

# 共分散行列と正則化
lambda <- 0.01
Sigma <- (cov(data0) + cov(data1)) / 2
Sigma <- Sigma + lambda * diag(length(neurons))

# 判別軸 w
w <- solve(Sigma) %*% df
w_scaled <- as.numeric(w) * 2   # 見やすいようにスケーリング

# 平均ベクトル
mean0 <- colMeans(data0)
mean1 <- colMeans(data1)

# ggplot用データフレーム
arrows_df <- data.frame(
  x = c(mean0[1], mean0[1]),
  y = c(mean0[2], mean0[2]),
  xend = c(mean0[1] + df[1], mean0[1] + w_scaled[1]),
  yend = c(mean0[2] + df[2], mean0[2] + w_scaled[2]),
  type = c("df", "w")
)

ggplot(data, aes(x = neuron1, y = neuron2, color = Stimulus_bin)) +
  geom_point(size = 2) +
  geom_point(aes(x = mean0[1], y = mean0[2]), color = "blue", size = 4) +
  geom_point(aes(x = mean1[1], y = mean1[2]), color = "red", size = 4) +
  geom_segment(data = arrows_df, 
               aes(x = x, y = y, xend = xend, yend = yend, color = type),
               arrow = arrow(length = unit(0.15,"cm")), lwd = 1.2) +
  coord_fixed() +
  scale_color_manual(values = c("0" = "blue", "1" = "red", "df" = "green", "w" = "purple")) +
  labs(title = "判別軸 w と平均差ベクトル df の可視化",
       x = "Neuron 1", y = "Neuron 2", color = "Legend") +
  theme_minimal()