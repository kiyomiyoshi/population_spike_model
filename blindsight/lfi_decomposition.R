# 2dプロットに最適判別面の線を引くことで、それらニューロンの寄与を可視化できるかも

library(tidyverse)

compute_lfi_decomposition <- function(data, lambda = 1e-6) {
  
  neurons <- grep("^[0-9]+$", colnames(data), value = TRUE)
  p <- length(neurons)
  
  data0 <- data[data$Stimulus_bin == 0, neurons]
  data1 <- data[data$Stimulus_bin == 1, neurons]
  
  df <- colMeans(data1) - colMeans(data0)
  
  Sigma <- (cov(data0) + cov(data1)) / 2
  Sigma <- Sigma + lambda * diag(p) # 正則化
  
  # 1. 最適なデコーディング重み w = Sigma^-1 * df
  # (これはロジスティック回帰の重みベクトルと数学的に等価な方向を指します)
  w <- solve(Sigma) %*% df
  
  # 2. 各ニューロンの貢献度 (Element-wise contribution)
  # I = t(df) %*% w = sum(df * w) なので、各要素がその細胞の寄与分
  contribution <- as.numeric(df * w)
  
  # 3. 独立仮定時の情報量
  # 対角成分のみを取り出し、相関を0にする
  Sigma_diag <- diag(diag(Sigma))
  w_shuffled <- solve(Sigma_diag) %*% df
  I_shuffled <- as.numeric(t(df) %*% w_shuffled)
  
  return(list(
    total_I = sum(contribution),
    neuron_contributions = contribution,
    shuffled_I = I_shuffled,
    correlation_gain = sum(contribution) / I_shuffled
  ))
}

maxFR_values <- 10

df_sub_high <- df_high_wide[df_high_wide$Contrast == maxFR, ]
lfi_high <- compute_lfi_decomposition(df_sub_high)
lfi_high
plot(cumsum(sort(lfi_high$neuron_contributions, decreasing = TRUE)))

plot_data_high <- data.frame(
  Neuron_ID = 1:length(lfi_high$neuron_contributions),
  Contribution = lfi_high$neuron_contributions
)

gg1 <- ggplot(plot_data_high, aes(x = Neuron_ID, y = Contribution)) +
  geom_bar(stat = "identity", fill = "#E41A1C") +
  theme_classic() +
  labs(
    title = "High Gain Variability",
    subtitle = paste("Total LFI:", round(lfi_high$total_I, 2), 
                     " | Correlation Gain:", round(lfi_high$correlation_gain, 2)),
    x = "Neuron ID (1 - 180)",
    y = "Contribution (df * w)"
  ) +
  scale_x_continuous(breaks = seq(0, 180, by = 30)) +
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) 

df_sub_low <- df_low_wide[df_low_wide$Contrast == maxFR, ]
lfi_low <- compute_lfi_decomposition(df_sub_low)
lfi_low
plot(cumsum(sort(lfi_low$neuron_contributions, decreasing = TRUE)))

plot_data_low <- data.frame(
  Neuron_ID = 1:length(lfi_low$neuron_contributions),
  Contribution = lfi_low$neuron_contributions
)
gg1

gg2 <- ggplot(plot_data_low, aes(x = Neuron_ID, y = Contribution)) +
  geom_bar(stat = "identity", fill = "#377EB8") +
  theme_classic() +
  labs(
    title = "Low Gain Variability",
    subtitle = paste("Total LFI:", round(lfi_low$total_I, 2), 
                     " | Correlation Gain:", round(lfi_low$correlation_gain, 2)),
    x = "Neuron ID (1 - 180)",
    y = "Contribution (df * w)"
  ) +
  scale_x_continuous(breaks = seq(0, 180, by = 30)) +
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) 
gg2

combined_data <- data.frame(
  Neuron_ID = 1:length(lfi_high$neuron_contributions),
  High = lfi_high$neuron_contributions,
  Low = lfi_low$neuron_contributions
) %>%
  pivot_longer(
    cols = c(High, Low),
    names_to = "Condition",
    values_to = "Contribution"
  )

gg3 <- ggplot(combined_data, aes(x = Neuron_ID, y = Contribution, color = Condition)) +
  # 線画（geom_line）で重ね合わせる
  geom_line(size = 1, alpha = 0.8) +
  # カラーパレットの指定（ご提示の色を反映）
  scale_color_manual(values = c("High" = "#E41A1C", "Low" = "#377EB8")) +
  # テーマ設定
  theme_classic() +
  labs(
    title = "High vs Low Gain Variability",
    x = "Neuron ID (1 - 180)",
    y = "Contribution (df * w)",
    color = "Condition"
  ) +
  # 軸の設定（ご提示の設定を反映）
  scale_x_continuous(breaks = seq(0, 180, by = 30)) +
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) +
  # 凡例の位置調整
  theme(legend.position = "top")


plot_list <- list(gg1, gg2, gg3)
plots <- lapply(plot_list, function(p) {
  p + theme(
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )
})

lfi <- wrap_plots(plots, row = 1)
ggsave("lfi.png", lfi, width = 10, height = 3, dpi = 300)