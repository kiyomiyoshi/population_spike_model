library(tidyverse)

df_pred <- df_summary %>%
  arrange(Contrast) %>%
  mutate(
    mu = Mean_sum_spikes,
    sigma = sqrt(mu)
  )

dC <- diff(df_pred$Contrast)
dmu <- diff(df_pred$mu)

# μ'
mu_prime <- c(
  dmu[1] / dC[1],                                 　　　　　# 前方差分（最初）
  (dmu[-1]/dC[-1] + dmu[-length(dmu)]/dC[-length(dC)]) / 2, # 中央差分
  dmu[length(dmu)] / dC[length(dC)]                         # 後方差分（最後）
)

df_pred$mu_prime <- mu_prime

# 予測弁別閾（比例定数なし）
df_pred$DeltaC_pred <- df_pred$sigma / df_pred$mu_prime
df_pred

0.0462 * 221 / 10
0.0685 * 263 / 18

df_pred <- df_pred %>%
  mutate(
    Weber_ratio_pred = DeltaC_pred / Contrast) # ΔC / C

df_pred

ggplot(df_pred, aes(x = Contrast, y = DeltaC_pred)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  ylab("Predicted Delta C") +
  xlab("Contrast")

ggplot(df_pred, aes(x = Contrast, y = Weber_ratio_pred)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  ylab("Predicted Weber Ratio") +
  xlab("Contrast")

df_pred <- df_pred %>%
  mutate(
    Perceived = log(Mean_sum_spikes + 1)
  )

Perceived <- df_pred$Perceived
dPerceived <- diff(Perceived)

# コントラスト差
dC <- diff(df_pred$Contrast)

# 中央差分風に μ'perceived を計算
mu_prime_perceived <- c(
  dPerceived[1]/dC[1],
  (dPerceived[-1]/dC[-1] + dPerceived[-length(dPerceived)]/dC[-length(dC)]) / 2,
  dPerceived[length(dPerceived)]/dC[length(dC)]
)

df_pred$mu_prime_perceived <- mu_prime_perceived

df_pred$DeltaC_pred_perceived <- df_pred$sigma / df_pred$mu_prime_perceived

df_pred$Weber_ratio_perceived <- df_pred$DeltaC_pred_perceived / df_pred$Contrast

library(ggplot2)

ggplot(df_pred, aes(x = Contrast)) +
  geom_line(aes(y = Weber_ratio_pred, color = "Physical Contrast")) +
  geom_point(aes(y = Weber_ratio_pred, color = "Physical Contrast")) +
  geom_line(aes(y = Weber_ratio_perceived, color = "Perceived Contrast (log μ)")) +
  geom_point(aes(y = Weber_ratio_perceived, color = "Perceived Contrast (log μ)")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_classic() +
  ylab("Weber Ratio (ΔC / C)") +
  xlab("Contrast") +
  ggtitle("Physical vs Perceived Contrast Weber Ratio")