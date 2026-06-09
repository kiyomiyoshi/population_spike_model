library(tidyverse)
library(xgboost)

df_conf <- left_join(lfi$predictions, df_wide)

# ----------------------------------------
# 1. データの整理と前処理
# ----------------------------------------

# ニューロンの列名（"1", "2", ..., "180"）を自動で指定
neuron_cols <- as.character(1:180)

# 特徴量（ニューロンの発火数）を行列に変換
X_data <- as.matrix(df_conf[, neuron_cols])

# 目的変数（Correct）を TRUE->1, FALSE->0 に変換
Y_data <- as.integer(df_conf$Correct)

# XGBoost専用のデータ型（DMatrix）に変換
dtrain <- xgb.DMatrix(data = X_data, label = Y_data)


# ----------------------------------------
# 2. モデルのパラメータ設定（非線形を記述し尽くす設定）
# ----------------------------------------
# 前述の通り、シミュレーション等のデータを限界まで記述し尽くす（非線形 gradient）ための設定です
params <- list(
  objective = "binary:logistic",  # 二値分類で「確率（確信度）」を出力する設定
  eval_metric = "logloss",        # 評価指標（対数損失）
  max_depth = 12,                 # 木の深さ。180個のニューロンの複雑な同時発火パターンを捉える
  eta = 0.1,                      # 学習率
  gamma = 0                       # ブレーキを最小限にする
)


# ----------------------------------------
# 3. モデルの訓練
# ----------------------------------------
# nrounds（学習回数）を多めにして、完全にデータを学習させます
model_xgb <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 300,
  verbose = 1 # 学習の進捗を表示
)


# ----------------------------------------
# 4. 確信度（正当確率）の予測
# ----------------------------------------
# 各トライアルに対する予測確率（0.0 〜 1.0）が計算されます
# これが「モデルが判断した正解への確信度」になります
df_conf$Confidence <- predict(model_xgb, dtrain)


# ----------------------------------------
# 5. 結果の確認
# ----------------------------------------
# 元のデータに「Confidence」列が追加されているか確認
head(df_conf[, c("Trial", "Correct", "Confidence")])

g <- df_conf %>%
  ggplot(aes(x = `45`, y = `135`, z = Confidence)) +
  stat_summary_2d(
    fun = mean,
    bins = 15
  ) +
  scale_fill_viridis_c(
    limits = c(0.85, 0.95)
  ) +
  coord_equal() +
  labs(
    fill = "Proportion(45 deg)",
    x = "Neuron 45",
    y = "Neuron 135"
  ) +
  theme_minimal()
g

ggsave(
  filename = "confidence_plot_peb_test.png",
  plot = g,
  width = 6,
  height = 5,
  dpi = 300
)