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
library(mvtnorm)
library(viridis)
library(data.table)
library(mgcv)
library(xgboost)

cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

# Naka-Rushton function
Rmax <- 115
C50 <- 19.3
n   <-  2.9
R0 <- Rmax * 0.03

# parameters
n_neurons <- 180
orientations <- c(45, 135)
preferred_orientations <- seq(1, 180, length.out = n_neurons)
tuning_width <- 40
n_trials <- 10000

# stimulus 
contrasts <- seq(0, 10, 1)
stim_angles <- c(45, 135)

results_list <- list()

for (cc in contrasts) {
  
  max_firing_rate <- Rmax * cc^n / (cc^n + C50^n)
  
  tuning_curves <- matrix(0, nrow = n_neurons, ncol = length(orientations))
  
  for (i in 1:n_neurons) {
    tuning_curves[i, ] <- max_firing_rate * exp(
      -0.5 * (pmin(abs(orientations - preferred_orientations[i]),
                   180 - abs(orientations - preferred_orientations[i]))^2) /
        tuning_width^2
    )
  }
  
  for (stim_angle in stim_angles) {
    
    input_stim <- rnorm(n_trials, stim_angle, 0)
    
    responses <- foreach(i = 1:n_trials,
                         .combine = rbind,
                         .packages = "stats") %dopar% {
                           idx <- ifelse(round(input_stim[i]) == 45, 1, 2)
                           
                           g <- 1
                           mu <- g * tuning_curves[, idx]
                           mu[mu == 0] <- 1e-8
                           
                           resp <- rnbinom(n_neurons, size = 1000000, mu = mu) +
                             rnbinom(n_neurons, size = 1000000, mu = R0)
                           
                           cbind(resp,
                                 neuron = seq_len(n_neurons),
                                 stim = stim_angle,
                                 contrast = cc,
                                 trial = i)
                         }
    key <- paste(cc, stim_angle, sep = "_")
    results_list[[key]] <- as.data.frame(responses)
  }
}

df_responses <- do.call(rbind, results_list)
colnames(df_responses) <- c("Spikes", "Neuron", "Stimulus", "Contrast", "Trial")
df_responses <- df_responses %>%
  mutate(Target = ifelse(Stimulus == 45, 1, 0))
colSums(is.na(df_responses))
any(is.na(df_responses))
# fwrite(df_responses, "df_responses.csv")
# df_responses <- fread("df_responses.csv", header = T)

# LFI
compute_lfi_decomposition <- function(data, lambda = 1e-6) {
  
  neurons <- grep("^[0-9]+$", colnames(data), value = TRUE)
  p <- length(neurons)
  
  data0 <- data[data$Target == 0, neurons]
  data1 <- data[data$Target == 1, neurons]
  
  mu0 <- colMeans(data0)
  mu1 <- colMeans(data1)
  
  df <- mu1 - mu0
  
  Sigma <- (cov(data0) + cov(data1)) / 2
  Sigma <- Sigma + lambda * diag(p)
  
  # Fisher decoder
  w <- solve(Sigma) %*% df
  
  contribution <- as.numeric(df * w)
  
  # ------------------
  # Trial-wise prediction
  # ------------------
  
  X <- as.matrix(data[, neurons])
  
  # 判別スコア
  score <- as.numeric(X %*% w)
  
  # クラス平均の中点を閾値にする
  threshold <- as.numeric(((mu0 + mu1) / 2) %*% w)
  
  predicted <- ifelse(score > threshold, 1, 0)
  
  prediction_df <- data.frame(
    Trial = data$Trial,
    Target = data$Target,
    Score = score,
    Predicted = predicted,
    Correct = predicted == data$Target
  )
  
  # 独立仮定
  Sigma_diag <- diag(diag(Sigma))
  w_shuffled <- solve(Sigma_diag) %*% df
  I_shuffled <- as.numeric(t(df) %*% w_shuffled)
  
  list(
    total_I = sum(contribution),
    neuron_contributions = contribution,
    shuffled_I = I_shuffled,
    correlation_gain = sum(contribution) / I_shuffled,
    lfi_weights = w,
    predictions = prediction_df
  )
}

df_wide <- df_responses %>%
  pivot_wider(
    id_cols = c(Stimulus, Contrast, Trial, Target),
    names_from = Neuron,
    values_from = Spikes
  )

colSums(is.na(df_wide))
any(is.na(df_wide))

lfi <- compute_lfi_decomposition(df_wide)
lfi

# xgboost
df_temp <- left_join(lfi$predictions, df_wide)
neuron_cols <- as.character(1:180)
X_data <- as.matrix(df_temp[, neuron_cols])
Y_data <- as.integer(df_temp$Correct)
dtrain <- xgb.DMatrix(data = X_data, label = Y_data)

params <- list(
  objective = "binary:logistic",  # 二値分類で「確率（確信度）」を出力する設定
  eval_metric = "logloss",        # 評価指標（対数損失）
  max_depth = 15,                 # 木の深さ。180個のニューロンの複雑な同時発火パターンを捉える
  min_child_weight = 2,           # 1桁の試行数でしか成り立たない特殊すぎるパターンを排除
  eta = 0.1,                      # 学習率
  gamma = 0                       # ブレーキを最小限にする
)

model_xgb <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 500,
  verbose = 1
)

df_temp$Confidence <- predict(model_xgb, dtrain)
head(df_temp[, c("Trial", "Correct", "Confidence")])

temp <- df_temp %>%
  ggplot(aes(x = `45`, y = `135`, z = Confidence)) +
  stat_summary_2d(
    fun = mean,
    bins = 15
  ) +
  scale_fill_viridis_c(
    limits = c(0.7, 1)
  ) +
  coord_equal() +
  labs(
    fill = "Proportion(45 deg)",
    x = "Neuron 45",
    y = "Neuron 135"
  ) +
  theme_minimal()
temp


df_mean_conf <- df_temp %>%
  group_by(`45`, `135`) %>%
  summarise(
    Confidence = mean(Confidence),
    .groups = "drop"
  )
ggplot(df_mean_conf,
       aes(`45`, `135`, z = Confidence)) +
  geom_contour_filled()


df_prediction <- df_temp %>%
  group_by(`45`, `135`) %>%
  summarise(
    P_45 = mean(Predicted),
    .groups = "drop"
  )
ggplot(df_prediction,
       aes(`45`, `135`, z = P_45)) +
  geom_contour_filled()


df_mean_score <- df_temp %>%
  group_by(`45`, `135`) %>%
  summarise(
    Score = mean(Score),
    .groups = "drop"
  )
ggplot(df_mean_score,
       aes(`45`, `135`, z = Score)) +
  geom_contour_filled()


df_accuracy <- df_temp %>%
  group_by(`45`, `135`) %>%
  summarise(
    Accuracy = mean(Correct),
    .groups = "drop"
  )
ggplot(df_accuracy,
       aes(`45`, `135`, z = Accuracy)) +
  geom_contour_filled() 


df_mean <- df_temp %>%
  group_by(`45`, `135`) %>%
  summarise(
    Confidence = mean(Confidence),
    .groups = "drop"
  )
ggplot(df_mean,
       aes(`45`, `135`, z = Confidence)) +
  geom_contour_filled()


ggsave(
  filename = "confidence_template.png",
  plot = temp,
  width = 6,
  height = 5,
  dpi = 300
)