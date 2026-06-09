library(mgcv)
library(dplyr)
library(tidyr)
library(akima)
library(ggplot2)

###
# 平均Confidenceを計算
df_mean <- df_temp %>%
  group_by(`45`, `135`) %>%
  summarise(
    Confidence = mean(Confidence),
    .groups = "drop"
  )

# 補間
interp_res <- akima::interp(
  x = df_mean$`45`,
  y = df_mean$`135`,
  z = df_mean$Confidence,
  xo = seq(0, 45, length = 30),
  yo = seq(0, 45, length = 30)
)

interp_df <- expand.grid(
  `45` = interp_res$x,
  `135` = interp_res$y
)

interp_df$Confidence <- as.vector(interp_res$z)

ggplot(interp_df,
       aes(`45`, `135`, fill = Confidence)) +
  geom_raster(interpolate = TRUE) +
  coord_equal() +
  scale_fill_viridis_c(
    limits = c(0.85, 0.95)
  ) +
  theme_minimal()


###
fit <- gam(
  Confidence ~ s(`45`, `135`),
  data = df_temp
)

grid <- expand.grid(
  `45` = seq(0, 15, length = 100),
  `135` = seq(0, 15, length = 100)
)

grid$pred <- predict(fit, newdata = grid)

ggplot(grid,
       aes(`45`, `135`, fill = pred)) +
  geom_raster(interpolate = TRUE) +
  coord_equal() +
  scale_fill_viridis_c()


###
ggplot(df_mean,
       aes(`45`, `135`, z = Confidence)) +
  geom_contour_filled()
