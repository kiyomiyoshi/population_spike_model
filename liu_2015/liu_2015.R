library(tidyverse)

df_liu <- data.frame(max_firing = c(20, 30, 42, 51, 60, 70, 71,  73),
                     contrast   = c( 5, 16, 19, 31, 45, 50, 70, 100))

# power
fit <- nls(max_firing ~ k * contrast^n, data = df_liu, start = list(k = 1, n = 0.5))
summary(fit)

x_new <- seq(0, 100, length.out = 200)
params <- coef(fit)
y_new <- params["k"] * x_new^params["n"]
df_fit <- data.frame(contrast = x_new, max_firing = y_new)

ggplot(df_liu, aes(x = contrast, y = max_firing)) +
  geom_point(size = 3) +  
  geom_line(data = df_fit, aes(x = contrast, y = max_firing), color = "blue", size = 1) +
  coord_cartesian(xlim = c(0, 100)) +
  theme_minimal() +
  labs(title = "Stevens power law",
       x = "Contrast",
       y = "Max Firing Rate") -> g1

# log
fit_log <- nls(max_firing ~ a * log(contrast) + b,
               data = df_liu,
               start = list(a = 10, b = 0))
summary(fit_log)

x_new <- seq(1, 100, length.out = 200)
params_log <- coef(fit_log)
y_new <- params_log["a"] * log(x_new) + params_log["b"]
df_fit_log <- data.frame(contrast = x_new, max_firing = y_new)

ggplot(df_liu, aes(x = contrast, y = max_firing)) +
  geom_point(size = 3) +
  geom_line(data = df_fit_log, aes(x = contrast, y = max_firing), color = "green", size = 1) +
  coord_cartesian(xlim = c(0, 100)) +
  theme_minimal() +
  labs(title = "Log",
       x = "Contrast",
       y = "Max Firing Rate") -> g2

# log squared
fit_log_squared <- nls(max_firing ~ a * (log(contrast) + b)^2,
               data = df_liu,
               start = list(a = 0.2, b = 0))
summary(fit_log_squared)

x_new <- seq(1, 100, length.out = 200)
params_log <- coef(fit_log_squared)
y_new <- params_log["a"] * (log(x_new) + params_log["b"])^2 
df_fit_log_squared <- data.frame(contrast = x_new, max_firing = y_new)

ggplot(df_liu, aes(x = contrast, y = max_firing)) +
  geom_point(size = 3) +
  geom_line(data = df_fit_log_squared, aes(x = contrast, y = max_firing), color = "red", size = 1) +
  coord_cartesian(xlim = c(0, 100)) +
  theme_minimal() +
  labs(title = "Log squared",
       x = "Contrast",
       y = "Max Firing Rate") -> g3

# Naka-Rushton function
fit_naka <- nls(max_firing ~ Rmax * contrast^n / (contrast^n + C50^n),
                data = df_liu,
                start = list(Rmax = 75, C50 = 25, n = 2))
summary(fit_naka)
x_new <- seq(0, 100, length.out = 200)
params <- coef(fit_naka)
y_new <- params["Rmax"] * x_new^params["n"] / (x_new^params["n"] + params["C50"]^params["n"])
df_fit_naka <- data.frame(contrast = x_new, max_firing = y_new)

ggplot(df_liu, aes(x = contrast, y = max_firing)) +
  geom_point(size = 3) +
  geom_line(data = df_fit_naka, aes(x = contrast, y = max_firing), color = "navy", size = 1) +
  coord_cartesian(xlim = c(0, 100)) +
  theme_minimal() +
  labs(title = "Naka-Rushton",
       x = "Contrast",
       y = "Max Firing Rate") -> g4

g <- cowplot::plot_grid(g1, g2, g3, g4, nrow = 1)
ggsave("liu_2015.png", g, width = 8, height = 2, dpi = 300)