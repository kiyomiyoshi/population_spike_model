library(MASS)
library(mvtnorm)

n <- 1000
mu0 <- c(1, 0); Sigma0 <- matrix(c(1.5, 0, 0, 1), 2, 2)
mu1 <- c(0, 1); Sigma1 <- matrix(c(1,   0, 0, 1.5), 2, 2)

x0 <- mvrnorm(n, mu0, Sigma0)
x1 <- mvrnorm(n, mu1, Sigma1)

X <- rbind(x0, x1)
y <- c(rep(0, n), rep(1, n))
df <- data.frame(x1 = X[, 1], x2 = X[, 2], y = factor(y))

x_seq <- seq(-4, 6, length.out = 200)
y_seq <- seq(-4, 6, length.out = 200)
grid <- expand.grid(x1 = x_seq, x2 = y_seq)

prior0 <- 0.5
prior1 <- 0.5

grid$px_C0 <- dmvnorm(as.matrix(grid[, c("x1","x2")]), mean = mu0, sigma = Sigma0)
grid$px_C1 <- dmvnorm(as.matrix(grid[, c("x1","x2")]), mean = mu1, sigma = Sigma1)
grid$prob <- (grid$px_C1 * prior1) / (grid$px_C0*prior0 + grid$px_C1*prior1)
prob_matrix <- matrix(grid$prob, nrow = length(x_seq))

bayesian <- filled.contour(
  x_seq, y_seq, prob_matrix,
  color.palette = colorRampPalette(c("#c6dbef", "#f7f7f7", "#fcbba1")),
  xlab = "Neuron 1 firing",
  ylab = "Neuron 2 firing",
  main = "Logistic Regression Probability",
  plot.axes = {
    axis(1); axis(2)
    
    points(df$x1, df$x2,
           col = ifelse(df$y == 0, "blue", "red"),
           pch = 16, cex = 0.6)
    
    contour(x_seq, y_seq, prob_matrix,
            levels = 0.5, add = TRUE, lwd = 2)
  }
)

bayesian