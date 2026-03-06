library(ggplot2)
library(patchwork)

make_gabor <- function(size = 256, sigma = 30, freq = 0.03,
                       phase = 0, theta = 0, contrast = 1){
  
  x <- seq(-size/2, size/2, length.out = size)
  y <- seq(-size/2, size/2, length.out = size)
  grid <- expand.grid(x = x, y = y)
  
  x_theta <- grid$x*cos(theta) + grid$y*sin(theta)
  y_theta <- -grid$x*sin(theta) + grid$y*cos(theta)
  
  gaussian <- exp(-(x_theta^2 + y_theta^2)/(2*sigma^2))
  sinusoid <- cos(2*pi*freq*x_theta + phase)
  grid$gabor <- contrast * gaussian * sinusoid
  
  return(grid)
}

p1 <- ggplot(make_gabor(theta = pi/4, contrast = 0.1), aes(x, y, fill = gabor)) +
  geom_raster() +
  scale_fill_gradient2(low = "black", mid = "gray", high = "white", limits = c(-1, 1)) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none")

p2 <- ggplot(make_gabor(theta = -pi/180, contrast = 1), aes(x, y, fill = gabor)) +
  geom_raster() +
  scale_fill_gradient2(low = "black", mid = "gray", high = "white", limits = c(-1, 1)) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none")

g1 <- list(p1, p2)
g1 <- lapply(g1, function(p) {
  p + theme(
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )
})

gabor <- wrap_plots(g1, nrow = 1)
ggsave("gabor.png", gabor, width = 4, height = 2, dpi = 300, bg = "white")