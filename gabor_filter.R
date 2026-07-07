library(imager)
library(tidyverse)

make_gabor <- function(size     = 64,
                       theta    = 0,
                       sf       = 0.08,
                       sigma    = 12,
                       contrast = 1){
  
  x <- seq(-size/2, size/2, length.out = size)
  y <- x
  
  grid <- expand.grid(x = x,y = y)
  
  xr <- grid$x*cos(theta*pi/180)　+　grid$y*sin(theta*pi/180)
  
  gabor <- contrast *
    exp(-(grid$x^2　+　grid$y^2)/(2*sigma^2)) *
    cos(2*pi*sf*xr)
  
  matrix(gabor,　nrow　=　size)
}


add_noise <- function(img,
                      noise_sd　=　0.4){
  
  img + matrix(rnorm(length(img),　0,　noise_sd),
               nrow　=　nrow(img))
}


img <- make_gabor(theta　=　90)
img <- add_noise(img, noise_sd　=　0.3)
img <- pmin(pmax(img, -1), 1)

df <- expand.grid(
  x = 1:ncol(img),
  y = 1:nrow(img)
)
df$z <- as.vector(img)

ggplot(df, aes(x, y, fill = z)) +
  geom_raster() +
  scale_fill_gradient2(low = "black",
                       mid = "gray",
                       high = "white") +
  coord_equal() +
  theme_void()


angles <- seq(0, 170, 10)
filters <- lapply(angles, function(a){
  
  make_gabor(theta = a)
  
})


responses <- sapply(filters, function(f){
  
  sum(img*f)
  
})


# Tuning function
# 受容野（90°選択）
RF <- make_gabor(theta = 90)
# 刺激方位
stim_angles <- seq(0, 170, 10)
# 各刺激方位に対する応答
tuning_response <- sapply(stim_angles, function(a){
  
  stimulus <- make_gabor(theta = a)
  sum(stimulus * RF)
  
})

tuning_df <- data.frame(
  angle = stim_angles,
  response = tuning_response
)

ggplot(tuning_df,
       aes(x = angle,
           y = response)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  theme_classic() +
  xlab("Stimulus orientation (deg)") +
  ylab("Filter response") +
  ggtitle("Orientation tuning curve")


# Simulation
RF_narrow <- make_gabor(
  theta = 90,
  sigma = 5,
  sf = 0.08
)

RF_wide <- make_gabor(
  theta = 90,
  sigma = 20,
  sf = 0.08
)


calc_tuning <- function(RF){
  
  sapply(stim_angles,function(a){
    
    stimulus <- make_gabor(theta = a)
    sum(stimulus*RF)
    
  })
  
}

response_narrow <- calc_tuning(RF_narrow)
response_wide <- calc_tuning(RF_wide)

#########
calc_response <- function(img, filter){
  
  sum(img * filter) /
    sqrt(sum(img^2) * sum(filter^2))
  
}

response_narrow <- sapply(stim_angles,function(a){
  
  stimulus <- make_gabor(theta=a)
  
  calc_response(stimulus, RF_narrow)
  
})


response_wide <- sapply(stim_angles,function(a){
  
  stimulus <- make_gabor(theta=a)
  
  calc_response(stimulus, RF_wide)
  
})

#########


plot_df <- data.frame(
  
  angle = rep(stim_angles,2),
  response = c(response_narrow,
               response_wide),
  condition = rep(c("narrow",
                    "wide"),
                  each = length(stim_angles))
)



ggplot(plot_df,
       aes(angle,response,
           color = condition)) +
  geom_line(size = 1.2) +
  geom_point() +
  theme_classic() +
  xlab("Orientation") +
  ylab("Response")