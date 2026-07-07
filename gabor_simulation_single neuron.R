library(imager)
library(tidyverse)

set.seed(1)


# gabor function
make_gabor <- function(size=64,
                       theta=0,
                       sf=0.08,
                       sigma=12,
                       contrast=1){
  
  x <- seq(-size/2,size/2,length.out = size)
  y <- x
  
  grid <- expand.grid(x = x,y = y)
  
  xr <- grid$x*cos(theta*pi/180) + grid$y*sin(theta*pi/180)
  
  gabor <- contrast *
    exp(-(grid$x^2 + grid$y^2)/(2*sigma^2)) *
    cos(2*pi*sf*xr)
  
  matrix(gabor,nrow = size)
}


# orientation filter bank
angles <- seq(0, 170, 10)

filters <- lapply(angles, function(a){
  
  make_gabor(
    theta = a,
    sigma = 12,
    sf = 0.08
  )
  
})


# stimulus
stim_orientation <- 90

clean_img <- make_gabor(
  theta = stim_orientation,
  sigma = 12,
  sf = 0.08
)


# filter responses
calc_filter_response <- function(img, filter){
  
  sum(img * filter) 
  
}


# weights across different filters
make_weights <- function(width){
  
  w <- exp(
    -(angles - 90)^2 /
      (2*width^2)
  )
  
  w / sum(w) ##### 総入力を揃える #####

}


# trial
simulate_trial <- function(width,
                           noise_sd = 0.3){
  
  noise <- matrix(
    rnorm(length(clean_img),
          mean = 0,
          sd = noise_sd),
    nrow = nrow(clean_img)
  )
  
  img <- clean_img + noise
  img <- pmax(pmin(clean_img + noise, 1), -1) # クリップ
  
  channel_response <- sapply(filters, function(f){
    
    calc_filter_response(img, f)
    
  })
  
  # tuning width
  weights <- make_weights(width)
  # neuron response
  sum(channel_response * weights)
  
}


# simulations
n_trials <- 1000

responses_narrow <- replicate(
  n_trials,
  simulate_trial(width = 10)
)

responses_wide <- replicate(
  n_trials,
  simulate_trial(width = 20)
)

fano <- function(x){
  
  var(x) / mean(x)
  
}

FF_narrow <- fano(responses_narrow)
FF_wide <- fano(responses_wide)
FF_narrow
FF_wide

df <- data.frame(
  
  response = c(
    responses_narrow,
    responses_wide
  ),
  
  condition = rep(
    c("narrow",
      "wide"),
    each = n_trials
  )
  
)

g1 <- ggplot(df,
       aes(response,
           fill = condition)) +
  geom_histogram(
    alpha = 0.5,
    bins = 40,
    position = "identity"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0(
      "FF narrow = ",
      round(FF_narrow, 3),
      "\nFF wide = ",
      round(FF_wide, 3)
    ),
    hjust = 1.1,
    vjust = 1.5,
    size = 3.5
  ) +
  theme_classic()
g1


# tuning curves
calc_tuning_curve <- function(width){
  
  weights <- make_weights(width)
  
  responses <- sapply(angles, function(stim_angle){
    
    img <- make_gabor(
      theta = stim_angle,
      sigma = 12,
      sf = 0.08
    )
    
    channel_response <- sapply(filters, function(f){
      
      sum(img * f)
      
    })
    
    # pooling
    sum(channel_response * weights)
    
  })
  
  responses
}

curve_narrow <- calc_tuning_curve(width = 10)
curve_wide <- calc_tuning_curve(width = 20)

tuning_df <- data.frame(
  
  angle = rep(angles, 2),
  
  response = c(
    curve_narrow,
    curve_wide
  ),
  
  condition = rep(
    c(
      "narrow (sigma_theta = 10)",
      "wide (sigma_theta = 20)"
    ),
    each = length(angles)
  )
  
)

g2 <- ggplot(
  tuning_df,
  aes(
    x = angle,
    y = response,
    color = condition
  )
) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(
    breaks = seq(0, 170, 30)
  ) +
  theme_classic() +
  labs(
    x = "Stimulus orientation (degree)",
    y = "Neuron response",
    color = "Tuning width"
  )
g2

area_narrow <- sum(curve_narrow)
area_wide <- sum(curve_wide)
area_narrow
area_wide