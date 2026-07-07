library(imager)
library(tidyverse)
library(patchwork)

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

preferred_angles <- seq(0,170,10)

make_orientation_weights <- function(pref, width){
  
  w <- exp(
    -(angles-pref)^2 /
      (2*width^2)
  )
  
  w / sum(w)
  
}


make_population_weights <- function(width){
  
  lapply(
    preferred_angles,
    function(pref){
      
      make_orientation_weights(
        pref,
        width
      )
      
    }
  )
  
}


weights_narrow <- make_population_weights(10)
weights_wide <- make_population_weights(20)


simulate_population_trial <- function(
    weights,
    noise_sd=0.5
){
  
  noise <- matrix(
    rnorm(
      length(clean_img),
      mean=0,
      sd=noise_sd
    ),
    nrow=nrow(clean_img)
  )
  
  img <- clean_img + noise
  
  
  channel_response <- sapply(
    filters,
    function(f){
      
      sum(img*f)
      
    }
  )
  
  
  neuron_response <- sapply(
    weights,
    function(w){
      
      sum(channel_response*w)
      
    }
  )
  
  
  neuron_response
  
}



n_trials <- 3000


responses_narrow_pop <- replicate(
  n_trials,
  simulate_population_trial(
    weights_narrow
  )
)


responses_wide_pop <- replicate(
  n_trials,
  simulate_population_trial(
    weights_wide
  )
)

dim(responses_narrow_pop)
dim(responses_wide_pop)

fano <- function(x){
  var(x)/mean(x)
}


FF_narrow_each <- apply(
  responses_narrow_pop,
  1,
  fano
)


FF_wide_each <- apply(
  responses_wide_pop,
  1,
  fano
)


ff_df <- data.frame(
  
  angle = rep(
    preferred_angles,
    2
  ),
  
  FF = c(
    FF_narrow_each,
    FF_wide_each
  ),
  
  width = rep(
    c(
      "narrow (10)",
      "wide (20)"
    ),
    each=18
  )
  
)


g3 <- ggplot(
  ff_df,
  aes(
    x=angle,
    y=FF,
    color=width
  )
)+
  geom_line(size=1)+
  geom_point(size=2)+
  scale_x_continuous(
    breaks=seq(0,170,30)
  )+
  scale_y_log10() +
  theme_classic()+
  labs(
    x="Preferred orientation",
    y="Fano factor (log scale)"
  )
g3

population_sum_narrow <- colSums(
  responses_narrow_pop
)

FF_population_narrow <- fano(
  population_sum_narrow
)

population_sum_wide <- colSums(
  responses_wide_pop
)

FF_population_wide <- fano(
  population_sum_wide
)

FF_population_narrow
FF_population_wide

pop_df <- data.frame(
  
  response=c(
    population_sum_narrow,
    population_sum_wide
  ),
  
  width=rep(
    c(
      "narrow (10)",
      "wide (20)"
    ),
    each=n_trials
  )
  
)


g4 <- ggplot(
  pop_df,
  aes(
    response,
    fill=width
  )
)+
  geom_histogram(
    bins=40,
    alpha=0.5,
    position="identity"
  )+
  annotate(
    "text",
    x=Inf,
    y=Inf,
    label=paste0(
      "FF narrow = ",
      round(FF_population_narrow,3),
      "\nFF wide = ",
      round(FF_population_wide,3)
    ),
    hjust=1.1,
    vjust=2,
    size=3
  )+
  theme_classic()
g4

g_combined <- g3 + g4 + 
  plot_layout(ncol = 2)

g_combined

ggsave(
  "FF_comparison.png",
  g_combined,
  width = 8,
  height = 3.5,
  dpi = 300
)