sigma_g <- 0.2
g <- rgamma(10000, shape = 1 / sigma_g^2, scale = sigma_g^2)

summary(g)
sd(g)

hist(g,
     breaks = 50,
     col = "skyblue",
     border = "white")