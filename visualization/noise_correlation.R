# 2つのニューロンが「同じ方向に揺れる」ため、総和方向（信号方向）での変動が増える
# ニューロンが独立なら、一方が大きく揺れても、もう一方が逆方向に揺れることもあるので、総和方向では揺れがある程度打ち消される
# これによって総和方向の分散は抑えられる

# ======== データ生成 ========
set.seed(42)

# ニューロン数
n_neurons <- 2

# 試行数
n_trials <- 1000

# クラス平均 (2ニューロンの総和方向)
mu1 <- c(1,1)
mu2 <- c(-1,-1)

# 共分散行列の例
rho <- 0.2  # ノイズ相関
sigma <- 1
Sigma <- matrix(c(sigma^2, rho*sigma^2,
                  rho*sigma^2, sigma^2), nrow=2)

# データ生成
library(MASS)
class1 <- mvrnorm(n_trials, mu = mu1, Sigma = Sigma)
class2 <- mvrnorm(n_trials, mu = mu2, Sigma = Sigma)

# クラスラベル
labels <- c(rep(1, n_trials), rep(2, n_trials))
data <- rbind(class1, class2)

# ======== 散布図 ========
plot(data[,1], data[,2], col=ifelse(labels==1,"red","blue"),
     pch=16, xlab="Neuron 1", ylab="Neuron 2",
     main=paste("Scatter plot with rho =", rho))
legend("topright", legend=c("Class 1", "Class 2"), col=c("red","blue"), pch=16)

# ======== 線形判別情報 (LFI) と d' ========
# LFI = Δμ^T Σ^{-1} Δμ
delta_mu <- mu1 - mu2
LFI <- t(delta_mu) %*% solve(Sigma) %*% delta_mu

# d' = sqrt(LFI)
d_prime <- sqrt(LFI)

cat("Linear Fisher Information (LFI):", LFI, "\n")
cat("d':", d_prime, "\n")