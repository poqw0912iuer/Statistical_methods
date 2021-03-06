#Problem D1
MVrnorm <- function(n, mu, sigma) {
    ncols <- ncol(sigma)
    m <- rep(mu, each = n) 
    m + matrix(rnorm(n * ncols), ncol = ncols) %*% chol(sigma)
}

sigma <- matrix(c(1, -0.2, 0.7, -0.2, 1, -0.1, 0.7, -0.1, 1), ncol = 3)
mu <- c( 7, -12, 4)

x <- mvrnorm(100000, mu,  sigma)
cor(x)
colMeans(x)

#comparision of integrate multivariate normal distribution
x = mvrnorm(100000, mu, sigma)
var(x)
colMeans(x)

sigma <- matrix(c(10,0,0,10),2,2)
mu <- c(1,2)
y <- MVrnorm(100000, mu,  sigma)
x = mvrnorm(100000, mu, sigma)

# use a kernel density estimator to plot the distribution
library(gplots)
y.kde = kde2d(y[,1], y[,2], n=100)
x.kde = kde2d(x[,1], x[,2], n=100)

par(mfrow=c(1,2))
image(y.kde, xlim=c(-10,10), ylim=c(-10,10))
contour(y.kde, add = T)
image(x.kde, xlim=c(-10,10), ylim=c(-10,10))
contour(x.kde, add = T)

#PRoblem D2 b)
#Direchlet dist.
rdirichlet1 = function(n, alpha) {
  k = length(alpha)
  r = matrix(0, nrow=n, ncol=k) 
  for (i in 1:k) {
    r[,i] = rGammaGeneral(n, alpha[i], 2) 
  }
  r <- matrix(mapply(function(r, s) {return (r/s)}, r, rowSums(r)), ncol=k)
  return(r)
}
r <- rdirichlet1(10000, c(1,2,3)) 
plot(r, main="alpha=(1, 2, 3)")

mean(r[,1])
mean(r[,2])
mean(r[,3])
var(r)