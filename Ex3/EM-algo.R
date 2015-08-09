y <- c(5, 10, 8, 0, 7, 12)
hist <- c(0, 0, 0, 0, 0, 0, 0)
mu <- mean(y) * 6 / 5
y22 <- 1
while(abs(y[4] - y22) > 1e-5){
   alpha1 <- (y[1] + y[3] + y[5]) / 3 - mu
   beta1 <- (y[1] + y[2]) / 2 - mu
   beta3 <- (y[5] + y[6]) / 2 - mu
   alpha2 <- -alpha1
   beta2 <- -beta1 - beta3
   y22 <- y[4]
   y[4] <- mu +alpha2 + beta2 
   rbind(hist, c(mu, alpha1, alpha2, beta1, beta2, beta3, y[4])) -> hist
   mu <- mean(y)
}


yrange<-range(c(hist[,1], hist[,2], hist[,3], hist[,4], hist[,5], hist[,6]))
plot(hist[2:length(hist[,1]), 1], type="l", ylim=yrange, col=1, lwd = 2, xlab = "EM iteration number", ylab = "Parameter values")
lines(hist[2:length(hist[,1]), 2], col=2, lwd = 2)
lines(hist[2:length(hist[,1]), 3], col=3, lwd = 2)
lines(hist[2:length(hist[,1]), 4], col=4, lwd = 2)
lines(hist[2:length(hist[,1]), 5], col=5, lwd = 2)
lines(hist[2:length(hist[,1]), 6], col=6, lwd = 2)
legend(0, 7, lty = c(1,1),lwd=c(3,3),col=c(1,2,3,4,5,6), c(expression(mu*"1"), expression(alpha*"1"), 
expression(alpha*"2"),expression(beta*"1"),expression(beta*"2"),expression(beta*"3")), horiz = TRUE)

plot(hist[2:length(hist[,1]), 7], type="l", col=1, lwd = 2, xlab = "EM iteration number", ylab = "y22")

Y = matrix(c(5, 10, 8, 7, 12), ncol = 1)
X = matrix(c(1, 1, 1, 1, 1,
	     1,-1, 1, 1,-1,
             1, 1,-1, 0, 0,
             0, 0,-1, 1, 1), nrow = 5) 
solve(t(X)%*%X)%*%(t(X)%*%Y)

