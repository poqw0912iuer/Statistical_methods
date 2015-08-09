source("probBhelp.R")
source("probBdata.R")
residualSampling <- function(data, B){
   betas <- matrix(nrow = B, ncol = 4)
   for(i in 1:B){
      ARp.beta.est(data, 2)$LS -> betaLS
      ARp.beta.est(data, 2)$LA -> betaLA
      ARp.resid(data, betaLS) -> e.observedLS
      ARp.resid(data, betaLA) -> e.observedLA
      sample(e.observedLS, size=100, replace=TRUE) -> eLS
      sample(e.observedLA, size=100, replace=TRUE) -> eLA
      idx = sample(1:99, 1)
      x0 <- c(data[idx], data[idx+1])
      ARp.filter(x0, betaLS, eLS) -> xLS
      ARp.filter(x0, betaLA, eLA) -> xLA 
      ARp.beta.est(xLS, 2)$LS -> betaLS
      ARp.beta.est(xLA, 2)$LA -> betaLA
      betas[i, 1] = betaLS[1]
      betas[i, 2] = betaLS[2]
      betas[i, 3] = betaLA[1]
      betas[i, 4] = betaLA[2]
   }
   betas
}

ARp.beta.est(data3A$x, 2)$LS -> betaLS_Original
ARp.beta.est(data3A$x, 2)$LA -> betaLA_Original

set.seed(9925)
residualSampling(data3A$x, 10000) -> betas
par(mfrow=c(2,2))

truehist(betas[,1], prob=TRUE, ylab="Density", 
xlab=expression(beta*"1 LS"), col="lightblue")
abline(v=betaLS_Original[1], col=2, lwd=3)

truehist(betas[,2], prob=TRUE, ylab="Density", 
xlab=expression(beta*"2 LS"), col="lightblue")
abline(v=betaLS_Original[2], col=2, lwd=3)

truehist(betas[,3], prob=TRUE, ylab="Density", 
xlab=expression(beta*"1 LA"), col="lightblue")
abline(v=betaLA_Original[1], col=2, lwd=3)

truehist(betas[,4], prob=TRUE, ylab="Density", 
xlab=expression(beta*"2 LA"), col="lightblue")
abline(v=betaLA_Original[2], col=2, lwd=3)

mean(betas[,1] - betaLS_Original[1])
mean(betas[,2] - betaLS_Original[2])
mean(betas[,3] - betaLA_Original[1])
mean(betas[,4] - betaLA_Original[2])
var(betas[,1])
var(betas[,2])
var(betas[,3])
var(betas[,4])



x101LS <- rep(0, 10000)
x101LA <- rep(0, 10000)
for(i in 1:10000){
   x101LS[i] = betas[i,1]*data3A$x[99]+betas[i,2]*data3A$x[100]
   x101LA[i] = betas[i,3]*data3A$x[99]+betas[i,4]*data3A$x[100]
}
# Plotting routine
par(mfrow=c(1,2))
truehist(x101LS, prob=TRUE, ylab="Density", 
xlab="x101_LS", col="lightblue")
abline(v=mean(x101LS), col=2, lwd=3)

truehist(x101LA, prob=TRUE, ylab="Density", 
xlab="x101_LA", col="lightblue")
abline(v=mean(x101LA), col=2, lwd=3)

quantile(x101LS,  probs = c(0.025, 0.975))
quantile(x101LA,  probs = c(0.025, 0.975))



# Not necessary
x1 <- rep(0, 101)
x2 <- rep(0, 101)
x1[1] = data3A$x[1]
x1[2] = data3A$x[2]
x2[1] = data3A$x[1]
x2[2] = data3A$x[2]
for(i in 1:98){
   x1[i+2] = mean(betas[,1])*data3A$x[i]+mean(betas[,2])*data3A$x[i+1]
   x2[i+2] = mean(betas[,3])*data3A$x[i]+mean(betas[,4])*data3A$x[i+1]
}
x1[101] = mean(x101LS)
x2[101] = mean(x101LA)
par(mfrow=c(2,1))
plot(data3A$x, type = "l", col=1)
lines(x1, col=2, lwd = 1)
plot(data3A$x, type = "l", col=1)
lines(x2, col=3, lwd = 1)

