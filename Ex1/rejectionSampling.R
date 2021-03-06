#Problem C1
rnorm1 <- function(n){
  u = runif(n)
  v = runif(n)
  x=rep(0,n)

  for (i in 1:n){
    x[i] = sqrt(-2*log(u[i]))*cos(2*pi*v[i])  
  }
}
truehist(rnorm(100000))
t <- -100:100/10
lines(t,dnorm(t,mean=0,sd=1),lwd=2)

#Proble C2 a)
truehist(rInvertedG(10000, 0.7))
c <- 2
t <- 0:1000/100
lines(t,dgamma(t,1,1),lwd=2) 
lines(t, c * g_pdf(t, 0.7), lwd=2)

#visualisation of target and proposal denisities

c<-2
t <- 0:1000/100
plot (t,dgamma(t,1,1), type="l", col="red")
lines(t, c*g_pdf(t, 0.7), col="blue", lwd=1)
legend("topright", lty = c(1,1),lwd=c(3,3), col=c("blue","red"), c("Target denisity", "c * Proposal denisity"))

#Proble C2 b)
#Rejection sampling from Gamma
rGammaReject <- function(n, alpha){
    c <- 2
    y <- numeric(n) 
    accepted.samples <- 0 
    while (accepted.samples<n) {
       x <- rInvertedG(1, 0.7) 
       f.at.x <- dgamma(x,alpha,1) 
       g.at.x <- g_pdf(x,0.7) 
       accept <- f.at.x / (c*g.at.x) 
       if (runif(1)<accept) { 
          accepted.samples <- accepted.samples+1
          y[accepted.samples] <- x  
       }
    }
    return(y)
}

#Test comparition
y <- rGammaReject(30000, 0.6)
truehist(y) 
t <- 0:700/100
lines(t,dgamma(t,0.6,1), lwd=2)
mean(y)
var(y)

#Proble C3 b)
#Ratio of uniform method for Gamma
rGamma3 <- function(n, alpha){
   y <- numeric(n)
   accepted.samples <- 0
   a <- ((alpha - 1)/exp(1))^(0.5*alpha - 0.5)
   b_plus <- ((alpha+1)/exp(1))^(0.5*alpha +0.5)
   while(accepted.samples < n){
      u <- a*runif(1)
      v <- b_plus*runif(1)
      x <- v/u
      accepted <- (alpha - 1)*log(x) - x
      if(2*log(u) < accepted){
            accepted.samples <- accepted.samples + 1
            y[accepted.samples] <- x
      } 
   }
   return(y)
}
#Test comparition
library(MASS) 
y <- rGamma3(1000, 4) 
truehist(y)   
t <- 0:1500/100
lines(t,dgamma(t,4,1), lwd=2)
mean(y)
var(y)

#Proble C4
#Gamma distribution (alpha > 0, beta > 0)
rGammaGeneral <- function(n, alpha, beta){
   smaple <- numeric(n)
   if(alpha > 0){
      if(alpha < 1)   smaple <- 1/beta*rGamma1(n, alpha)
      else   smaple <- 1/beta*rGamma3(n, alpha)
   }
   return(smaple)
}

y <- rGammaGeneral(30000, .5, 1) 
truehist(y) 
t <- 0:10000/10
lines(t,dgamma(t,.7, 1), lwd=3)
legend("topright",  "Ga(0.7, 11)")