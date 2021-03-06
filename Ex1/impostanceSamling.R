#Problem B1
MassFunction <- function(x){
   x <- (2+x)^125*(1-x)^(18+20)*x^34
   return(x)
}
optimize(function(x)(2+x)^125*(1-x)^(18+20)*x^34, lower=0,upper=1, maximum= T)
plot(function(x)(2+x)^125*(1-x)^(18+20)*x^34)

#Rajection sampling algorithm
RejectionSampling <- function(n){
   y <- numeric(n)
   c <- 1.838839e+29 
   accepted.samples <- 0
   while(accepted.samples < n){
      u <- runif(1)
      if(c* runif(1) < MassFunction(u)){  
         accepted.samples <- accepted.samples+1
         y[accepted.samples] <- u
      }
   }
   return(y)
}

#Problem B2
integrate(function(x)(2+x)^125*(1-x)^(18+20)*x^34, lower=0, upper=1)

const <- 2.357695e+28 
y <- RejectionSampling(10000)
truehist(y)
t <- 0:500/100
lines(t, MassFunction(t)/const,lwd = 2)
abline(v = mean(y), col = "red", lwd = 3)

y <- RejectionSampling(10000)
mean(y)
var(y)
integrate(function(x)x*(2+x)^125*(1-x)^(18+20)*x^34/const, lower=0, upper=1)
0.6228062 with absolute error < 1.7e-07
integrate(function(x)(x-0.6228062 )^2*(2+x)^125*(1-x)^(18+20)*x^34/const, lower=0, upper=1)
0.002594932 with absolute error < 5e-08

#Problem B4
y <- RejectionSampling(10000)
mean(y)
w <- MassFunction(y)/dbeta(y,2,5)
sum <- sum(w)
posterior_mean <- 1/sum*sum(y*w)
posterior_mean