#PROBLEM A
#Exponential distribution
rExpo <- function(n,lambda){ 
  u <- runif(n)
  -log(u)/lambda  
}

#Schape of Exp. Dist.
x <- rExpo(n=10000,lambda=0.5)
hist(x,freq=FALSE)
t <- 0:1000/1
lines(t,dexp(t,rate=0.5),lwd=2)

#Problem A2 c)
#Function that generate from pdf "f"
rF <- function(n,alpha){ 
   u <- runif(n)
   -log(1/u-1)/alpha
}

#Probability denisity function
f <- function(t, alpha){
    alpha*exp(alpha*t)/(1+exp(alpha*t))^2
}

#Plots of denisities
x <- rF(100000,3.2)
hist(x,freq=FALSE)
t <- -100:100/10
lines(t,f(t,3.2),lwd=2)

y <- f(t, 2)

Var<- function(t, alpha){  
    v=0
    v = sum(t^2*f(t, alpha))
    return(v)
}
Var(t, 2)


rG <- function(alpha, n)
{
   x = 1:n
   u = 1:n
   u[i] <- runif(n)
   for(i in 1:n){
   if(any(u < (exp(1)/(exp(1) + alpha)))){ 
      x[i] <- ((exp(1) + alpha) / exp(1) * u)^(1/alpha)
   }
   else
       x[i] <- -log( 1/(1 - u), base = exp(1)) + log(exp(1) * alpha/(exp(1) + alpha), base = exp(1)) 
   }   
   hist(x, freq = FALSE)
}
rG(0.3, 100)
hist(x,freq=FALSE)

g <- function(alpha){
   c = (alpha+exp(1))/(alpha*exp(1))
   curve( (x<1)*(c*x^(alpha-1)) + (!(x<1))*c*exp(-x), 0, 10)
}