par(mfrow=c(1,1))
bilirubin <- read.table("bilirubin.txt",header=T)
boxplot(meas~pers,data=bilirubin, main="Bilirubin log-concentration Data", 
  	 xlab="Individual", ylab="bilirubin log-concentration (mg/dL)")


x1 <- c(rep(1, 11), rep(0, 18))
x2 <- c(rep(0, 11), rep(1, 10), rep(0, 8))
x3 <- c(rep(0, 21), rep(1, 8))
summary(lm(log(bilirubin[,1])~x1+x2+x3-1)) -> fit
fit$fstatistic[1] -> Fval
fit


permTest <- function(){
   N = length(bilirubin[,1])
   sample(1:N, N, replace = F) -> idx
   y <- bilirubin[idx,1]
   summary(lm(log(y)~x1+x2+x3-1))$fstatistic[1]
}

B = 999
tot = 0
F_data = rep(0, B)
for(i in 1:B){
   F <- permTest()
   if ( F > Fval)
   {
      tot = tot + 1
   }
   F_data[i] = F
}
tot = tot/B
tot
library(MASS)
truehist(F_data)
points(Fval,0,col=2, pch=20)

