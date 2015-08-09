data(kyphosis)
scatterplot3d(kyphosis[, 3], kyphosis[, 2], kyphosis[, 4], color = c("blue", "red")[kyphosis$Kyphosis], pch = 15,
main="blue color - absent disease;  reed color - present disease", xlab="Number", ylab="Age", zlab="Start")

lda(Kyphosis ~ Age + Number + Start, kyphosis) -> zlda
qda(Kyphosis ~ Age + Number + Start, kyphosis) -> zqda


Kfold_CV <- function(K){
   tot = c(0, 0)
   names(tot) <- c("lda", "qda")
   n = length(kyphosis$Kyphosis)
   train <- seq(1:n)
   T = n %/% K
   rest = n %% K
   for(i in 1:T){
      lda(Kyphosis ~ ., kyphosis, subset = train[-(((i-1)*K+1):i*K)]) -> zlda
      qda(Kyphosis ~ ., kyphosis, subset = train[-(((i-1)*K+1):i*K)]) -> zqda
      for(j in 1:K){
         if(predict(zlda, kyphosis[(i-1)*K+j,])$class != kyphosis$Kyphosis[(i-1)*K+j]){
            tot[1] = tot[1] + 1
         }
         if(predict(zqda, kyphosis[(i-1)*K+j,])$class != kyphosis$Kyphosis[(i-1)*K+j]){
            tot[2] = tot[2] + 1
         }
      }
   }
   if(rest != 0){
      lda(Kyphosis ~ ., kyphosis, subset = train[-((T*K+1):(T*K+rest))]) -> zlda
      qda(Kyphosis ~ ., kyphosis, subset = train[-((T*K+1):(T*K+rest))]) -> zqda
      for(i in 1:rest){
         if(predict(zlda, kyphosis[T*K+i,])$class != kyphosis$Kyphosis[T*K+i]){
            tot[1] = tot[1] + 1
         }
         if(predict(zqda, kyphosis[T*K+i,])$class != kyphosis$Kyphosis[T*K+i]){
            tot[2] = tot[2] + 1
         }              
      }
   }   
   tot = tot / n
   tot
}

Kfold_CV(10)



Kfold_CV <- function(K, neighbours){
   n = length(kyphosis$Kyphosis)
   T = n %/% K
   rest = n %% K
   tot = 0
   for(i in 1:T){
      test <- cbind(kyphosis[((i-1)*K+1):(i*K),2], kyphosis[((i-1)*K+1):(i*K),3], kyphosis[((i-1)*K+1):(i*K),4])
      train <- cbind(kyphosis[-(((i-1)*K+1):(i*K)),2], kyphosis[-(((i-1)*K+1):(i*K)),3], kyphosis[-(((i-1)*K+1):(i*K)),4])
      control_data <- knn(train, test, kyphosis[-(((i-1)*K+1):(i*K)),1], k = neighbours)
      for(j in 1:K){
         if(control_data[j] != kyphosis$Kyphosis[(i-1)*K+j]){
            tot = tot + 1
         }
      }
   }
   if(rest != 0){
      test <- cbind(kyphosis[(T*K+1):(T*K+rest),2], kyphosis[(T*K+1):(T*K+rest),3], kyphosis[(T*K+1):(T*K+rest),4])
      train <- cbind(kyphosis[-((T*K+1):(T*K+rest)),2], kyphosis[-((T*K+1):(T*K+rest)),3], kyphosis[-((T*K+1):(T*K+rest)),4])
      control_data <- knn(train, test, kyphosis[-((T*K+1):(T*K+rest)),1], k = neighbours)
      for(i in 1:rest){
         if(control_data[i] != kyphosis$Kyphosis[T*K+i]){
            tot = tot + 1
         }             
      }
   }    
   tot = tot / n
   tot
}

ten_fold_cv <- rep(0,20)
for(i in 1:20){
   ten_fold_cv[i] <- Kfold_CV(10, i)
}
plot(ten_fold_cv, type = "o", col = "green", lwd = 2,
     ylim = c(0,0.3), xlab = "Number of neighbors", ylab = "Misclassification Error")
legend(14, 0.05, col ="green", lty = 1,lwd=2, "10-fold CV")

