GI <- function(Y, E, R, N, thin, burn_in, memory=list(nnzcolindices=6467)){
    n <- 544
    yes <- 0
    no <- 0
    alpha_u <- 1
    alpha_v <- 1
    beta_u <- 0.01
    beta_v <- 0.01
    Y = matrix(Y, ncol = 1)
    E = matrix(E, ncol = 1)
    #u <- rep(.2, n)
    #u = matrix(u, ncol = 1)
    #eta <- rep(1, n)
    #eta = matrix(eta, ncol = 1)
    
    Ku <- matrix(nrow = 1, ncol = N);
    Kv <- matrix(nrow = 1, ncol = N);
    u   <- matrix(data=1, nrow = n, ncol = N);
    eta <- matrix(data=1, nrow = n, ncol = N);
    Ku[,1]  <- 20;
    Kv[,1]  <- 400;
    u[,1]   <- 0;
    eta[,1] <- 1;
    
    for(i in 2:N){
       for(j in 1:thin){
           
      Ku[,i] <- rgamma(1, (n-2)/2+alpha_u, beta_u+t(u[,i-1])%*%R%*%u[,i-1]*0.5)
	    Kv[,i] <- rgamma(1, n/2+alpha_v, beta_v+t(eta[,i-1]-u[,i-1])%*%(eta[,i-1]-u[,i-1])*0.5)
	    u[,i] <- rmvnorm.canonical(1, Kv[,i]*eta[,i-1], Ku[,i]*R+diag(Kv[,i], n))
	    u[,i] <- t(u[,i])
	    b <- Y+E*exp(eta[,i-1])*(eta[,i-1]-1)
	    cv <- E*exp(eta[,i-1])
      new <- c(t(rmvnorm.canonical(1, Kv[,i]*u[,i]+b,  diag.spam(Kv[,i] + c(cv)), memory=list(nnzcolindices=6467))))
           # new <- t(new)
	    log.target.ratio <- {(-t(new)%*%new*Kv[,i]*0.5+t(new)%*%u[,i]*Kv[,i]+t(new)%*%Y-t(exp(new))%*%E)
				-(-t(eta[,i-1])%*%eta[,i-1]*Kv[,i]*0.5+t(eta[,i-1])%*%u[,i]*Kv[,i]+t(eta[,i-1])%*%Y-t(exp(eta[,i-1]))%*%E)}
	    b_new <- Y+E*exp(new)*(new-1)
	    c_new <- E*exp(new)
	    Q_num <- diag.spam(Kv[,i] + c(c_new)); #Kv*diag.spam(1,n)+diag(c_new)
	    Q_den <- diag.spam(Kv[,i] + c(cv)); #Kv*diag.spam(1,n)+diag(c)
	    Epsilon_num <- Kv[,i]*u[,i]+b_new
	    Epsilon_den <- Kv[,i]*u[,i]+b 

            log.proposal.ratio <- {dmvnorm.canonical(eta[,i-1], Epsilon_num, Q_num, log=TRUE)
                                  -dmvnorm.canonical(new, Epsilon_den, Q_den, log=TRUE)}
	   
	    alpha <- min(0, log.target.ratio+log.proposal.ratio)
      print(exp(alpha))  
			print(i);
	    if(runif(1) <= exp(alpha)){
	        print("ACCEPT");
	        eta[,i] <- new
	       	yes <- yes + 1	
	    }
	    else{
	        print("REJECT");
          eta[,i] <- eta[,i-1]
 	      	no <- no + 1
 	    }

        }
	acceptance <- yes/(yes + no) * 100
#	container[i,] <- c(Ku, Kv, u[387,i], eta[387,i], acceptance)
    }    
   # container
return(list(
  "kv"  = Kv[(burn_in+1):N],
  "ku"  = Ku[(burn_in+1):N],
  "u"   = u[,(burn_in+1):N],
  "eta" = eta[,(burn_in+1):N],
  "v"   = eta[,(burn_in+1):N] - u[,(burn_in+1):N],
  "accepted" = yes/N,
  "rejected" = no/N
));

}
t1 <- Sys.time()
GI(Y, E, R, 2000,1,1)->vol
print(Sys.time() - t1)
vol[1,]
