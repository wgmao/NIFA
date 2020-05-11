#K*K*N
S_2_expect_update <- function(mu_S, sigma_S){
  N <- dim(sigma_S)[3]
  res <- array(0, dim = dim(sigma_S))
  for (i in 1:N){
    res[,,i] <- mu_S[,i]%*%t(mu_S[,i])+sigma_S[,,i]
  }#for i
  return(res)
}#S_2_expect_update


####################################################################################################
lambda_A_update <- function(dim_lambda_A, lambda_A_0, beta_expect, S_2_expect){
  #dim_lambda_A = (1, K)
  res <- matrix(0, nrow = dim_lambda_A[1], ncol = dim_lambda_A[2])
  
  for (i in 1:dim_lambda_A[1]){
    for (j in 1:dim_lambda_A[2]){
      res[i,j] <- lambda_A_0[i,j]+beta_expect*sum(S_2_expect[j,j,])
    } #for j 
  }#for i
  return(res)
}#lambda_A_update


sigma_A_update <- function(dim_sigma_A, eta_A, lambda_A){
  #dim_sigma_A = (K,K,P)
  res <- array(0, dim = dim_sigma_A) 
  sigma <- sqrt(1/lambda_A)
  
  for (i in 1:dim_sigma_A[2]){
    for (j in 1:dim_sigma_A[3]){
      Z <- 1-pnorm(-eta_A[j,i]/sigma[i], mean = 0, sd = 1, lower.tail = T)
      if (Z==0){
        res[i,i,j] <- sigma[i]^2 
      }else{
        num <- dnorm(-eta_A[j,i]/sigma[i], mean = 0, sd = 1)
        res[i,i,j] <- sigma[i]^2*(1-eta_A[j,i]/sigma[i]*num/Z-(num/Z)^2)  
      }
    }#for j   
  }#for i
  
  return(res)
}#sigma_A_update


eta_A_update <- function(dim_eta_A, lambda_A_0, eta_A_0, lambda_A, beta_expect, S_expect, mean_A_expect, X){
  #dim_eta_A = (P,K)
  #lambda_A_0 = (1,K)
  #S_expect = (K, N)
  res <- matrix(0, nrow = dim_eta_A[1], ncol = dim_eta_A[2])
  
  P <- dim_eta_A[1]
  K <- dim_eta_A[2]
  N <- dim(S_expect)[2]
  
  prod <- t(mean_A_expect)%*%S_expect
  for (i in 1:P){
    for (j in 1:K){
      num <- 0
      out_prod <- mean_A_expect[j,i]*S_expect[j,]
        
      for (n in 1:N){
       num <- num+ (prod[i,n]-out_prod[n]-X[i,n])*S_expect[j,n]
    }#for n
       tmp <- lambda_A_0[1,j]*eta_A_0[i,j] - beta_expect*num
       res[i,j] <- tmp/lambda_A[1,j]
    }#for j
  }#for i
  
  return(res)
}#eta_A_update

  

#K*P
mu_A_update <- function(dim_mu_A, eta_A, lambda_A){
  #dim_mu_A = (K,P) !!!
  #eta_A = (P, K)
  res <- matrix(0,nrow=dim_mu_A[1],ncol=dim_mu_A[2])
  sigma <- sqrt(1/lambda_A)
  
  for (i in 1:dim_mu_A[1]){
    for (j in 1:dim_mu_A[2]){
      Z <- 1-pnorm(-eta_A[j,i]/sigma[i],mean = 0, sd = 1, lower.tail = T)
      if (Z==0){
        res[i,j] <- 0
      }else{
        num <- dnorm(-eta_A[j,i]/sigma[i],mean = 0, sd = 1)
        res[i,j] <- eta_A[j,i]+num/Z*sigma[i]  
      }#else
      
    }#for j
  }#for i

  return(res)
}#mu_A_update


eta_A_update_cycle <- function(dim_eta_A, lambda_A_0, eta_A_0, lambda_A, beta_expect, S_expect, mean_A_expect, X){
  #dim_eta_A = (P,K)
  #lambda_A_0 = (1,K)
  #S_expect = (K, N)
  res <- matrix(0, nrow = dim_eta_A[1], ncol = dim_eta_A[2])
  
  P <- dim_eta_A[1]
  K <- dim_eta_A[2]
  N <- dim(S_expect)[2]
  
  mean_A_expect_cycle <- mean_A_expect
  sigma <- sqrt(1/lambda_A)
  
  prod <- t(mean_A_expect_cycle)%*%S_expect
  for (i in 1:P){
    for (j in 1:K){
      num <- 0
      out_prod <- mean_A_expect_cycle[j,i]*S_expect[j,]
      
      for (n in 1:N){
        num <- num+ (prod[i,n]-out_prod[n]-X[i,n])*S_expect[j,n]
      }#for n
      tmp <- lambda_A_0[1,j]*eta_A_0[i,j] - beta_expect*num
      res[i,j] <- tmp/lambda_A[1,j]
      
      #update mu_A, K*P
      #######################################
      Z <- 1-pnorm(-res[i,j]/sigma[j],mean = 0, sd = 1, lower.tail = T)
      if (Z==0){
        mu_A_cycle <- 0
      }else{
        num <- dnorm(-res[i,j]/sigma[j],mean = 0, sd = 1)
        mu_A_cycle <- res[i,j]+num/Z*sigma[j]  
      }#else
        
      #update mean_A_expect_cycle
      #######################################
      mean_A_expect_cycle[j,i] <- mu_A_cycle
      
    }#for j
  }#for i
  
  return(res)
}#eta_A_update_cycle



mean_A_expect_update <- function(dim_mu_A, mu_A){
  return(mu_A)
}#mean_A_expect_update


var_A_expect_update <- function(mu_A, sigma_A){
  res <- sigma_A
  for (j in 1:dim(sigma_A)[3]){
    res[,,j] <- mu_A[,j]%*%t(mu_A[,j])+res[,,j]
  }#for i
  return(res)
}#var_A_expect_update


####################################################################################################
#sigma_expect: K*M
sigma_S_update <- function(dim_sigma_S, epsilon_expect,sigma_expect, beta_expect, var_A_expect){
  res <- array(0, dim=dim_sigma_S)
  M <- dim(sigma_expect)[2]
  left <- beta_expect*rowSums(var_A_expect,dims=2)
  for (n in 1:dim_sigma_S[3]){
    right <- 0
    for (m in 1:M){
      right <- right+diag(epsilon_expect[,m,n]*sigma_expect[,m])
    }#for m
    res[,,n] <- solve(left+right)
  }#for n
  return(res)
}#sigma_S_update


#K*N
mu_S_update <- function(dim_mu_S, sigma_S, X, mu_S_expect, epsilon_expect, sigma_expect, beta_expect, mean_A_expect){
  res <- matrix(0,nrow = dim_mu_S[1], ncol = dim_mu_S[2])
  P <- dim(X)[1]
  for (n in 1:dim_mu_S[2]){
    left <- X[1,n]*mean_A_expect[,1]
    for (j in 2:P){
      left <- left+X[j,n]*mean_A_expect[,j]
    }#for j
    left <- beta_expect*left
    
    right <- c()
    for (i in 1:dim_mu_S[1]){
      right <- c(right, sum(epsilon_expect[i,,n]*sigma_expect[i,]*mu_S_expect[i,]))
    }#for i
    right <- matrix(right, ncol=1)
    res[,n] <- sigma_S[,,n]%*%(left+right)
  }#for n
  return(res)
}#mu_S_update



phi_S_update <- function(dim_phi_S, phi_S, epsilon_expect){
  res <- matrix(0, nrow = dim_phi_S[1], ncol = dim_phi_S[2])
  for (i in 1:dim_phi_S[1]){
    for (m in 1:dim_phi_S[2]){
      res[i,m] <- phi_S[i,m]+sum(epsilon_expect[i,m,])
    }#for m
  }#for i
  return(res)
}#phi_S_update


rho_S_update <- function(dim_rho_S, phi_S, phi_S_tmp, rho_S_tmp, sigma_expect, epsilon_expect, S_expect){
  res <- matrix(0, nrow = dim_rho_S[1], ncol = dim_rho_S[2])
  for (i in 1:dim_rho_S[1]){
    for (m in 1:dim_rho_S[2]){
      res[i,m] <- 1/(phi_S[i,m]*sigma_expect[i,m])*( sigma_expect[i,m]*sum( epsilon_expect[i,m,]*S_expect[i,]) +phi_S_tmp[i,m]*rho_S_tmp[i,m]*sigma_expect[i,m])
    }#for m
  }#for i
  return(res)
}#rho_S_update



lambda_S_update <- function(dim_lambda_S,sigma_log_S_expect, pi_S,sigma_expect,S_expect,S_2_expect, mu_S_expect,mu_S_2_expect){
  res <- array(0, dim=dim_lambda_S)
  res_tmp <- array(0, dim=dim_lambda_S)
  for (i in 1:dim_lambda_S[1]){
    for (m in 1:dim_lambda_S[2]){
      for (n in 1:dim_lambda_S[3]){
        res_tmp[i,m,n] <- 0.5*sigma_log_S_expect[i,m]+log(pi_S[i,m])-0.5*log(2*pi)-0.5*sigma_expect[i,m]*(S_2_expect[i,i,n]-2*S_expect[i,n]*mu_S_expect[i,m]+mu_S_2_expect[i,m])
      }#for n
    }#for m
  }#for i
  
  for (i in 1:dim_lambda_S[1]){
    for (n in 1:dim_lambda_S[3]){
      sum_along_m <- sum(exp(res_tmp[i,,n]))
      if ((sum_along_m==0)|(sum_along_m==Inf)){
        for (m in 1:dim_lambda_S[2]){
          res[i,m,n] <- 1/sum(exp(res_tmp[i,,n]-res_tmp[i,m,n]))
        }#for m
      }else{
        for (m in 1:dim_lambda_S[2]){
          res[i,m,n] <- exp(res_tmp[i,m,n])/sum_along_m
        }#for m  
      }#else
    }#for n
  }#for i
  return(res)
}#lambda_S_update


pi_S_update <- function(dim_pi_S, lambda_S){
  res <- matrix(0, nrow = dim_pi_S[1], ncol = dim_pi_S[2])
  for (i in 1:dim_pi_S[1]){
    for (m in 1:dim_pi_S[2]){
      res[i,m] <- mean(lambda_S[i,m,])
    }#for m
  }#for i
  res <- pmax(res, 1e-10)
  return(res/rowSums(res))
}#pi_S



a_S_update <- function(dim_a_S,a_S,epsilon_expect){
  res <- matrix(0, nrow = dim_a_S[1], ncol = dim_a_S[2])
  for (i in 1:dim_a_S[1]){
    for (m in 1:dim_a_S[2]){
      res[i,m] <- a_S[i,m]+0.5*sum(epsilon_expect[i,m,])+0.5
    }#for m
  }#for i
  return(res)
}#a_S_update

#K*M
b_S_update <- function(dim_b_S,b_S, epsilon_expect, S_expect, S_2_expect, mu_S_expect, mu_S_2_expect, phi_S, rho_S){
  res <- matrix(0, nrow = dim_b_S[1], ncol = dim_b_S[2])
  for (i in 1:dim_b_S[1]){
    for (m in 1:dim_b_S[2]){
      res[i,m] <- b_S[i,m]+0.5*sum( epsilon_expect[i,m,]*(S_2_expect[i,i,]-2*S_expect[i,]*mu_S_expect[i,m]+mu_S_2_expect[i,m]))+0.5*phi_S[i,m]*( mu_S_2_expect[i,m]-2*rho_S[i,m]*mu_S_expect[i,m]+rho_S[i,m]^2 )
    }#for m
  }#for i
  return(res)
}#b_S_update



a_noise_update <- function(a_noise,N,P){
  return(a_noise+N*P/2)
}#a_noise_update


b_noise_update  <- function(b_noise, X, mean_A_expect, var_A_expect, S_2_expect, S_expect){
  res <- 0
  for (j in 1:nrow(X)){
    for (n in 1:ncol(X)){
      res <- res+0.5*(X[j,n]^2-2*X[j,n]*t(mean_A_expect[,j])%*%S_expect[,n]+sum(diag(var_A_expect[,,j]%*%S_2_expect[,,n])))
    }#for n
  }#for j
  return(b_noise+as.numeric(res))
}#b_noise_update



