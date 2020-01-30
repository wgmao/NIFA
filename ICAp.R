tscale=function(mat){
  t(scale(t(mat)))
}

normF <- function(x){
  return(sum(x^2))
}#normF


fstats=function (dat, mod) 
{
  mod0=cbind(rep(1,ncol(dat)))
  n = dim(dat)[2]
  m = dim(dat)[1]
  df1 = dim(mod)[2]
  df0 = dim(mod0)[2]
  p = rep(0, m)
  Id = diag(n)
  resid = dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% 
                     t(mod))
  resid0 = dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*% 
                      t(mod0))
  rss1 = resid^2 %*% rep(1, n)
  rss0 = resid0^2 %*% rep(1, n)
  fstats = ((rss0 - rss1)/(df1 - df0))/(rss1/(n - df1))
  return(fstats)
} #fstats

one_hot_encode <- function(x){
  nrow <- length(x)
  cat <- names(table(x))
  ncol <- length(cat)
  res <- matrix(0,nrow=nrow,ncol=ncol)
  for (i in 1:ncol){
    res[which(x==cat[i]),i] <- 1
  }#for i
  return(res)
}#one_hot_encode


ICAp <- function(X, K = 6, M = 4, init="sd", S.init=NULL, A.init=NULL, L1.sd = NULL, L2.sd = NULL , verbose=F, beta_threshold = 2e-5, S_threshold = 6e-5, max.iter = 1000, ref=NULL, rho_S_diff=1.5,phi_S_prior=1, a_S_prior=1,b_S_prior=1, a_noise_prior=1, b_noise_prior=1,lambda_A_prior=1, eta_A_prior=0, ELBO=F, A.immediate =F, mixture.inner = 1, beta_expect_flag = NULL){
  P <- nrow(X)
  N <- ncol(X)
  
  #initialize variational variables randomly
  if (init=="random"){#################################################################################################
    set.seed(1)
    sigma_S <- array(abs(rnorm(K*K*N)), dim = c(K, K, N))
    if (is.null(S.init)){
      mu_S <- matrix(runif(K*N),nrow = K, ncol = N)
    }else{
      mu_S <- S.init
    }#else
      
    #mu_A <- svdres$u
    #mu_S <- t(svdres$v)
    
    lambda_S <- array(runif(K*M*N), dim = c(K, M, N))
    for (i in 1:K){
      for (j in 1:N){
        lambda_S[i,,j] <- lambda_S[i,,j]/sum(lambda_S[i,,j])
      }#for j
    }#for i
    pi_S <- pi_S_update(c(K,M),lambda_S)
    rho_S <- matrix(-2, nrow = K, ncol = M)
    for (i in 2:M){
      rho_S[,i] <- rho_S[,i-1]+rho_S_diff
    }#for i
    
    phi_S <- matrix(phi_S_prior, nrow = K, ncol = M)
    a_S <- matrix(rep(a_S_prior,K*M), nrow = K, ncol = M) #positive
    b_S <- matrix(rep(b_S_prior,K*M), nrow = K, ncol = M) #positive
    
    a_noise <- a_noise_prior #positive
    b_noise <- b_noise_prior #positive
    lambda_A <- matrix(rep(lambda_A_prior,K), nrow = 1, ncol = K) 
    eta_A <- matrix(rep(eta_A_prior,K*P), nrow = P, ncol = K)  #positive
    
    
    mu_A <- mu_A_update(c(K,P),eta_A, lambda_A)
    sigma_A <- sigma_A_update(c(K,K,P), eta_A, lambda_A)
  }else if (init=="ICA"){#################################################################################################
    set.seed(1)
    ICA.res <- fastICA(X, n.comp = K)
    
    sigma_S <- array(rnorm(K*K*N), dim = c(K, K, N))
    lambda_S <- array(runif(K*M*N), dim = c(K, M, N))
    for (i in 1:K){
      for (j in 1:N){
        lambda_S[i,,j] <- lambda_S[i,,j]/sum(lambda_S[i,,j])
      }#for j
    }#for i
    rho_S <- matrix(rho_S_diff, nrow = K, ncol = M)
    phi_S <- matrix(phi_S_prior, nrow = K, ncol = M)
    pi_S <- pi_S_update(c(K,M),lambda_S)
    a_S <- matrix(abs(rnorm(K*M)), nrow = K, ncol = M) #positive
    b_S <- matrix(abs(rnorm(K*M)), nrow = K, ncol = M) #positive
    a_noise <- a_noise_prior #positive
    b_noise <- b_noise_prior #positive
    lambda_A <- matrix(rep(lambda_A_prior,K), nrow = 1, ncol = K) 
    eta_A <- matrix(rep(eta_A_prior,K*P), nrow = P, ncol = K)  #positive
    
    if (is.null(S.init)){
      mu_S <- ICA.res$A
    }else{
      mu_S <- S.init
    }#else
    
    if (is.null(A.init)){
      mu_A <- t(ICA.res$S)  
    }else{
      mu_A <- A.init
    }#else
    
    sigma_A <- sigma_A_update(c(K,K,P), eta_A, lambda_A)
  }else if (init=="sd"){#################################################################################################
    #simple decomposition
    set.seed(1)
    if ( (is.null(S.init))| (is.null(A.init))){
      sd.res <- simpleDecomp(Y=X,k=K,L1=L1.sd,L2= L2.sd,max.iter=50,Zpos = T)
    }#if
    
    
    if (is.null(S.init)){
      mu_S <- sd.res$B
    }else{
      mu_S <- S.init
    }#else
    
    if (is.null(A.init)){
      mu_A <- t(sd.res$Z)  
    }else{
      mu_A <- A.init
    }#else
    
    
    a_S <- matrix(abs(rnorm(K*M)), nrow = K, ncol = M) #positive
    b_S <- matrix(abs(rnorm(K*M)), nrow = K, ncol = M) #positive
    a_noise <- a_noise_prior #positive
    b_noise <- b_noise_prior #positive
    
    sigma_S <- array(rnorm(K*K*N), dim = c(K, K, N))
    lambda_S <- array(runif(K*M*N), dim = c(K, M, N))
    rho_S <- matrix(rho_S_diff, nrow = K, ncol = M)
    phi_S <- matrix(phi_S_prior, nrow = K, ncol = M)
    
    
    for (i in 1:K){
      if (sd(mu_S[i,])==0){
        for (j in 1:N){
          lambda_S[i,,j] <- lambda_S[i,,j]/sum(lambda_S[i,,j])
        }#for j
        rho_S[i,] <- runif(M, min = min(mu_S[i,]) , max=max(mu_S[i,]))
      }else{
        gm.res <- Mclust(mu_S[i,], G=M)
        lambda_S[i,,] <- t(gm.res$z)
        rho_S[i,] <- gm.res$parameters$mean 
      }#else
    }#for i
    
    
    #simple decomposition with random means
    ###!!!
    #for (i in 1:K){
    #  for (j in 1:N){
    #    lambda_S[i,,j] <- lambda_S[i,,j]/sum(lambda_S[i,,j])
    #  }#for j
    #}#for i
    
    #for (i in 1:K){
    #  rho_S[i,] <- runif(M, min = min(sd.res$B[i,]) , max=max(sd.res$B[i,]))
    #}#for i
    ###!!!

    pi_S <- pi_S_update(c(K,M),lambda_S)
    lambda_A <- matrix(rep(lambda_A_prior,K), nrow = 1, ncol = K) 
    eta_A <- matrix(t(mu_A), nrow = P, ncol = K)  #positive
    
    sigma_A <- sigma_A_update(c(K,K,P), eta_A, lambda_A)
  }else{#################################################################################################
    #simple decomposition hybrid
    set.seed(1)
    print(K)
    sd.res <- simpleDecompHybrid(Y=X,k=K, k.max = 50, L1=10, L2=10, max.iter=50, Zpos = T)
    
    if (is.null(S.init)){
      mu_S <- sd.res$B
    }else{
      mu_S <- S.init
    }#else
    
    mu_A <- t(sd.res$Z)
    
    a_S <- matrix(abs(rnorm(K*M)), nrow = K, ncol = M) #positive
    b_S <- matrix(abs(rnorm(K*M)), nrow = K, ncol = M) #positive
    a_noise <- a_noise_prior #positive
    b_noise <- b_noise_prior #positive
    
    sigma_S <- array(rnorm(K*K*N), dim = c(K, K, N))
    lambda_S <- array(runif(K*M*N), dim = c(K, M, N))
    rho_S <- matrix(rho_S_diff, nrow = K, ncol = M)
    phi_S <- matrix(phi_S_prior, nrow = K, ncol = M)
    
    
    for (i in 1:K){
      if (length(table(sd.res$B[i,])) <= M){
        count <- length(table(sd.res$B[i,]))
        lambda_S[i,,] <- rbind(t(one_hot_encode(sd.res$B[i,])), matrix(0, nrow = M-count, ncol =ncol(sd.res$B) ))
        rho_S[i,] <- c(table(sd.res$B[i,]),rep(0,M-count))
      }else{
        gm.res <- Mclust(sd.res$B[i,], G=M)
        lambda_S[i,,] <- t(gm.res$z)
        rho_S[i,] <- gm.res$parameters$mean
      }#else
    }#for i
  
    
    pi_S <- pi_S_update(c(K,M),lambda_S)
    lambda_A <- matrix(rep(lambda_A_prior,K), nrow = 1, ncol = K) 
    eta_A <- matrix(sd.res$Z, nrow = P, ncol = K)  #positive
    
    sigma_A <- sigma_A_update(c(K,K,P), eta_A, lambda_A)
  }#else
  
  
  #record the initial values
  rho_S_0 <- rho_S
  phi_S_0 <- phi_S
  a_S_0 <- a_S
  b_S_0 <- b_S
  lambda_S_0 <- lambda_S
  pi_S_0 <- pi_S
  a_noise_0 <- a_noise
  b_noise_0 <- b_noise
  lambda_A_0 <- lambda_A
  eta_A_0 <- eta_A
  
  
  #########################################################################################
  #loop
  mu_S_accu <- c()
  if (ELBO){
    ELBO_accu <- c()  
  }else{
    ELBO_accu <- NULL
  }#else
  
  
  if (is.null(beta_expect_flag)){
    beta_expect_old <- 1e-2  
  }else{
    beta_expect <- beta_expect_flag
  }#else
  
  mu_S_old <- matrix(0,nrow = K, ncol = N)
  loop <- 0
  
  while (  (normF(mu_S-mu_S_old)/normF(mu_S_old) >=S_threshold) & (loop <= max.iter)){
    loop <- loop+1
    if (verbose){
      print(loop)
      print(paste("mu_S relative change:",normF(mu_S-mu_S_old)/normF(mu_S_old)))
      
      mu_S_accu <- c(mu_S_accu, normF(mu_S-mu_S_old)/normF(mu_S_old))
      if (length(mu_S_accu)>2){
        mu_S_accu.len <- length(mu_S_accu)
        #print( (mu_S_accu[mu_S_accu.len]-mu_S_accu[mu_S_accu.len-1])/mu_S_accu[mu_S_accu.len-1])
        second.der <- (mu_S_accu[mu_S_accu.len]-mu_S_accu[mu_S_accu.len-1])/mu_S_accu[mu_S_accu.len-1]
        if (abs(second.der) < 1e-2){
          break
        }#
      }#if 
    }#if verbose
    
    #update moment
    S_expect <- mu_S 
    S_2_expect <- S_2_expect_update(mu_S, sigma_S)
    mu_S_expect <- rho_S
    epsilon_expect <- lambda_S
    sigma_log_S_expect <- digamma(a_S)-log(b_S)
    
    if (is.null(beta_expect_flag)){
      beta_expect <- a_noise/b_noise  
    }#if
    
    sigma_expect <- a_S/b_S
    mu_S_2_expect <- 1/(phi_S*sigma_expect)+rho_S^2
    
    mean_A_expect <- mean_A_expect_update(dim(mu_A), mu_A)
    var_A_expect <- var_A_expect_update(mu_A, sigma_A)
    
    
    #simple decomposition only
    #!!! break
    
    #update variational variables using moments
    sigma_S <- sigma_S_update(dim(sigma_S), epsilon_expect,sigma_expect, beta_expect, var_A_expect) 
    mu_S_old <- mu_S
    #mu_S <- mu_S_update(dim(mu_S), sigma_S, X, mu_S_expect, epsilon_expect, sigma_expect, beta_expect, mean_A_expect)
    mu_S <- mu_S_update_fast(dim(mu_S), sigma_S, X, mu_S_expect, epsilon_expect, sigma_expect, beta_expect, mean_A_expect)
    
    S_expect <- mu_S ##@
    S_2_expect <- S_2_expect_update(mu_S, sigma_S) ##@
    
    
    phi_S_tmp <- phi_S_0
    rho_S_tmp <- rho_S_0
    phi_S <- phi_S_update(dim(phi_S), phi_S_tmp, epsilon_expect)
    rho_S <- rho_S_update(dim(rho_S), phi_S, phi_S_tmp, rho_S_tmp, sigma_expect, epsilon_expect, S_expect)
    
    mu_S_expect <- rho_S
    mu_S_2_expect <- 1/(phi_S*sigma_expect)+rho_S^2
    
    
    lambda_A <- lambda_A_update(dim(lambda_A), lambda_A_0, beta_expect,S_2_expect)
    
    if (A.immediate){
      #eta_A <- eta_A_update_cycle(dim(eta_A), lambda_A_0, eta_A_0, lambda_A, beta_expect, S_expect, mean_A_expect, X)
      eta_A <- eta_A_update_cycle_fast(dim(eta_A), lambda_A_0, eta_A_0, lambda_A, beta_expect, S_expect, mean_A_expect, X)  
    }else{
      #eta_A <- eta_A_update(dim(eta_A), lambda_A_0, eta_A_0, lambda_A, beta_expect, S_expect, mean_A_expect, X)
      eta_A <- eta_A_update_fast(dim(eta_A), lambda_A_0, eta_A_0, lambda_A, beta_expect, S_expect, mean_A_expect, X)  
    }#else
    

    sigma_A <- sigma_A_update(dim(sigma_A), eta_A, lambda_A)
    mu_A <- mu_A_update(dim(mu_A), eta_A, lambda_A)
    
    mean_A_expect <- mean_A_expect_update(dim(mu_A), mu_A)
    var_A_expect <- var_A_expect_update(mu_A, sigma_A)
    
    
    #membership
    #lambda_S <- lambda_S_update(dim(lambda_S),sigma_log_S_expect, pi_S,sigma_expect,S_expect,S_2_expect, mu_S_expect, mu_S_2_expect)
    lambda_S <- lambda_S_update_fast(dim(lambda_S),sigma_log_S_expect, pi_S,sigma_expect,S_expect,S_2_expect, mu_S_expect, mu_S_2_expect)
    
    epsilon_expect <- lambda_S ##@
    pi_S_old <- pi_S
    pi_S <- pi_S_update(dim(pi_S),lambda_S)
    if (verbose){
      print(paste("pi_S max change:",max(abs(pi_S-pi_S_old))))  
    }#verbose
    
    
    a_S <- a_S_update(dim(a_S),a_S_0,epsilon_expect)
    b_S <- b_S_update(dim(b_S),b_S_0, epsilon_expect, S_expect, S_2_expect, mu_S_expect, mu_S_2_expect, phi_S, rho_S)
    
    if (is.null(beta_expect_flag)){
      if (verbose){
        print(paste("beta_expect relative change:" ,abs( (beta_expect-beta_expect_old)/beta_expect_old )))
        print(beta_expect)
      }#if
      
      if (abs( (beta_expect-beta_expect_old)/beta_expect_old ) >= beta_threshold){
        a_noise <- a_noise_update(a_noise_0,N,P)
        #b_noise <- b_noise_update(b_noise_0, X, mean_A_expect, var_A_expect, S_expect)
        b_noise <- b_noise_update_fast(b_noise_0, X, mean_A_expect, var_A_expect, S_2_expect, S_expect)
        beta_expect_old <- beta_expect
        a_noise <- max(a_noise, 1e-10) #positive
        b_noise <- max(b_noise, 1e-10) #positive
      }#if  
    }else{
      #a_noise <- a_noise_update(a_noise_0,N,P)
      #b_noise <- b_noise_update_fast(b_noise_0, X, mean_A_expect, var_A_expect, S_2_expect, S_expect)
      #print(b_noise)
      #print(normF(X-t(mean_A_expect)%*%mu_S))
      #print(a_noise/normF(X-t(mean_A_expect)%*%mu_S))
    }#else
    
    
    
    
    #positive constratint
    a_S <- pmax(a_S,1e-10) #positive
    b_S <- pmax(b_S,1e-10) #positive
    
    
    if (verbose & (!is.null(ref))){
      print( apply(abs(cor(t(mu_S),ref)),2,max))
      print( apply(abs(cor(t(mu_S),ref)),2, which.max))
    }#if
    
    if (ELBO){
      ELBO.res <- list(a_noise=a_noise, b_noise=b_noise, beta_expect = beta_expect, mean_A_expect = mean_A_expect, S_expect= S_expect, var_A_expect = var_A_expect, epsilon_expect = epsilon_expect, sigma_log_S_expect = sigma_log_S_expect, sigma_expect = sigma_expect, S_2_expect = S_2_expect, mu_S_expect= mu_S_expect, mu_S_2_expect= mu_S_2_expect,sigma_S=sigma_S, lambda_A=lambda_A, eta_A=eta_A, lambda_S=lambda_S,phi_S=phi_S, rho_S=rho_S, a_S=a_S, b_S=b_S)
      
      ELBO.initialize.res <- list(X=X, lambda_A_0=lambda_A_0, eta_A_0=eta_A_0, pi_S_0=pi_S_0,phi_S_0=phi_S_0,rho_S_0=rho_S_0, a_S_0=a_S_0,b_S_0=b_S_0, a_noise_0=a_noise_0,b_noise_0=b_noise_0)
      ELBO_accu <- c(ELBO_accu, ELBO_cal(ELBO.res, ELBO.initialize.res))  
    }#ELBO

    
}#while

 
  if (is.null(beta_expect_flag)){
    beta_change <- abs( (beta_expect-beta_expect_old)/beta_expect_old )  
  }else{
    beta_change <- NULL
  }#else
  

return(list(X = X, K= K, M=M, beta_threshold = beta_threshold, S_threshold = S_threshold, rho_S_diff=rho_S_diff,phi_S_prior=phi_S_prior, a_S_prior=a_S_prior,b_S_prior=b_S_prior, a_noise_prior=a_noise_prior, b_noise_prior=b_noise_prior,lambda_A_prior=lambda_A_prior, eta_A_prior=eta_A_prior, loop = loop, S_expect = S_expect, S_2_expect = S_2_expect, mu_S_expect = mu_S_expect, epsilon_expect = epsilon_expect, sigma_log_S_expect = sigma_log_S_expect, beta_expect = beta_expect, sigma_expect = sigma_expect, mu_S_2_expect = mu_S_2_expect, mean_A_expect = mean_A_expect, var_A_expect = var_A_expect, sigma_S = sigma_S, mu_S = mu_S, lambda_S = lambda_S, pi_S = pi_S, rho_S = rho_S, phi_S = phi_S, a_S = a_S, b_S = b_S, a_noise = a_noise, b_noise = b_noise, lambda_A = lambda_A, eta_A = eta_A, mu_A = mu_A, sigma_A = sigma_A, rho_S_0 = rho_S_0, phi_S_0 = phi_S_0, a_S_0 = a_S_0, b_S_0 = b_S_0, lambda_S_0 = lambda_S_0, pi_S_0 = pi_S_0, a_noise_0 = a_noise_0, b_noise_0 = b_noise_0, lambda_A_0 = lambda_A_0, eta_A_0 = eta_A_0, mu_S_change = normF(mu_S-mu_S_old)/normF(mu_S_old), beta_change = beta_change, ELBO_accu = ELBO_accu))
}#ICAp





ELBO_cal <- function(res, initialize.res){
  accu <- c()
  
  a_noise <- res$a_noise
  b_noise <- res$b_noise
  beta_expect <- res$beta_expect
  mean_A_expect <- res$mean_A_expect
  S_expect <- res$S_expect
  var_A_expect <- res$var_A_expect
  epsilon_expect <- res$epsilon_expect
  sigma_log_S_expect <- res$sigma_log_S_expect
  sigma_expect <- res$sigma_expect
  S_2_expect <- res$S_2_expect
  mu_S_expect <- res$mu_S_expect
  mu_S_2_expect <- res$mu_S_2_expect
  sigma_S <- res$sigma_S
  lambda_A <- res$lambda_A
  eta_A <- res$eta_A
  lambda_S <- res$lambda_S
  phi_S <- res$phi_S
  rho_S <- res$rho_S
  a_S <- res$a_S
  b_S <- res$b_S
  
  beta_log_expect <- digamma(a_noise)-log(b_noise)
  
  X <- initialize.res$X
  lambda_A_0 <- initialize.res$lambda_A_0
  eta_A_0 <- initialize.res$eta_A_0
  pi_S_0 <- initialize.res$pi_S_0
  phi_S_0 <- initialize.res$phi_S_0
  rho_S_0 <- initialize.res$rho_S_0
  a_S_0 <- initialize.res$a_S_0
  b_S_0 <- initialize.res$b_S_0
  a_noise_0 <- initialize.res$a_noise_0
  b_noise_0 <- initialize.res$b_noise_0
  
  
  
  accu <- c(accu, lower_bound_1(beta_log_expect, beta_expect, X, mean_A_expect, S_expect, var_A_expect))
  accu <- c(accu, lower_bound_2(epsilon_expect, sigma_log_S_expect, sigma_expect,S_2_expect, S_expect, mu_S_expect, mu_S_2_expect))
  accu <- c(accu, lower_bound_3(lambda_A_0, var_A_expect, eta_A_0))
  accu <- c(accu, lower_bound_4(epsilon_expect,pi_S_0))
  accu <- c(accu, lower_bound_5(phi_S_0, sigma_expect, sigma_log_S_expect, mu_S_2_expect, mu_S_expect, rho_S_0))
  accu <- c(accu, lower_bound_6(a_S_0, b_S_0, sigma_log_S_expect, sigma_expect))
  accu <- c(accu, lower_bound_7(a_noise_0, b_noise_0, beta_expect, beta_log_expect))
  accu <- c(accu, -lower_bound_8(sigma_S))
  accu <- c(accu, -lower_bound_9(lambda_A, var_A_expect, eta_A))
  accu <- c(accu, -lower_bound_10(lambda_S))
  accu <- c(accu, -lower_bound_11(phi_S, sigma_expect, sigma_log_S_expect, mu_S_2_expect, rho_S))
  #accu <- c(accu, -lower_bound_12(a_S, b_S, sigma_log_S_expect, sigma_expect))
  #accu <- c(accu, -lower_bound_13(a_noise, b_noise, beta_expect, beta_log_expect))
  
  return(sum(accu))
}#ELBO_cal

