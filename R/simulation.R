normF <- function(mat){
  return(sum(mat^2))
}#normF


#the output is non-negative
rprops <- function(nsamples = 100, k = 10, M = 3){
  
  #hyper-parameter
  rho <- 2
  phi <- 0.1
  a_sigma <- 10
  b_sigma <- 1
  S <- matrix(0, nrow = k, ncol = nsamples)
  
  #generate mean, variance
  mean_kM <- matrix(rnorm(k*M,mean = rho,sd = 1/sqrt(phi)), nrow = k, ncol = M)
  var_kM <- matrix(rgamma(k*M,shape = a_sigma, scale = b_sigma), nrow = k, ncol = M)
  
  #seperate mean_kM
  for (i in 1:k){
    mean_kM[i,1] <- min(mean_kM[i,])
    for (j in 2:M){
      mean_kM[i,j] <- mean_kM[i,j-1]+runif(1,min=2,max=4)
    }#for j
  }#for i
  
  #generate proportion and indication variable
  prop_kM <- matrix(sample(10,k*M,replace = T),nrow = k, ncol = M)
  #for (i in 2:M){
  #prop_kM[,i] <- prop_kM[,i-1]+10
  #prop_kM[,i] <- prop_kM[,i-1]+0
  #}#for i
  prop_kM <- prop_kM/apply(prop_kM,1,sum)
  #categorical based on the prop
  ind_kM <- c()
  for (i in 1:k){
    ind_kM <- rbind(ind_kM, sample(1:M, nsamples, replace = T, prob = prop_kM[i,]))
  }#for i
  
  #generate samples and assign to each entries
  for (i in 1:k){
    for (m in 1:M){
      #see how many m are there
      m_index <- which(ind_kM[i,]==m)
      m_count <- length(m_index)
      #generate samples and assign
      m_samples <- rnorm(m_count,mean = mean_kM[i,m],sd = 1/sqrt(var_kM[i,m]))
      S[i,m_index] <- m_samples
    }#for m
  }#for i
  
  #single norm check
  #S <- matrix(rnorm(k*nsamples), nrow = k, ncol = nsamples)
  #uniform check
  S <- matrix(runif(k*nsamples,min = -1,max=1), nrow = k, ncol = nsamples)
  
  #return all information
  return(list(S = S, mean_kM = mean_kM, var_kM = var_kM, prop_kM = prop_kM, ind_kM = ind_kM))
}#rprops


#the output is a gene by cell.type matrix, all non-negative entries
rcoef <- function(ngenes, k, shape, scale, colinear.sd = 3){
  W <- matrix(nrow = ngenes, ncol = k)
  
  W[,1] <- rgamma(ngenes, shape = shape, scale = scale)
  
  for (i in 2:k){
    W[,i] <- W[,1]+abs(rnorm(ngenes, mean = 0, sd = colinear.sd))
  }#for i
  return(W)
}#rcoef



simulateData <- function(ngenes = 300, nsamples = 10, k = 5, M = 3, sd = 1, colinear.sd = 3, shape = 5, scale = 1, PC.frac = 1, threshold = seq(0.01, 0.1, by = 0.01)){
  S_list <- rprops(nsamples, k, M)
  S <- S_list$S
  A <- rcoef(ngenes, k, shape = shape, scale = scale, colinear.sd = colinear.sd)
  
  #create I from A, pick up genes based on a random threshold
  I <- matrix(0,nrow = ngenes, ncol = k)
  
  for (i in 1: (ncol(I)*PC.frac)){
    frac <-sample(threshold,1)
    index <- order(A[,i],decreasing = T)[1:(ngenes*frac)]
    I[index,i] <- 1
    #A[setdiff(1:nrow(A),index),i] <- abs(rnorm(length(setdiff(1:nrow(A),index))))
    A[setdiff(1:nrow(A),index),i] <- 0
  }#for i
  
  #noise needs to be non-negative
  X <- A%*%S
  noise <- matrix(rnorm(prod(dim(X)), 0, sd), nrow = ngenes)
  print(normF(X)/normF(noise))
  X <- X+noise
  
  #add some additional noise pathways
  #currently no
  
  return(list(X = X, A = A, S = S, I = I,mean_kM = S_list$mean_kM, var_kM = S_list$var_kM, prop_kM = S_list$prop_kM, ind_kM = S_list$ind_kM))
}#simulateData






