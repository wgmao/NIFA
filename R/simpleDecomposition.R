#require(rsvd)
require(mclust)
require(MASS)


#Y is gene by sample (zscores across rows)
#svd will be computed if not supplied
simpleDecomp=function(Y, k, svdres=NULL, L1=NULL, L2=NULL, Zpos=T,max.iter=200, tol=5e-6, trace=F){
  ng=nrow(Y)
  ns=ncol(Y)
  
  Bdiff=-1
  BdiffTrace=double()
  BdiffCount=0
  if(is.null(svdres)){
    message("Computing SVD")
    #svdres=rsvd(Y, k = k)
    svdres <- svd(Y, nu = k, nv = k)
  }
  
  if (is.null(L2)){
    L2=svdres$d[k]
  }#if
  print(paste0("L2 is set to ",L2))
  
  if (is.null(L1)){
    L1=L2/2
  }
  print(paste0("L1 is set to ",L1))

  #initialize B with svd
  if (k>1){
    B=t(svdres$v[1:ncol(Y), 1:k]%*%diag((svdres$d[1:k])))  
  }else{
    B=matrix(t(svdres$v[1:ncol(Y), 1:k]), nrow=1)*svdres$d[1:k]
  }#else
  
  
  round2=function(x){signif(x,4)}
  for ( i in 1:max.iter){
    #main loop    
    Zraw=Z=(Y%*%t(B))%*%ginv(tcrossprod(B)+L1*diag(k))
    
    if(Zpos){
      Z[Z<0]=0
    }
    oldB=B
   # B=solve(crossprod(Z)+L2*diag(k))%*%t(Z)%*%Y
    B=ginv(t(Z)%*%Z+L2*diag(k))%*%t(Z)%*%Y
    
    
    #update error
    Bdiff=sum((B-oldB)^2)/sum(B^2)
    BdiffTrace=c(BdiffTrace, Bdiff)
    err0=sum((Y-Z%*%B)^2)+sum((Z)^2)*L1+sum(B^2)*L2
    if(trace){
    message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))))
    }
    
    #check for convergence
    if(i>52&&Bdiff>BdiffTrace[i-50]){
      BdiffCount=BdiffCount+1
    }
    else if(BdiffCount>1){
      BdiffCount=BdiffCount-1
    }
    
    if(Bdiff<tol){
      message(paste0("converged at  iteration ", i))
      break
    }
    if( BdiffCount>5){
      message(paste0("stopped at  iteration ", i, " Bdiff is not decreasing"))
      break
    }
  }
  return(list(B=B, Z=Z, Zraw=Zraw, Z2=ginv(tcrossprod(B)+L1*diag(k))))
}#simpleDecomp





simpleDecompGMM=function(Y, k,svdres=NULL, L1=1, L2=1, Zpos=T,max.iter=200, tol=5e-6, trace=F, rseed=NULL, gm=F, B=NULL, max.ggm.iter=10,  var.ind=F){
  ng=nrow(Y)
  ns=ncol(Y)
  
  Bdiff=-1
  BdiffTrace=double()
  BdiffCount=0
  
  if(is.null(svdres)){
    message("Computing SVD")
    #svdres=rsvd(Y, k = k) 
    svdres <- svd(Y, nu = k, nv = k)
  }
  
  if(is.null(B)){
    #initialize B with svd
    if (k>1){
      B=t(svdres$v[1:ncol(Y), 1:k]%*%diag((svdres$d[1:k])))  
    }else{
      B=matrix(t(svdres$v[1:ncol(Y), 1:k]), nrow=1)*svdres$d[1:k]
    }#else
  }
  else{
    message("B given")
  }
  Bmean=B
  if (!is.null(rseed)) {
    message("using random start")
    set.seed(rseed)
    B = t(apply(B, 1, sample))
  }
  round2=function(x){signif(x,4)}
  
  #variables to hold the mixture params
  mgL=list(k)
  varinv=double(k)
  allM=double(k)
  
  for ( i in 1:max.iter){
    #main loop    
    Zraw=Z=(Y%*%t(B))%*%solve(tcrossprod(B)+L1*diag(k))
    
    if(Zpos){
      Z[Z<0]=0
    }
    oldB=B
    # B=solve(crossprod(Z)+L2*diag(k))%*%t(Z)%*%Y
    
    #get mixture params for B
    if(gm){
      Bmean[]=0
      
      allM[]=1
      for(j in 1:k){
        
        mgL[[j]]=res=Mclust(data=B[j,], modelNames = c("V"), G = 1:6, control=emControl(itmax=max.ggm.iter), verbose = F, prior=priorControl())
        #get a pooled variance, not  actually used if var.ind=F
        varpool=sum(res$parameters$pro*res$parameters$variance$sigmasq)
        varinv[j]=1/varpool
        allM[j]=length(res$parameters$mean)
        
        Bmean[j,]=rowSums(sweep(res$z, 2, res$parameters$pro*(varinv[j])*res$parameters$mean, "*"))
      }#for j
      #the variance is a single parameter folded into L2
      if(!var.ind){
        B=solve(t(Z)%*%Z+L2*diag(k))%*%(t(Z)%*%Y+L2*Bmean)
      }else{ #use the varinv, may not converge
        B=solve(t(Z)%*%Z+L2*diag(varinv))%*%(t(Z)%*%Y+L2*sweep(Bmean,1,varinv, "*"))
      }#else
      
    }else{ 
      B=solve(t(Z)%*%Z+L2*diag(k))%*%(t(Z)%*%Y)
    }#no gmm
    
    Bsum=rowMeans(B^2)
    #deal with factors that drop out
    ii=which(Bsum<1e-3)
    if(length(ii)>0){
      k=k-length(ii)
      message(paste("new k=",k))
      
      oldB=oldB[-ii,]
      B=B[-ii,]
      Bmean=Bmean[-ii,]
      Z=Z[,-ii]
      allM=list()
    }#if
    #update error
    Bdiff=sum((B-oldB)^2)/sum(B^2)
    show(Bdiff)
    BdiffTrace=c(BdiffTrace, Bdiff)
    err0=sum((Y-Z%*%B)^2)+sum((Z)^2)*L1+sum(B^2)*L2
    if(trace){
      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))))
    }
    
    #check for convergence
    if(i>52&&Bdiff>BdiffTrace[i-50]){
      BdiffCount=BdiffCount+1
    }else if(BdiffCount>1){
      BdiffCount=BdiffCount-1
    }#else if
    
    if(Bdiff<tol){
      message(paste0("converged at  iteration ", i))
      break
    }#if (Bdiff<tol)
    
    if( BdiffCount>5){
      message(paste0("stopped at  iteration ", i, " Bdiff is not decreasing"))
      break
    }#if( BdiffCount>5)
  }#for i
  return(list(B=B, Z=Z, Zraw=Zraw, Z2=solve(tcrossprod(B)+L1*diag(k)), allM=allM, mgL=mgL, L1=L1, L2=L2))
}#simpleDecompGMM




simpleDecompHybrid <- function(Y, k, k.max, L1=10,L2=10,max.iter=50,Zpos = T){
  #SVD, first 50 components, kmeans clustering
  svd.dim <- 20
  svdres <- svd(Y, nu = svd.dim, nv = svd.dim)
  
  tol.dis <- c()
  tol.median <- c()
  for (i in 1:k.max){
    kmeans.res <- kmeans(svdres$v, centers = i, nstart = 10, iter.max = 20)
    tol.dis <- c(tol.dis, kmeans.res$tot.withinss)
    tol.median <- c(tol.median, median(table(kmeans.res$cluster)))
  }#for i
  
  k.kmeans <- quick.elbow(tol.dis[1:max(which(tol.median>1))])
  k.kmeans <- min(k,k.kmeans)
  kmeans.res <- kmeans(svdres$v, centers = k.kmeans, nstart = 10, iter.max = 20)
  
  #regress out
  lm.res <- lm(t(svdres$v)~t(kmeans.res$centers)+0)
  Yp <- svdres$u%*%diag(svdres$d[1:svd.dim])%*%lm.res$residuals
  

  #simple decomposition on the leftover
  if (k > k.kmeans){
    sd.res <- simpleDecomp(Y=Yp,k=k-k.kmeans,L1=L1,L2=L2,max.iter=max.iter, Zpos = Zpos)  
    return(list(B=rbind(sd.res$B, t(one_hot_encode(kmeans.res$cluster))), Z= cbind(sd.res$Z, abs(svdres$u[,1:k.kmeans] )), k.kmeans = k.kmeans))
  }else{
    return(list(B= t(one_hot_encode(kmeans.res$cluster))[1:k,], Z = abs(svdres$u[,1:k]), k.kmeans = k.kmeans))
  }#else
}#simpleDecompHybrid


simpleDecompHybrid.ridge <- function(Y, k, k.max, L1=10,L2=10,max.iter=50,Zpos = T){
  #SVD, first 50 components, kmeans clustering
  svd.dim <- 20
  svdres <- svd(Y, nu = svd.dim, nv = svd.dim)
  
  tol.dis <- c()
  tol.median <- c()
  for (i in 1:k.max){
    kmeans.res <- kmeans(svdres$v, centers = i, nstart = 10, iter.max = 20)
    tol.dis <- c(tol.dis, kmeans.res$tot.withinss)
    tol.median <- c(tol.median, median(table(kmeans.res$cluster)))
  }#for i
  
  k.kmeans <- quick.elbow(tol.dis[1:max(which(tol.median>1))])
  kmeans.res <- kmeans(svdres$v, centers = k.kmeans, nstart = 10, iter.max = 20)
  
  #regress out
  lm.res <- lm(t(svdres$v)~t(kmeans.res$centers)+0)
  Yp <- svdres$u%*%diag(svdres$d[1:svd.dim])%*%lm.res$residuals

  
  #simple decomposition on the leftover
  if (k > k.kmeans){
    sd.res <- simpleDecomp(Y=Yp,k=k-k.kmeans,L1=L1,L2=L2,max.iter=max.iter, Zpos = Zpos) 
    B <- rbind(sd.res$B, t(one_hot_encode(kmeans.res$cluster)))
    lm.res <- lm(t(Y)~t(B)+0)
    Z <- lm.res$coefficients
    return(list(B=B, Z= Z, k.kmeans = k.kmeans))
  }else{
    B <- t(one_hot_encode(kmeans.res$cluster))[1:k,]
    lm.res <- lm(t(Y)~t(B)+0)
    Z <- lm.res$coefficients
    return(list(B=B, Z = Z, k.kmeans = k.kmeans))
  }#else
}#simpleDecompHybrid.ridge





quick.elbow <- function(varpc,low=.08,max.pc=.9) {
  ee <- varpc/sum(varpc) # ensure sums to 1
  #print(round(log(ee),3))
  while(low>=max(ee)) { low <- low/2 } # when no big components, then adjust 'low'
  lowie <- (ee<low) ; highie <- ee>low/8
  low.ones <- which(lowie & highie)
  others <- length(which(!lowie))
  if(length(low.ones)>0) {
    if(length(low.ones)==1) {
      elbow <- low.ones 
    } else {
      set <- ee[low.ones]
      pc.drops <- abs(diff(set))/(set[1:(length(set)-1)])
      infz <- is.infinite(pc.drops)
      #print(pc.drops)
      elbow <- which(pc.drops==max(pc.drops[!infz],na.rm=T))[1]+others
    }
  } else { 
    # if somehow there are no small eigenvalues, just choose the elbow as the second last
    cat("no eigenvalues were significantly smaller than the previous\n")
    elbow <- length(ee) 
  }
  if(tail(cumsum(ee[1:elbow]),1)>max.pc) {
    elbow <- which(cumsum(ee)>max.pc)[1]-1
  }
  if(elbow<1) {
    warning("elbow calculation failed, return zero")
    return(0)
  }
  names(elbow) <- NULL
  return(elbow)
}#quick.elbow
