NMF_generate_param <- function(M,N,H,shape,scale){
  LIMIT <- 100
  for(i in 1:LIMIT){
    A <- matrix(rgamma(H*M,shape = shape,scale = scale),nrow=N,ncol=H)
    B <- matrix(rgamma(N*H,shape = shape,scale = scale),nrow=H,ncol=M)
    if( ( min(dim(A)) == qr(A)$rank ) && ( min(dim(B)) == qr(B)$rank ) ){
      return(list(A=A,B=B))
    }
  }
  print("missing to create true matrix until 100 loops")
  return(list())
}

rNMF <- function(n,true.param){
  N <- dim(true.param$A)[1]
  M <- dim(true.param$B)[2]
  x <- array(0,dim=c(N,M,n))
  for(i in 1:n){
    x[,,i] <- matrix(rpois(N*M,lambda=true.param$A%*%true.param$B),nrow=N,ncol=M)
  }
  return(x)
}

disp.param <- function(param){
  print(param$A%*%param$B)
}

NMF.lambda <- function(M,N,H,H0){
  m <- 1
  lambda <- ((H-H0)*min(c(M,N))+H0*(M+N-1))/2
  return(list(lambda=lambda,m=m))
}

nmf.poisson.em <- function(data, learning.seed=1, H, 
                   init.shape=8, iteration=100, seed.num = 50){
  N <- dim(data)[1]
  M <- dim(data)[2]
  n <- dim(data)[3]
  
  set.seed(learning.seed)

  local.mle <- -Inf
  est.param <- list()

  ##change initial parameters seed.num times
  for(seed.ind in 1:seed.num){
    ##initialize parameters
    A <- matrix(rgamma(N*H,shape = init.shape),nrow=N,ncol=H)
    B <- matrix(rgamma(H*M,shape = init.shape),nrow=H,ncol=M)
    
    s <- array(0,dim=c(N,H,M,n))
    
    ##EM algorithm
    for(ite in 1:iteration){
      #E-step
      for(nu in 1:N){
        for(tau in 1:M){
          L <- log(A[nu,])+log(B[,tau])
          L <- L - rep(1,H) * max(L)
          factor_ratio <- exp(L)/sum(exp(L))
          s[nu,,tau,] <- factor_ratio%*%t(data[nu,tau,])
        }
      }
      
      #M-step
      A <- apply(s,c(1,2),sum)/(n*rep(1,N)%*%t(apply(B,1,sum)))
      B <- apply(s,c(2,3),sum)/(n*apply(A,2,sum)%*%t(rep(1,M)))
      
      current.param <- list(A=A,B=B)
      current.log.like <- dnmflikelihood(data=x, param=current.param, log=T)
      
      if(is.nan(current.log.like)){
        break
      }
      
      previous.loglike <- current.log.like
      previous.param <- current.param
      
      # print(log.like)
    }
    #calc likelihood
    candidate.param <- current.param
    candidate.log.like <- current.log.like
      
    if(candidate.log.like>local.mle){
      print(candidate.log.like)
      local.mle <- candidate.log.like
      est.param <- candidate.param
    }
  }
  return(list(mle=local.mle,param=est.param))
}

nmf.poisson.vb <- function(data, learning.seed=1, H, prior.alpha=1, prior.beta=1, prior.gamma=1, prior.delta=1,
                           init.shape=8, iteration=100, seed.num = 50){
  ###alpha,beta,gamma,delta is hyperparameter of prior dist
  N <- dim(data)[1]
  M <- dim(data)[2]
  n <- dim(data)[3]
  
  set.seed(learning.seed)
  
  est.energy <- Inf
  est.param <- list()
  
  ##change initial parameters seed.num times
  for(seed.ind in 1:seed.num){
    ##initialize hyperparameter of variational posterior
    est.alpha <- matrix(rgamma(N*H,shape = init.shape),nrow=N,ncol=H)
    est.beta <- matrix(rgamma(N*H,shape = init.shape),nrow=N,ncol=H)
    est.gamma <- matrix(rgamma(H*M,shape = init.shape),nrow=H,ncol=M)
    est.delta <- matrix(rgamma(H*M,shape = init.shape),nrow=H,ncol=M)
    
    s <- array(0,dim=c(N,H,M,n))
    
    ##EM algorithm
    for(ite in 1:iteration){
      #E-step
      factor.ratio <- array(0,dim=c(N,H,M))
      for(nu in 1:N){
        for(tau in 1:M){
          L <- digamma(est.alpha[nu,])-log(est.beta[nu,])+digamma(est.gamma[,tau])-log(est.delta[,tau])
          L <- L - rep(1,H) * max(L)
          factor.ratio[nu,,tau] <- exp(L)/sum(exp(L))                    
          s[nu,,tau,] <- factor.ratio[nu,,tau]%*%t(data[nu,tau,])
        }
      }
      
      ##update of variational posterior of parameter A
      est.alpha <- apply(s,c(1,2),sum)+prior.alpha
      est.beta <- n*rep(1,N)%*%t(apply(est.gamma/est.delta,1,sum)+prior.beta)
      
      ##update of variational posterior of parameter B
      est.gamma <- apply(s,c(2,3),sum)+prior.gamma
      est.delta <- n*apply(est.alpha/est.beta,2,sum)%*%t(rep(1,M))
      
      ##calc mean param
      est.mean.param <- list(A=est.alpha/est.beta, B=est.gamma/est.delta)
      
      ##energy calculation
      dquasi.pois.apply <- function(x) sum(dquasi.pois(x,lambda1=est.mean.param$A%*%est.mean.param$B,lambda2=apply(factor.ratio,c(1,3),sum),log=T))
      current.energy <- sum((est.alpha-prior.alpha)*(digamma(est.alpha)-log(est.beta)) - (est.beta-prior.beta)*est.alpha/est.beta -
        lgamma(est.alpha) + lgamma(prior.alpha) + est.alpha*log(est.beta) - prior.alpha*log(prior.beta))
      current.energy <- current.energy + sum((est.gamma-prior.gamma)*(digamma(est.gamma)-log(est.gamma)) - (est.delta-prior.delta)*est.gamma/est.delta -
                                               lgamma(est.gamma) + lgamma(prior.gamma) + est.gamma*log(est.delta) - prior.gamma*log(prior.delta))
      current.energy <- current.energy -(sum(apply(data,3,dquasi.pois.apply)))
      
      # print(current.energy)
    }
    #calc likelihood
    candidate.param <- est.mean.param
    candidate.energy <- current.energy
    
    if(candidate.energy<est.energy){
      print(current.energy)
      est.energy <- candidate.energy
      est.param <- candidate.param
    }
  }
  return(list(energy=est.energy,param=est.param))  
}

dquasi.pois <- function(x,lambda1,lambda2,log=T){
  exp.arg <- -lambda1+x*log(lambda2)-lgamma(x+1)
  if(log==T){
    return(exp(exp.arg))
  }
  return(exp.arg)
}

dnmflikelihood <- function(data,param,log=T){
  if(log==T){
    dpois.apply <- function(x) sum(dpois(x,lambda=param$A%*%param$B,log=T))
    return(sum(apply(data,3,dpois.apply)))
  }else{
    dpois.apply <- function(x) prod(dpois(x,lambda=param$A%*%param$B, log=F))
    return(prod(apply(data,3,dpoi)))
  }
}