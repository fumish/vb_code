rdirichlet <- function(phi){
  gamma_val <- rgamma(length(phi),shape = phi)
  return(gamma_val/sum(gamma_val))
}

# logdet <- function(X){
#   return( 2*sum(log(diag(chol(X)))) )
# }
# 
# mvdnorm <- function(x,mu,Sigma,log=F){
#   exp_term <- t(x-mu)%*%solve(Sigma)%*%(x-mu)/2-dim(Sigma)[1]/2*log(2*pi)-determinant(Sigma,logarithm=T)$modulus/2
#   if(log == T){
#     return(exp_term)
#   }else{
#     return(exp(exp_term))
#   }
# }

gmm.loglikelihood <- function(data,param){
  library(mvtnorm)
  n <- dim(data)[1]
  K <- dim(param$mu)[2]
  L <- matrix(1,nrow=n,ncol=1)%*%log(param$ratio)
  for(k in 1:K){
    L[,k] <- L[,k] + dmvnorm(data,param$mu[,k],param$Sigma[,,k],log=T)
  }
  maxL <- apply(L,1,max)
  L <- L - maxL%*%matrix(1,nrow=1,ncol=K)

  val <- sum(log(apply(exp(L),1,sum))+maxL)
  
}

GMM_generate_param <- function(K0,phi,beta,Cov_df,Cov_Sigma){
  library(mvtnorm)
  M <- dim(Cov_Sigma)[1]
  
  zero_vec <- rep(0,times=M)
  
  ratio <- rdirichlet(phi)
  
  Sigma <- array(0,dim=c(M,M,K0))
  mu <- matrix(0,nrow=M,ncol=K0)
  
  for(k in 1:K0){
    Sigma[,,k] <- solve(rWishart(1,Cov_df,Cov_Sigma)[,,1])
    mu[,k] <- rmvnorm(1,mean=zero_vec, sigma=Sigma[,,k]/beta)
  }
  
  return(list(ratio=ratio, mu=mu, Sigma=Sigma))
}

rgmm <- function(n,param){
  library(mvtnorm)
  M <- dim(param$mu)[1] #dim of data
  K <- dim(param$mu)[2] #num of components
  
  data <- matrix(0, nrow=n, ncol=M)
  
  latent_val <- t(rmultinom(n,1,prob = param$ratio))
  label <- numeric(n)
  
  for(i in 1:n){
    label[i] <- which(latent_val[i,]==1)
    data[i,] <- rmvnorm(1,mean=param$mu[,label[i]], sigma=param$Sigma[,,label[i]])
  }

  return(list(data=data, label=label))
  
}

gmm.cov.em <- function(data,K, learning_num=100, seed=1, seed_num=10){
  library(mvtnorm)
  
  set.seed(seed)
  
  n <- dim(data)[1]
  M <- dim(data)[2]
  
  mle <- -Inf
  param <- list()
  
  for(seed_iteration in 1:seed_num){
    ##initialize
    latent_prob <- matrix(0,nrow=n,ncol=K)
    for(i in 1:n){
      latent_prob[i,] <- rdirichlet(rep(1,K))
    }
    
    Sigma <- array(0,dim=c(M,M,K))
    for(learning_iteration in 1:learning_num){
      #M-step
      n_k <- apply(latent_prob,2,sum)
      ratio <- n_k/n
      mu <- (t(data) %*% latent_prob) / (matrix(1,nrow=M,ncol=1) %*% n_k)
      for(k in 1:K){
        for(j in 1:M){
          for(l in 1:M){
            Sigma[j,l,k] <- t(latent_prob[,k]*(data[,j]-mu[j,k])) %*% (data[,l]-mu[l,k]) / n_k[k]
          }
        }
      }
      
      #E-step
      L <- matrix(1,nrow=n,ncol=1)%*%log(ratio)
      for(k in 1:K){
        L[,k] <- L[,k] + dmvnorm(data,mu[,k],Sigma[,,k],log=T)
      }
      maxL <- apply(L,1,max)
      L <- L - maxL%*%matrix(1,nrow=1,ncol=K)
      latent_prob <- exp(L) / (apply(exp(L),1,sum)%*%matrix(1,nrow=1,ncol=K))
      
      candidate_mle <- sum(log(apply(exp(L),1,sum))+maxL)
      
      if(is.nan(candidate_mle)){
        break
      }
      
      previous_mle <- candidate_mle
      previous_param <- list(ratio=ratio,mu=mu,Sigma=Sigma)
    }
    # print(previous_mle)
    if(previous_mle > mle){
      mle <- previous_mle
      param <- previous_param
    }
    
  }
  if(length(param)==0){
    print("miss to estimate")
  }
  return(list(param=param,mle=mle))
}

gmm.cov.vb <- function(data, K, phi=10, beta=0.001, nu.a = 1, nu.b = 1, learning_num=100, seed=1, seed_num=10){
  library(mvtnorm)
  
  log.det <- function(mat){
    return(unlist(determinant(mat),use.names = F)[1])
  }
  
  set.seed(seed)
  
  n <- dim(data)[1]
  M <- dim(data)[2]
  
  energy <- Inf
  param <- list()
  
  if(nu.a <= M-1){
    print("error for hyperparameter at the degrees of freedom of wishart distribution!")
    print("note: degrees of freedom nu.a is larger than M-1")
    return(NULL)
  }
  nu.b <- array(diag(rep(nu.b,M)),dim=c(M,M,K))
  
  for(seed_iteration in 1:seed_num){
    ##initialize
    latent_prob <- matrix(0,nrow=n,ncol=K)
    for(i in 1:n){
      latent_prob[i,] <- rdirichlet(rep(1,K))
    }
    
    ##learning part
    for(learning_iteration in 1:learning_num){
      #M-step
      n_k <- apply(latent_prob,2,sum)
      est_phi <- n_k+phi
      est_beta <- n_k+beta
      est_mean <- t(data)%*%latent_prob/(matrix(1,nrow=M,ncol=1)%*%est_beta)
      est_nu.a <- n_k + nu.a
      est_nu.b <- array(0,dim=c(M,M,K))
      for(k in 1:K){
        inv_nu.b <- solve(nu.b[,,k])
        inv_est_nu.b <- matrix(0,nrow=M,ncol=M)
        for(j in 1:M){
          for(l in 1:M){
            inv_est_nu.b[j,l] <- sum((latent_prob[,k]*(data[,j]-est_mean[j,k]))*(data[,l]-est_mean[l,k])) + beta*est_mean[j,k]*est_mean[l,k]+inv_nu.b[j,l]
          }
        }
        est_nu.b[,,k] <- solve(inv_est_nu.b)
      }
      
      L <- matrix(1,nrow=n,ncol=1)%*%(digamma(est_phi)-digamma(sum(est_phi)))
      for(k in 1:K){
        L[,k] <- L[,k] + dmvnorm(data,est_mean[,k],sigma = solve(est_nu.a[k]*est_nu.b[,,k]),log=T) 
          +M/2*log(2)-M*log(est_nu.a[k])/2
          + sum(digamma((est_nu.a[k]+(1-1:M))/2)/2)
      }
      maxL <- apply(L,1,max)
      L <- L - maxL%*%matrix(1,nrow=1,ncol=K)
      latent_prob <- exp(L) / (apply(exp(L),1,sum)%*%matrix(1,nrow=1,ncol=K))
      
      # candidate_energy <- sum((est_phi-phi)*(digamma(est_phi)-digamma(sum(est_phi)))-lgamma(est_phi))+lgamma(sum(est_phi))-lgamma(K*phi)+K*lgamma(phi)
      # candidate_energy <- candidate_energy + sum(log(est_beta/beta)/2+M*beta/(2*est_beta)-M/2)
      # for(k in 1:K){
      #   candidate_energy <- candidate_energy + beta*t(est_mean[,k])%*%(est_nu.a[k]*est_nu.b[,,k])%*%est_mean[,k]/2
      #   candidate_energy <- candidate_energy + sum((est_nu.a[k]-nu.a)*digamma((est_nu.a[k]+(1-1:M))/2)/2+lgamma((nu.a+(1-1:M))/2)-lgamma((est_nu.a[k]+(1-1:M))/2))
      #   candidate_energy <- candidate_energy + est_nu.a[k]*(sum(diag(solve(nu.b[,,k])%*%est_nu.b[,,k]))-M)/2 - nu.a*(log.det(est_nu.b[,,k])+log.det(solve(nu.b[,,k])))/2
      #   # candidate_energy <- candidate_energy + est_nu.a[k]*(sum(diag(solve(nu.b[,,k])%*%est_nu.b[,,k]))-M)/2 -nu.a*(determinant(est_nu.b[,,k])$modulus+determinant(solve(nu.b[,,k]))$modulus)/2
      # }
      # candidate_energy <- candidate_energy - sum(log(apply(exp(L),1,sum))+maxL)
      # # calc mean parameters
      # Sigma <- array(0,dim=c(M,M,K))
      # for(k in 1:K){
      #   Sigma[,,k] <-solve(est_nu.a[k]*est_nu.b[,,k])
      # }
      # candidate_param <- list(ratio=est_phi/sum(est_phi),mu=est_mean,Sigma=Sigma)      
      # print(candidate_energy)
    }
    
    #calc energy
    candidate_energy <- sum((est_phi-phi)*(digamma(est_phi)-digamma(sum(est_phi)))-lgamma(est_phi))+lgamma(sum(est_phi))-lgamma(K*phi)+K*lgamma(phi)
    candidate_energy <- candidate_energy + sum(log(est_beta/beta)/2+M*beta/(2*est_beta)-M/2)
    for(k in 1:K){
      candidate_energy <- candidate_energy + beta*t(est_mean[,k])%*%(est_nu.a[k]*est_nu.b[,,k])%*%est_mean[,k]/2
      candidate_energy <- candidate_energy + sum((est_nu.a[k]-nu.a)*digamma((est_nu.a[k]+(1-1:M))/2)/2+lgamma((nu.a+(1-1:M))/2)-lgamma((est_nu.a[k]+(1-1:M))/2))
      candidate_energy <- candidate_energy + est_nu.a[k]*(sum(diag(solve(nu.b[,,k])%*%est_nu.b[,,k]))-M)/2 - nu.a*(log.det(est_nu.b[,,k])+log.det(solve(nu.b[,,k])))/2
      # candidate_energy <- candidate_energy + est_nu.a[k]*(sum(diag(solve(nu.b[,,k])%*%est_nu.b[,,k]))-M)/2 -nu.a*(determinant(est_nu.b[,,k])$modulus+determinant(solve(nu.b[,,k]))$modulus)/2
    }
    candidate_energy <- candidate_energy - sum(log(apply(exp(L),1,sum))+maxL)
    # calc mean parameters
    Sigma <- array(0,dim=c(M,M,K))
    for(k in 1:K){
      Sigma[,,k] <-solve(est_nu.a[k]*est_nu.b[,,k])
    }
    candidate_param <- list(ratio=est_phi/sum(est_phi),mu=est_mean,Sigma=Sigma)
    
    # print(candidate_energy)
    ##evaluation among obtained parameters
    if(candidate_energy < energy){
      energy <- candidate_energy
      param <- candidate_param
      # print(candidate_energy)
      # print(candidate_param)
    }
    
  }
  return(list(param=param,energy=energy))
}

gmm.em <- function(data,K, learning_Sigma, learning_num=100, seed=1, seed_num=10){
  library(mvtnorm)
  
  set.seed(seed)
  
  n <- dim(data)[1]
  M <- dim(data)[2]
  
  mle <- -Inf
  param <- list()
  
  for(seed_iteration in 1:seed_num){
    ##initialize
    latent_prob <- matrix(0,nrow=n,ncol=K)
    for(i in 1:n){
      latent_prob[i,] <- rdirichlet(rep(1,K))
    }
    
    for(learning_iteration in 1:learning_num){
      #M-step
      n_k <- apply(latent_prob,2,sum)
      ratio <- n_k/n
      mu <- (t(data) %*% latent_prob) / (matrix(1,nrow=M,ncol=1) %*% n_k)
      
      #E-step
      L <- matrix(1,nrow=n,ncol=1)%*%log(ratio)
      for(k in 1:K){
        L[,k] <- L[,k] + dmvnorm(data,mu[,k],diag(rep(learning_Sigma[k],M)),log=T)
      }
      maxL <- apply(L,1,max)
      L <- L - maxL%*%matrix(1,nrow=1,ncol=K)
      latent_prob <- exp(L) / (apply(exp(L),1,sum)%*%matrix(1,nrow=1,ncol=K))
      
      candidate_mle <- sum(log(apply(exp(L),1,sum))+maxL)
      
      if(is.nan(candidate_mle)){
        break
      }
      
      Sigma <- array(0,dim=c(M,M,K))
      for(k in 1:K){
        Sigma[,,k] <-diag(rep(learning_Sigma[k],M)) 
      }      
      
      previous_mle <- candidate_mle
      previous_param <- list(ratio=ratio,mu=mu,Sigma=Sigma)
    }
    # print(previous_mle)
    if(previous_mle > mle){
      mle <- previous_mle
      param <- previous_param
    }
    
  }
  return(list(param=param,mle=mle))
}

gmm.vb <- function(data,K, learning_Sigma, phi=5, beta=0.001, learning_num=100, seed=1, seed_num=10){
  library(mvtnorm)
  
  set.seed(seed)
  
  n <- dim(data)[1]
  M <- dim(data)[2]
  
  energy <- Inf
  param <- list()
  
  for(seed_iteration in 1:seed_num){
    ##initialize
    latent_prob <- matrix(0,nrow=n,ncol=K)
    for(i in 1:n){
      latent_prob[i,] <- rdirichlet(rep(1,K))
    }
    
    for(learning_iteration in 1:learning_num){
      #M-step
      est_phi <- apply(latent_prob,2,sum)+phi
      est_beta <- apply(latent_prob,2,sum)/learning_Sigma+beta
      est_mean <- t(data)%*%latent_prob/(matrix(1,nrow=M,ncol=1)%*%(learning_Sigma*est_beta))
      
      # print(est_mean)
      
      #E-step
      # log_N_nume <- -diag(data%*%t(data))%*%matrix(1,nrow=1,ncol=K)/2
      #   +data%*%est_mean
      #   -matrix(1,nrow=n,ncol=1)%*%diag(t(est_mean)%*%est_mean)/2
      
      L <- matrix(1,nrow=n,ncol=1)%*%(digamma(est_phi)-digamma(sum(est_phi))-M/(2*learning_Sigma*est_beta))
      for(k in 1:K){
        L[,k] <- L[,k] + dmvnorm(data,est_mean[,k],diag(rep(learning_Sigma[k],M)),log=T)
      }
      maxL <- apply(L,1,max)
      L <- L - maxL%*%matrix(1,nrow=1,ncol=K)
      latent_prob <- exp(L) / (apply(exp(L),1,sum)%*%matrix(1,nrow=1,ncol=K))
      
      candidate_energy <- sum((est_phi-phi)*(digamma(est_phi)-digamma(sum(est_phi)))-lgamma(est_phi))+lgamma(sum(est_phi))-lgamma(K*phi)+K*lgamma(phi) 
      candidate_energy <- candidate_energy + sum(M*log(est_beta/beta)/2+M*beta/(2*est_beta)+beta*diag(t(est_mean)%*%est_mean)/2-M/2)
      candidate_energy <- candidate_energy - sum(log(apply(exp(L),1,sum))+maxL)
      
      # print(candidate_energy)
      
      Sigma <- array(0,dim=c(M,M,K))
      for(k in 1:K){
        Sigma[,,k] <-diag(rep(learning_Sigma[k],M)) 
      }
      candidate_param <- list(ratio=est_phi/sum(est_phi),mu=est_mean,Sigma=Sigma)
    }
    
    # print(candidate_energy)
    if(candidate_energy < energy){
      energy <- candidate_energy
      param <- candidate_param
      # print(candidate_energy)
      # print(candidate_param)
    }
    
  }
  return(list(param=param,energy=energy))
}

gmm_nondirichlet_lambda <- function(M, K0, K){
  lambda <- (M*K0+K-1)/2
  m <- 1
  return(list(lambda=lambda, m=m))
}

gmm_dirichlet_lambda <- function(M,K0,K,phi=(M+1)/2){
  quasi_d <- c(K0*(M+1)+(K-K0-1)*phi,K*M+K0)

  lambda <- min(quasi_d)/2
  m <- 1
  
  return(list(lambda=lambda, m=m))
}

gmm_dirichlet_lambda2 <- function(M,K0,K,phi=(M+1)/2){
  lambda <- (K0*M+K0-1+(K-K0)*phi)
  m <- 1
  
  return(list(lambda=lambda, m=m))
}

sample_var <- function(x) return(var(x)*(length(x)-1)/length(x))
