source("nmf_functions.R")

n <- 100 #number of samples
N <- 5 
M <- 5
H <- 2 #true inner dim

##generate true parameters
true.seed <- 1
set.seed(true.seed)
Tshape <- 0.5
Tscale <- 1
true.param <- NMF_generate_param(M,N,H,shape=Tshape,scale=Tscale)


#generate data
sample.seed <- 2
set.seed(sample.seed)
x <- rNMF(n,true.param)
##true parameters
# TA <- matrix(c(0.5,1.5, 0,2, 1,2), nrow=N, ncol=H, byrow=T)
# TB <- matrix(c(1,0, 2,0, 3,7), nrow=H, ncol=M, byrow=T)

TA <- matrix(c(0.5,1.5, 0,2, 1,2), nrow=N, ncol=H, byrow=T)
TB <- matrix(c(1,0, 2,0, 3,7), nrow=H, ncol=M, byrow=T)

# TA <- matrix(0, nrow=N, ncol=H, byrow=T)
# TB <- matrix(0, nrow=H, ncol=M, byrow=T)

## sample is taken from poisson dist
x <- array(0, dim=c(N,M,n))
for(i in 1:n){
  for(j in 1:N){
    for(k in 1:M){
      x[j,k,i] <- rpois(1, (TA%*%TB)[j,k])
    }
  }
}

## learner inner dim
learning_H <- 1
learning.seed <- 2

result <- nmf.em(data=x,learning.seed=learning.seed, H=learning_H)
log.like <- result$mle
param <- result$param
print(log.like)
print(param$A%*%param$B)