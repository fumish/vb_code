sample_seed <- 1
set.seed(sample_seed)

source("sbic_functions.R")
source("gmm.R")

#problem setting
M <- 2
n <- 1000 #num of training samples

#true setting
K0 <- 3
phi <- rep(5,K0)
df <- 4
Sigma <- diag(M)
beta <- 0.01
true_param <- GMM_generate_param(K0,phi,beta,df,Sigma)
true_info <- rgmm(n,true_param)

# true_sigma <- rep(1,M)
# true_param$Sigma <- array(diag(true_sigma),dim=c(M,M,K0))

#learner setting
K <- 3
dloglikelihood <- gmm.loglikelihood
learning_num <- 100
seed_num <- 10
learning_seed <- 1
phi <- 10
beta <- 0.001
nu.a <- 4
nu.b <- 5

data <- true_info$data

result <- gmm.cov.em(data,K, learning_num=learning_num, seed=learning_seed, seed_num=seed_num)

# result <- gmm.cov.vb(data, K, phi=phi, beta = beta, nu.a=nu.a, nu.b=nu.b, learning_num = learning_num, seed=learning_seed,seed_num=seed_num)
