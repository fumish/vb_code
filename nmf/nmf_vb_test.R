source("nmf_functions.R")

n <- 50 #number of samples
N <- 5 
M <- 5
H <- 1 #true inner dim

##generate true parameters
true.seed <- 1
set.seed(true.seed)
Tshape <- 1
Tscale <- 1
true.param <- NMF_generate_param(M,N,H,shape=Tshape,scale=Tscale)

#generate data
sample.seed <- 2
set.seed(sample.seed)
x <- rNMF(n,true.param)

## learner inner dim
learning_H <- 3
learning.seed <- 2
seed.num <- 1

result <- nmf.poisson.vb(data=x,learning.seed=learning.seed, H=learning_H, seed.num=seed.num)
energy <- result$energy
est.param <- result$param
print(energy)
disp.param(true.param)
disp.param(est.param)
