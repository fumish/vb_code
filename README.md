# vb_code
This repository includes scripts for variational bayes method, which is a method
to obtain approximated posterior distribution for bayes posterior distribution
I put the following vb algorithms and implementation for the algorithms:
1. Guassian Mixture Model(GMM)  
GMM is in the folder gmm.
GMM is mixture of gaussian distribution.
As the vb algorithm, I don't only assume multi-dimensional mean estimation, but also covariance of the elements. 
gmm.R includes the functions to work the vb algorithm.
gmm_cov_test.R is the script to implement them.
2. Non-negative Factorization(NMF)  
NMF is in the folder nmf.
NMF assume that a given non-negative matrix can be decomposed by two non-negative matrices. 
NMF_functions.R includes the functions to work the vb algorithm.
NMF_em_test.R is the script to implement Expectation-Maximization algorithm, which searches a local mixima of the likelihood.
NMF_vb_test.R is the script to implement vb algorithm.
