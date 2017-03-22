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
NMF assumes that a given non-negative matrix can be decomposed by two non-negative matrices. 
NMF_functions.R includes the functions to work the vb algorithm.
NMF_em_test.R is the script to implement Expectation-Maximization algorithm, which searches a local mixima of the likelihood.
NMF_vb_test.R is the script to implement vb algorithm.
3. Latent Dirichlet Allocation(LDA)  
LDA is in the folder lda.
LDA assumes that there is some document and each document has different mixing ratio to generate their hidden topics,
while same topic generate observal data by same probability. Thus parameter of the model is mixing ratio for each document and parameter of each component.
In this folder, I assume that each component obeys Totally Asymmetric Simple Exclusion Process(TASEP) and Zero Range Process(ZRP), which are the model in the field of traffic flow.
VB algorithm for LDA with TASEP is in tasep_artificial_data and the ZRP one is in zrp_artificial_data.
