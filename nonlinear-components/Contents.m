%
% NEIGHBOURHOOD COMPONENT ANALYSIS (NCA)
%
% NCA finds a linear transformation of data that maximises classification
% performance under a k-Nearest Neighbour model. The transformation
% can be thought of as finding the optimal (linear) distance measure for KNN.
% 
% It can find a low-dimensional projection of data that uses class labels
% (unlike PCA which finds data subspaces but not necessarily ones useful
% for classification). Unlike Linear Discriminant Analysis, the
% classifications it makes (being based on kNN) are nonlinear.
%
% This code includes Maximum Likelihood, Bayesian Regularisation and Bayesian
% Pruning options [2].
%
% The implementation here allows for a full linear transformation e.g. arbitrary linear 
% combination of original data. Mathworks' implementation (fscnca.m) assumes 
% a diagonal approximation, giving each feature a "score" [3]. 
%
% -------------------------------------------------------------------------
% nca_org.m                     Bayesian NCA
% nca_org_cost_grad.m           NCA Cost function and gradient
% tune_alpha.m                  NCA Cost as a function of step size
% 
% nca_prune.m                   NCA with Bayesian pruning of features
% nca_log_prior.m               Log of NCA prior
% nca_org_cost_row.m            Row-specific contribution to cost
% nca_remove_row.m              Remove row (feature)
%
% demo_nca_synth.m              Demo of ML, Bayes Reg and Pruning options
% demo_nca_pima.m               Demo on classic data set
%
% -------------------------------------------------------------------------
% REFERENCES
% 
% [1] Goldberger et al. Neighborhood Component Analysis, NIPS 2004
%
% [2] W Penny (2020) Sparse Neighborhood Component Analysis. TechReport, UEA.
%
% [3] Yang, W et al. Neighborhood Component Feature Selection for High-Dimensional Data. 
% Journal of Computers. Vol. 7, Number 1, January, 2012.
