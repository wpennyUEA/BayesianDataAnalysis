%
% PRINCIPAL COMPONENT ANALYSIS (PCA)
%
% These routines can automatically select an appropriate 
% number of principal components, shrinking irrelevant ones (variational)
% or using an explicit model selection criterion (Laplace)
%
% -------------------------------------------------------------------------
% spm_vpca.m            Variational Bayes (VB) for PCA [1]
% spm_vpca_f.m          Compute model evidence
% spm_vpca_init.m       Initialise
% spm_vpca_update.m     Update parameters
% 
% demo_vpca_small.m     Demo on low-dimensional data set
% demo_vpca_big.m       Demo on high-dimensional data set
% demo_svd_big.m        Reconstruct image using Singular Value
%                       Decomposition (SVD)
%
% spm_pca_order.m       Laplace approximation to model evidence for PCA [2,3]
% demo_pca_order.m      Demo
% 
% -------------------------------------------------------------------------
% REFERENCES
% 
% [1] C. Bishop. Variational Principal Components, ANN, 1999.
%
% [2] T.P. Minka. Automatic choice of dimensionality for PCA. Technical Report
% 514, MIT Media Lab, Perceptual Computing Section, 2000.
%
% [3] W. Penny, S. Roberts and R. Everson (2000) ICA: model order selection
% and dynamic source models. ICA: Principles and Practice, pages 299-314. 
% Cambridge University Press.