%
% GAUSSIAN MIXTURE MODELS (GMMs)
%
% Automatic clustering of multivariate data allowing number of clusters,
% their means and covariances to be identified. Offline and online 
% algorithms.
% 
% --------------------------------------------------------------------------------------
% spm_mix.m                 Variational Bayes (VB) multivariate mixture modelling [1,2]
% demo_mix_1d.m             Demo
% demo_mix_2d.m             Demo
%
% --------------------------------------------------------------------------------------
% spm_kmeans.m              Multivariate K-means clustering 
% spm_kmeans1.m             Univariate K-means clustering
% spm_MNpdf.m               PDF of Multivariate Normal
% spm_samp_gauss.m          Sample from Gausian density
% spm_samp_mix.m            Sample from GMM
% spm_mix_plot1d.m          Plot probability densities for 1D data
% spm_mix_plot2d.m          Plot probability contours for 2D data
%
% --------------------------------------------------------------------------------------
% igmm_init.m               Set up Incremental (online) GMM learning [3]
% igmm_create.m             Create new latent variable
% igmm_seq.m                Incremental learning (i.e. once through data set)
% igmm_update.m             Update means, (co)variances, proportions
% igmm_update_diag.m        Assuming diagonal covariances
% demo_igmm.m               Demo
% 
% -------------------------------------------------------------------------
% REFERENCES
%
% [1] H. Attias (2000) A Variational Bayesian framework for Graphical Models, 
%     NIPS 12, 209-215, MIT press, 2000.
%
% [2] W.D. Penny (2001) Variational Bayes for d-dimensional Gaussian mixture models 
%     Wellcome Department of Imaging Neuroscience, University College London.
%
% [3] Pinto & Engel (2015) A Fast Incremental Gaussian Mixture Model. PLoS One.

