%
% MULTIVARIATE LINEAR MODELLING (MLM)
%
% i.e. multiple independent and multiple dependent variables.
% a.k.a. Bayesian Ridge Regression
%
% -------------------------------------------------------------------------
% spm_mlm_bayes.m           Variational Bayes for MLM
% spm_mlm_posthoc.m         Compute Bayes Factors using Savage-Dickey
% spm_mlm_makecon.m         Make contrasts for spm_mlm_posthoc.m
% demo_mlm_basic.m          MLM Demo
% demo_mlm_posthoc.m        Demo of posthoc inference
% demo_mlm_priors.m         Show effect of priors
%
% -------------------------------------------------------------------------
% spm_cva_prob.m            Probabilistic Canonical Variates Analysis [2]
% spm_cva_compare.m         Get AIC, BIC for allowable latent dimensions
% demo_cva_params.m         Compare true and estimated parameters
% demo_cva_order.m          Compare LogLikelihood, AIC and BIC
% 
% -------------------------------------------------------------------------
%
% spm_mlm_bayes.m is the same as spm_mar.m [1] but where the independent
% variables are specified by the user (in spm_mar.m they are lagged time
% series vectors).
%
% -------------------------------------------------------------------------
% REFERENCES
% 
% [1] W.D. Penny and S.J. Roberts. Bayesian Multivariate Autoregresive Models 
% with structured priors. IEE Proceedings on Vision, Image and Signal Processing, 
% 149(1):33-41, 2002
%
% [2] F. Bach and M. Jordan (2005) A probabilistic interpretation of canonical
% correlation analysis. Dept. Stats, Univ California, Berkeley CA. 
% Tech Rep 688.
%
