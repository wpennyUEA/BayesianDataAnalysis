%
% AUTOREGRESSIVE (AR) MODELLING
%
% Predict the next value in single/multiple time series as a linear
% combination of the P previous samples (AR-P, MAR-P models).
% Make inferences about which time series are related to each other.
% The AR coefficients can be used to estimate power spectra (see /spectral)
%
% -------------------------------------------------------------------------
% spm_ar.m              AR modelling for single time series
% demo_ar.m             Demo
% spm_ar_pred.m         Predict next time series value
%
% -------------------------------------------------------------------------
% spm_rar.m             As spm_ar.m but with mixture noise process
% demo_rar.m            Demo
%
% -------------------------------------------------------------------------
% spm_mar.m             AR modelling for multiple time series
% demo_mar.m            Demo
% spm_mar_conn.m        Chi^2 test for significance of connections
% spm_mar_gen.m         Generate data from known model
% spm_mar_pred.m        Predict next time series values
% spm_mar_prior.m       Set up global, lag, interaction etc. priors
% test_spm_mar.m        Time fitting of large model
%
% -------------------------------------------------------------------------
% spm_get_omega.m       Get expected error of MAR model
% spm_kl_eig_normal.m   Kullback Liebler divergence of two normal densities
% 
% -------------------------------------------------------------------------
% REFERENCES
% 
% [1] W.D. Penny and S.J. Roberts. Bayesian Multivariate Autoregresive Models 
% with structured priors. IEE Proceedings on Vision, Image and Signal Processing, 
% 149(1):33-41, 2002
%
% [2] L. Harrison, W.D. Penny, and K.J. Friston. Multivariate Autoregressive 
% Modelling of fMRI time series. NeuroImage, 19(4):1477-1491, 2003
% 
% [3] W. Penny and L. Harrison. Multivariate autoregressive models. In K. Friston, 
% J. Ashburner, S. Kiebel, T. Nichols, and W. Penny, editors, Statistical 
% Parametric Mapping: The analysis of functional brain images. Elsevier, 
% London, 2006

