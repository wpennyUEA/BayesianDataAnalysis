%
% ROBUST GENERAL LINEAR MODELS (RGLMs)
%
% This model can automatically remove spiky artifacts from a time series.
%
% -------------------------------------------------------------------------
% spm_rglm.m            Variational Bayes for Robust GLM [1]
% spm_glm.m             Special case of spm_rglm.m with single error
%                       variance level
% spm_boxcars.m         Create boxcar-like variable (fMRI block design)
% demo_rglm.m           Demo
%
% -------------------------------------------------------------------------
% REFERENCES
%
% [1] W.D. Penny, J. Kilner and F. Blankenburg (2007) Robust Bayesian General Linear
% Models. Neuroimage 36(3): 661-671.
