function [mix] = igmm_init (x0,s0)
% Initialise Gaussian mixture model ready for online learning
% FORMAT [mix] = igmm_init (x0,s0)
%
% x0        [D x 1] vector corresponding to first data point
% s0        [1 x 1] initial standard deviation of component
%
% mix       mixture model data structure
%
% Pinto & Engel (2015) A Fast Incremental Gaussian Mixture 
% Model. PLoS One.

mix.D = length(x0); % Input dimension
mix.M=0;            % Number of components
mix.beta=0.001;       % Percentile of Chi^2 density

sig02=s0^2;
mix.lambda0=1/(sig02);  % inverse sigma_ini
mix.det0 = sig02^mix.D;

% Chi^2 Threshold for creating new component
mix.chi2 = spm_invXcdf(1-mix.beta,mix.D); 

mix = igmm_create(mix,x0);