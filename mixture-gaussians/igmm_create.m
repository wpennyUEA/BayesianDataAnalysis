function [mix] = igmm_create (mix,x)
% Create new GMM component
% FORMAT [mix] = igmm_create (mix,x)
%
% mix       mixture model data structure
% x         [D x 1] data vector
%
% Pinto & Engel (2015) A Fast Incremental Gaussian Mixture 
% Model. PLoS One, Oct7.

mix.M = mix.M+1;
M = mix.M;

mix.state(M).m = x;
mix.state(M).Lambda = mix.lambda0*eye(mix.D);
mix.state(M).detC = mix.det0;

mix.sp(M) = 1;
mix.v(M) = 1;
mix.prior(M) = 1/sum(mix.sp);


