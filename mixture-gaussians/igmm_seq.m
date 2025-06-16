function [mix] = igmm_seq (x,s0,diag)
% Online learning of Gaussian mixture mode
% FORMAT [mix] = igmm_seq (x,sx,diag)
%
% x         [D x N] time series of N data points
% s0        [1 x 1] initial standard deviation of component
% diag      1 for diagonal covariances (default = 0)
%
% mix       Mixture model data structure
%
% Pinto & Engel (2015) A Fast Incremental Gaussian Mixture 
% Model. PLoS One.

if nargin < 3 | isempty(diag), diag=0; end

[D,N] = size(x);

mix = igmm_init (x(:,1),s0);

if diag
    for n=2:N,
        mix = igmm_update_diag (mix,x(:,n));
    end
else
    for n=2:N,
        mix = igmm_update (mix,x(:,n));
    end
end