function [p] = vbcca_numparams (cca)
% Return number of parameters of CCA model (excluding obs noise matrices)
% FORMAT [p] = vbcca_numparams (cca)
%
% cca       model structure returned by vbcca.m

M = cca.M;

% Priors
p = M;

% Factor matrices
[d1,q] = size(cca.W1{1});
p = p+d1*q*M;
[d2,q] = size(cca.W2(1));
p = p+d2*q*M;

% Means
p = p + d1*M + d2*M;

