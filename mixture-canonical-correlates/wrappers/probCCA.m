function [cca,co] = probCCA (x1,x2,k)
% Standalone function for Probabilistic CCA [1]
% FORMAT [cca,co] = probCCA (x1,x2,k)
%
% Inputs:
% x1        [d1 x N] data matrix
% x2        [d2 x N] data matrix
% k         maximum latent dimension
%
% Outputs:
% cca   .W1{1} Factor loading matrix for x1   
%       .W2{1} Factor loading matrix for x2
%       .mu1{1} mean for x1
%       .mu2{1} mean for x2
%       .psi1{1} obs noise precision for x1
%       .psi2{1} obs noise precision for x2
%       .C1{1} obs noise covariance for x1
%       .C2{1} obs noise covariance for x2
% co    canonical correlations
%
% [1] Bach and Jordan (2006) A Probabilistic Interpretation of
% Canonical Correlation Analysis. Tech Rep 688, Dept Stats, Berkeley

cca.mu1{1} = mean(x1')';
cca.mu2{1} = mean(x2')';

[A,B,co] =  canoncorr(x1', x2');
codiag = diag(sqrt(co(1:k)));

xCov = cov(x1');
yCov = cov(x2');

cca.W1{1} = xCov * A(:,1:k) * codiag; 
cca.W2{1} = yCov * B(:,1:k) * codiag; 

cca.C1{1} = xCov - cca.W1{1}*cca.W1{1}';
cca.C2{1} = yCov - cca.W2{1}*cca.W2{1}';

cca.psi1{1} = inv(cca.C1{1});
cca.psi2{1} = inv(cca.C2{1});

cca.pi = 1;
cca.M = 1;