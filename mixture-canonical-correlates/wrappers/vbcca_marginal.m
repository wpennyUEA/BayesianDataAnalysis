function [p1,p2] = vbcca_marginal (cca,X1p,X2p)
% Compute marginal densities 
% FORMAT [p1,p2] = vbcca_marginal (cca,X1p,X2p)
%
% cca   data structure returned by vbcca.m
% X1p   [d1 x N] matrix of X1 data points over which to compute density
% X2p   [d1 x N] matrix of X2 data points over which to compute density
%
% p1    [1 x N] vector of probability densities for X1p
% p2    [1 x N] vector of probability densities for X2p

p=[];

for k=1:cca.M,
    
    mu = cca.mu1{k};
    C = cca.W1{k}*cca.W1{k}' + diag(diag(inv(cca.psi1{k})));  % Enforce diag obs noise cov
    pdf1(k,:) = spm_mvNpdf (X1p,mu,C);
    
    mu = cca.mu2{k};
    C = cca.W2{k}*cca.W2{k}' + diag(diag(inv(cca.psi2{k})));
    pdf2(k,:) = spm_mvNpdf (X2p,mu,C);
end

% Prior/Frequency
pr = cca.pi';

p1 = pr*pdf1;
p2 = pr*pdf2;
