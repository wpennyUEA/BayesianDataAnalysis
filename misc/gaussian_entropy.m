function [H] = gaussian_entropy (C)
% Compute entropy of multivariate Gaussian density
% FORMAT [H] = gaussian_entropy (C)
%
% C     Covariance matrix
% 
% H     Entropy (with natural logs)

LD = spm_logdet(C);
d = size(C,1);
H = 0.5*d*(1+log(2*pi))+0.5*LD;