function [y] = mci_bernoulli_rnd (g,N)
% Generate binary random variable
% FORMAT [y] = mci_bernoulli_rnd (g,N)
%
% g         p(y=1)
% N         Number of variables to generate
%
% y         [N x 1] binary vector

y = spm_multrnd([g 1-g],N);
ind = find(y==2);
y(ind) = 0;