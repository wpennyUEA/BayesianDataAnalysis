function [y] = mci_softmax (x,lambda)
% Softmax function
% FORMAT [y] = mci_softmax (x,lambda)
% 
% x         vector 
% lambda    inverse temperature
%
% y         softmax output

if nargin < 2 | isempty(lambda)
    lambda=1;
end

y=exp(lambda*x)+eps;
y=y./sum(y);

